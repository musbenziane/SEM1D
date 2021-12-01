program SEM1D

    !###################################################################################################################
    ! Created by mus on 28/07/2021.
    ! Spectral Elements 1D Solver with a regular mesh for the elastic wave equation 
    ! NOTE: This program is yet to be tested against analytical solutions.
    ! 
    ! Language: Fortran 90, with parralel impelementation using OpenMP API
    ! 
    ! Sources:  Igel 2017
    ! 
    ! The code used for arbitrary GLL points and weights was created by M.I.T departement of engineering. Link is hereafter
    ! https://geodynamics.org/cig/doxygen/release/specfem3d/    | file name: gll_library.f90
    ! 
    ! This is part of the Numerical Modelling Workshop.
    ! 
    ! Supervisor: Pr. Emmanual Chaljub
    ! Author    : Mus Benziane
    !
    ! Input file: [Example]
    ! testing              ! Model name prefix (model names: prefix_vp, prefix_rho)
    ! 3                    ! Polynomial order
    ! 1000                 ! Number of elements
    ! 20                   ! Element size
    ! 5.                   ! Wavelet's peak frequency
    ! 0.001                ! Time step
    ! 10000                ! Number of time steps
    ! 1                    ! Source location (element number)
    ! 1                    ! Source location (gll point)
    ! 1                    ! [1/2/3/4] 1: Free surface, 2: Rigid wall, 3: Periodic, 4: Sponge 
    ! 1                    ! Boundary condition only on the left side [not for periodic BC | on the right side for absorbing]
    ! 20                   ! Sponge layer width in gll points
    ! .65                  ! Att constant for sponge layer
    ! 
    ! -> Model files in C-Style binary floats [doubles]: Vs, Rho files are needed.
    !                                                  : For simple models, use create1Dmodel_files.f90
    ! 
    ! -> Outputs are created in OUTPUT/ if OUTPUT/ is not created by the user, the program will not handle it.
    !    Output files in OUTPUT directory:
    !
    !####################################################################################################################
    
    !$ USE OMP_LIB
    implicit none

    interface connectivity_matrix
        function connectivity_matrix(N,ne)
            integer, intent(in)                                :: N, ne
            integer connectivity_matrix(N+1,ne)
        end function
    end interface


    real (kind=8)                                              :: Jc, Jci, h, f0, dt, sum, CFL, mindist, lambdamin, taper
    real (kind=8)                                              :: time, t_cpu_0, t_cpu_1, t_cpu, tmp, tmpBC, attConst, sd
                                                            
    real (kind=8), dimension(:), allocatable                   :: xi, wi, v1D, rho1D, Me, M, xgll, rho1Dgll, v1Dgll, g,tauL 
    real (kind=8), dimension(:), allocatable                   :: mu1Dgll, u, uold, unew, F, src, temp1, temp2, temp3 ,&
                                                                  udot, udotnew, uddot, uddotnew, sigma, recsigma
    real (kind=8), dimension(:,:), allocatable                 :: lprime, Minv, Kg, Ke,Uout, as, Udotout, sigmaout
    integer, dimension(:,:), allocatable                       :: Cij
    integer                                                    :: N, ne, ngll, i, j, k, l, nt, t, el, isnap, reclsnaps
    integer                                                    :: ir, t0, t1, esrc, gsrc, bc, sbc, gWidth, IC
    character(len=40)                                          :: filename, filecheck, outname, modnameprefix
    !$ integer                                                 :: n_workers
    logical                                                    :: OMPcheck = .false.


    write(*,*) "##########################################"
    write(*,*) "############### OpenMP     ###############"
    !$ OMPcheck = .true.
    if (OMPcheck) then
        !$OMP PARALLEL
        !$ n_workers = OMP_GET_NUM_THREADS()
        !$OMP END PARALLEL
        !$ print '(3X,"Number of workers ->  ",i2)',n_workers
    else
        write(*,*) "Program has been compiled without OpenMP; Time marching will run in serial"
    end if


    call cpu_time(t_cpu_0)
    call system_clock(count=t0, count_rate=ir)

    outname = "OUTPUT/snapshots.bin"
    filename          = "parameters.in"


    write(*,*) "##########################################"
    write(*,*) "######## Reading parameters file #########"
    write(*,*) "##########################################"

    filename="parameters.in"

    print*,"Is the parameters input file (parameters.in) [Yes/no]"
    read(*,*) filecheck

    if (filecheck=="Yes" .or. filecheck=="yes" .or. filecheck=="y" .or. &
            filecheck=="Y") then
        write(*,*) "Reading simulation parameters..."

    elseif  (filecheck=="No" .or. filecheck=="no" .or. filecheck=="n" .or. &
                filecheck=="N") then
        write(*,*) "Enter simulation parameters text file name with extension"
        write(*,*) "40 characters max"
        read(*,*) filename

    else
        write(*,*) "Only: Yes/yes/Y/y & No/no/N/n are handled"
        write(*,*) "The program have been terminated, please star over"
        stop
    end if

    open (2, file=filename, status = 'old')
    read(2,*) modnameprefix
    read(2,*) N
    read(2,*) ne
    read(2,*) h
    read(2,*) f0
    read(2,*) dt
    read(2,*) nt
    read(2,*) esrc
    read(2,*) gsrc
    read(2,*) isnap
    read(2,*) bc
    read(2,*) sbc
    read(2,*) gWidth
    read(2,*) attConst
    read(2,*) IC
    read(2,*) sd
    close(2)

    print*,"Polynomial order          -> ",N
    print*,"Number of elements        -> ",ne
    print*,"Element size              -> ",h
    print*,"Wavelet's peak frequency  -> ",f0
    print*,"Time step                 -> ",dt
    print*,"Number of time steps      -> ",nt
    print*,"Source location [nel/ngll]-> ",esrc, gsrc
    print*,"Snapshot interval         -> ",isnap

    ngll = N * ne + 1                  ! Total GLL points
    Jc = h / 2                         ! Jacobian for structured 1D mesh
    Jci = 1 / Jc                       ! Jacobian inverse

    allocate(Cij(N+1,ne))               ! Connectivity matrix
    allocate(xi(N+1))                   ! GLL points
    allocate(wi(N+1))                   ! GLL Quadrature weights
    allocate(v1D(ne))                   ! 1D velocity model in elements
    allocate(rho1D(ne))                 ! Density velocity model in elements
    allocate(rho1Dgll(ngll))            ! 1D density model mapped
    allocate(v1Dgll(ngll))              ! 1D velocity mapped
    allocate(M(ngll))                   ! Global mass matrix in vector form
    allocate(Minv(ngll,ngll))           ! Inverse of the mass matrix
    allocate(Me(N+1))                   ! Elemental mass matrix
    allocate(Kg(ngll,ngll))             ! Global stifness matrix
    allocate(Ke(N+1,N+1))               ! Elemental stifness matrix
    allocate(lprime(N+1,N+1))           ! Dervatives of Lagrange polynomials
    allocate(u(ngll))                   ! Displacement vector at time t
    allocate(unew(ngll),uddotnew(ngll), udotnew(ngll))                ! displacement vecotr at time t - dt
    allocate(uddot(ngll), udot(ngll), sigma(ngll))
    allocate(uold(ngll))                ! displacement vector at time t + dt
    allocate(src(nt))                   ! Source time function
    allocate(F(ngll))                   ! External force
    allocate(temp1(ngll),temp2(ngll),temp3(ngll))
    allocate(Uout(NINT(REAL(nt/isnap)),ngll))             ! Snapshots
    allocate(Udotout(NINT(REAL(nt/isnap)),ngll),sigmaout(NINT(REAL(nt/isnap)),ngll))
    allocate(mu1Dgll(ngll))             ! Shear modulus mapped
    allocate(xgll(ngll))                ! Array for global mapping
    allocate(g(ngll))
    allocate(as(NINT(nt/REAL(isnap)),ngll))
    allocate(tauL(ngll), recsigma(nt))

    Cij  = connectivity_matrix(N,ne)

    open(28,file="DBC.bin",access="direct",recl=nt*8)
    read(28,rec=1) recsigma
    close(28)

    call lagrangeprime(N,lprime)                         ! Lagrange polynomials derivatives
    call zwgljd(xi,wi,N+1,0.,0.)                         ! Getting GLL points and weights
    call readmodelfiles1D(v1D, rho1D, ne,modnameprefix)  ! Reading model files
    call shapefunc(N,h,ne,Cij,xgll)                      ! Global domain mapping
    call mapmodel(N,ne,rho1D,v1D,rho1Dgll,v1Dgll)        ! Mapping models
    call ricker(nt,f0,dt,src)                            ! Source time function
    

    write(*,*)"##########################################"
    write(*,*)"############### CFL Check ################"
    write(*,*)"##########################################"

    mindist = xgll(2) - xgll(1)
    CFL = (dt/mindist) * maxval(v1D(:))

    if (CFL > .8) then
        print"(a14,f6.3)","CFL value is ",CFL
        print*,"Decrease time step, the program has been terminated"
        stop

    else
        print"(a14,f6.3)","CFL value is ",CFL
        print*,"Simulation is stable"
    end if


    write(*,*)"##########################################"
    write(*,*)"########## Space Sampling check ##########"
    write(*,*)"##########################################"
    lambdamin = minval(v1D)/(f0*2.5)

    print"(a32,f3.1)", " Elements per minimum wavelength ->", lambdamin/h

    if ((lambdamin/h)<1) then
        print*,"Element size is too large"
        print*,"Numerical dispersion might be present"
        print*,"Do you wish to continue anyways? Yes/no"
        read*, filecheck

        if (filecheck=="Yes" .or. filecheck=="yes" .or. filecheck=="y" .or. &
                filecheck=="Y") then
            print*,"Proceeding..."

        elseif  (filecheck=="No" .or. filecheck=="no" .or. filecheck=="n" .or. &
                filecheck=="N") then
            write(*,*) "The program had been terminated"
            write(*,*) "Reduce element size or frequency as such"
            write(*,*) "you have at least 1 element per min wavelength"
            write(*,*) "This holds for N=5 to N=10 as per"
            write(*,*) "Komatitsch and Tromp 1999"
            stop
        else
            write(*,*) "Only: Yes/yes/Y/y & No/no/N/n are handled"
            write(*,*) "The program have been terminated, please star over"
            stop
        end if
    else
        print*, "Spatial sampling as OK!"
    end if


    !##########################################
    !####### Construct the mass matrix ########
    !##########################################

    do i=1,ne
        do j=1,N+1
            Me(j) = rho1Dgll(Cij(j,i)) * wi(j) * Jc   ! Elemental mass matrix construction
            M(Cij(j,i)) =  M(Cij(j,i)) + Me(j)        ! Filling the global mass matrix
        end do
    end do

    Minv(:,:) = 0
    ! Invert mass matrix
    do i=1,ngll
        Minv(i,i) = 1 / M(i)
    end do


    !###############################################
    !####### Construct the Stiffness matrix ########
    !###############################################
    mu1Dgll(:) =  rho1Dgll(:) * v1Dgll(:)**2.          ! Shear modulus
    Kg(:,:) = 0
    !$omp parallel do private(el,i,j,k,sum) shared(Kg,Ke,Cij,lprime,wi,Jc,Jci,N) schedule(static)
    do el=1,ne
        do i=1,N+1                                    ! Elemental stifness matrix
            do j=1,N+1
                sum = 0.
                do k=1,N+1
                    sum = sum + mu1Dgll(Cij(k,el)) * wi(k) * lprime(j,k) * lprime(i,k) * Jc * Jci**2
                end do
                Ke(i,j) = sum
            end do
        end do

        do i=1,N+1                                      ! Global assembly
            do j=1,N+1
                Kg(Cij(i,el),Cij(j,el)) = Kg(Cij(i,el),Cij(j,el)) + Ke(i,j)
            end do
        end do
    end do
    !$omp end parallel do


    g(:)  = 1

    do k=1,gWidth ! sponge layer
        taper       = exp(-(attConst*(gWidth-k)/gWidth)**2)
        g(k)        = taper
        if (sbc .ne. 1) then
            g(ngll-k+1) = taper 
            g(k)        = 1       
        end if
    end do

    write(*,*) "###########################################"
    write(*,*) "############ BEGIN TIME LOOP ##############"
    write(*,*) "###########################################"
    
    uold(:)     = 0.
    u(:)        = 0.
    unew(:)     = 0.
    uddot(:)    = 0.
    uddotnew(:) = 0.
    udot(:)     = 0.
    udotnew(:)  = 0.
    F(:)        = 0.
    tauL        = 0.


    if (IC==1) then
        u(:)    = exp(-1./sd**2*(xgll(:)-xgll(Cij(gsrc,esrc)))**2)
        uold(:) = u(:)
    end if


    k = 0;
    do t=1,nt

        tauL(Cij(N+1,200)) = recsigma(t)

        !if (IC .ne. 1) then
        !    F(Cij(gsrc,esrc)) =  src(t) * wi(gsrc) * Jc
        !end if

        unew(:)     = u(:) + dt * udot(:) + (dt**2)/2 * uddot(:)

        temp1(:) = 0
        !$OMP PARALLEL DO PRIVATE(k) SHARED(Kg,unew,temp1,ngll) SCHEDULE(static) 
        do l=1, ngll
            temp1(l) =  DOT_PRODUCT(Kg( l , : ), unew(:))
        end do
        !$OMP END PARALLEL DO

        !$OMP PARALLEL WORKSHARE
        temp2 = F(:) - temp1(:) + tauL(:)
        !$OMP END PARALLEL WORKSHARE

        !$OMP PARALLEL WORKSHARE
        uddotnew(:)  = (1/M(:)) * temp2(:)
        !$OMP END PARALLEL WORKSHARE

        !$OMP PARALLEL WORKSHARE
        udotnew(:)  = udot(:) + (dt/2) * (uddot(:) + uddotnew(:))
        !$OMP END PARALLEL WORKSHARE

        do el=1,ne
            do i=1,N+1
                tmp = 0.
                do j=1,N+1
                    tmp = tmp + mu1Dgll(Cij(i,el)) * u(Cij(j,el)) * lprime(j,i) * Jci    
                end do
                sigma(Cij(i,el)) = tmp
            end do 
        end do
        recsigma(t)  = sigma(Cij(N+1,200))
        !##########################################
        !##### Boundary Conditions            #####
        !##### 1: Free surface (implicit)     #####
        !##### 2: Rigid wall                  #####
        !##### 3: Periodic                    #####
        !##### 4: Sponge layer                #####
        !##########################################


        if (bc .eq. 2) then ! Rigid BC
            unew(1)    =   0

            if (sbc .ne. 1) then
                unew(ngll) = 0
            end if
        end if

        if (bc .eq. 3) then ! Periodic
                tmpBC      = unew(ngll)
                unew(ngll) = unew(1)
                unew(1)    = tmpBC
        end if

        if (bc .eq. 4) then ! Absorbing
            unew = unew * g
        end if


        u     = unew
        udot  = udotnew
        uddot = uddotnew

        if (mod(t,isnap) == 0) then
            k = k + 1
            Uout(k,:)    = u
            Udotout(k,:) = udot
            sigmaout(k,:)= sigma

            if (mod(t,NINT(nt/10.))==0) then
                print*,"At time sample ->",t, "/",nt

            end if
        end if
    end do


    !##########################################
    !##### Analytical Solution            #####
    !##########################################
    if (IC==1) then
        k = 1
        
        do t=1,nt,isnap
            !$OMP WORKSHARE
            as(k,:) = 1./2.*(exp(-1./sd**2 * (xgll(:)-xgll(Cij(gsrc,esrc)) + MAXVAL(v1D)*t*dt)**2)+ &
                             exp(-1./sd**2 * (xgll(:)-xgll(Cij(gsrc,esrc)) - MAXVAL(v1D)*t*dt)**2))
            !$OMP END WORKSHARE
            k = k + 1
        end do
    end if


    write(*,*) "##########################################"
    write(*,*) "######### Write solution binary ##########"
    write(*,*) "######### Solution in OUTPUT/   ##########"
    write(*,*) "##########################################"

    outname = "OUTPUT/SEM_snapshots_U.bin"

    inquire(iolength=reclsnaps) Uout
    open(15,file=outname,access="direct",recl=reclsnaps)
    write(15,rec=1) Uout
    close(15)


    outname = "OUTPUT/SEM_snapshots_V.bin"

    open(17,file=outname,access="direct",recl=reclsnaps)
    write(17,rec=1) Udotout
    close(17)


    outname = "OUTPUT/SEM_snapshots_Sigma.bin"

    open(18,file=outname,access="direct",recl=reclsnaps)
    write(18,rec=1) sigmaout
    close(18)


    !open(28,file="DBC.bin",access="direct",recl=nt*8)
    !write(28,rec=1) recsigma
    !close(28)


    if (IC==1) then
        reclsnaps = 0
        inquire(iolength=reclsnaps) as
    
        open(16,file="OUTPUT/SEM_snapshots_analytical_V.bin",access="direct",recl=reclsnaps)
        write(16,rec=1) as
        close(16)
    end if


   
    ! Temps elapsed final.
    call system_clock(count=t1, count_rate=ir)
    time = real(t1 - t0,kind=8) / real(ir,kind=8)

    call cpu_time(t_cpu_1)
    t_cpu = t_cpu_1 - t_cpu_0

    write(*,*) "##########################################"
    write(*,*) "######### TIME:                 ##########"
    print '(//3X,"Elapsed Time        : ",1PE10.3," [s]",/ &
            &,3X,"CPU Time            : ",1PE10.3," [s]",//)', &
            & time,t_cpu
    write(*,*) "##########################################"

    deallocate(xi,wi,v1D,v1Dgll,rho1D,rho1Dgll)
    deallocate(M,Me,Minv)
    deallocate(Cij)
    deallocate(mu1Dgll,xgll)
    deallocate(u,uold,unew,lprime,Kg,Ke)
    deallocate(Uout,Udotout,udot,udotnew,uddot,uddotnew)
    deallocate(sigma,sigmaout)
    deallocate(tauL, recsigma)

end program