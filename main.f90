! Created by mus on 28/07/2021.
! Spectral Elements 1D Solver with a regular mesh

program SEM1D
    !$ USE OMP_LIB

    implicit none

    interface connectivity_matrix
        function connectivity_matrix(N,ne)
            integer, intent(in)                                :: N, ne
            integer connectivity_matrix(N+1,ne)
        end function
    end interface


    real (kind=8)                                              :: Jc, Jci, h, f0, dt, sum, CFL, mindist, lambdamin
    real(kind=8)                                               :: time, t_cpu_0, t_cpu_1, t_cpu, tmp
    real (kind=8), dimension(:), allocatable                   :: xi, wi, v1D, rho1D, Me, M, xgll, rho1Dgll, v1Dgll
    real (kind=8), dimension(:), allocatable                   :: mu1Dgll, u, uold, unew, F, src, temp1, temp2, temp3
    real (kind=8), dimension(:,:), allocatable                 :: lprime, Minv, Kg, Ke,Uout
    integer, dimension(:,:), allocatable                       :: Cij
    integer                                                    :: N, ne, ngll, i, j, k, nt, isrc, t, el, isnap, reclsnaps
    integer                                                    :: ir, t0, t1
    character(len=40)                                          :: filename, filecheck, outname
    !$ integer                                                 :: n_workers


    !$OMP PARALLEL
    !$ n_workers = OMP_GET_NUM_THREADS()
    !$OMP END PARALLEL
    !$ print '(//,3X,"Number of workers ->  ",i2)',n_workers

    call cpu_time(t_cpu_0)
    call system_clock(count=t0, count_rate=ir)

    outname = "OUTPUT/snapshots.bin"

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
    read(2,*) N
    read(2,*) ne
    read(2,*) h
    read(2,*) f0
    read(2,*) dt
    read(2,*) nt
    read(2,*) isrc
    read(2,*) isnap
    close(2)

    print*,"Polynomial order         -> ",N
    print*,"Number of elements       -> ",ne
    print*,"Element size             -> ",h
    print*,"Wavelet's peak frequency -> ",f0
    print*,"Time step                -> ",dt
    print*,"Number of time steps     -> ",nt
    print*,"Source location          -> ",isrc
    print*,"Snapshot interval        -> ",isnap

    ngll = N * ne + 1                  ! Total GLL points
    Jc = h / 2                         ! Jacobian for structured 1D mesh
    Jci = 1 / Jc                       ! Jacobian inverse

    allocate(xi(N+1))                   ! GLL points
    allocate(wi(N+1))                   ! GLL Quadrature weights
    allocate(v1D(ne))                   ! 1D velocity model in elements
    allocate(rho1D(ne))                 ! Density velocity model in elements
    allocate(rho1Dgll(ngll))            ! 1D density model mapped
    allocate(v1Dgll(ngll))              ! 1D velocity mapped
    allocate(mu1Dgll(ngll))             ! Shear modulus mapped
    allocate(xgll(ngll))                ! Array for global mapping
    allocate(Cij(N+1,ne))               ! Connectivity matrix
    allocate(M(ngll))                   ! Global mass matrix in vector form
    allocate(Minv(ngll,ngll))           ! Inverse of the mass matrix
    allocate(Me(N+1))                   ! Elemental mass matrix
    allocate(Kg(ngll,ngll))             ! Global stifness matrix
    allocate(Ke(N+1,N+1))               ! Elemental stifness matrix
    allocate(lprime(N+1,N+1))           ! Dervatives of Lagrange polynomials
    allocate(u(ngll))                   ! Displacement vector at time t
    allocate(uold(ngll))                ! displacement vector at time t + dt
    allocate(unew(ngll))                ! displacement vecotr at time t - dt
    allocate(src(nt))                   ! Source time function
    allocate(F(ngll))                   ! External force
    allocate(temp1(ngll),temp2(ngll),temp3(ngll))
    allocate(Uout(NINT(REAL(nt/isnap)),ngll))             ! Snapshots


    call lagrangeprime(N,lprime)                   ! Lagrange polynomials derivatives
    call zwgljd(xi,wi,N+1,0.,0.)                   ! Getting GLL points and weights
    call readmodelfiles1D(v1D, rho1D, ne)          ! Reading model files
    call shapefunc(N,h,ne, xgll)                   ! Global domain mapping
    call mapmodel(N,ne,rho1D,v1D,rho1Dgll,v1Dgll)  ! Mapping models
    call ricker(nt,f0,dt,src)                      ! Source time function


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
    Cij  = connectivity_matrix(N,ne)

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
    mu1Dgll(:) = rho1Dgll(:) * v1Dgll(:)**2.          ! Shear modulus
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

    write(*,*) "!##########################################"
    write(*,*) "############ BEGIN TIME LOOP ##############"
    write(*,*) "###########################################"

    uold(:) = 0.
    u(:)    = 0.
    unew(:) = 0.
    F(:)    = 0.
    k = 0;
    do t=1,nt

        ! Injecting source
        F(isrc) = src(t)

        !$omp parallel do private(i,k,tmp) shared(Kg,u,temp1,ngll) schedule(static)
        do i=1,ngll
            tmp = 0.0
            do k=1,ngll
                tmp = tmp + kg(i,k) * u(k)
            enddo
            temp1(i) = tmp
        enddo
        !$omp end parallel do

        !$OMP PARALLEL WORKSHARE
        temp2 = F(:) - temp1(:)
        !$OMP END PARALLEL WORKSHARE

        !$omp parallel do private(i) shared(Minv,temp2,temp3,ngll) schedule(static)
        do i=1,ngll
            temp3(i) = Minv(i,i) * temp2(i)
        enddo
        !$omp end parallel do


        !$OMP PARALLEL WORKSHARE
        unew(:) = (dt**2.) * temp3(:)  + 2. * u(:) - uold(:)
        !$OMP END PARALLEL WORKSHARE
        uold = u
        u = unew

        if (mod(t,isnap) == 0) then
            k = k + 1
            Uout(k,:) = u
            if (mod(t,NINT(nt/100.))==0) then
                print*,"At time sample ->",t, "/",nt
            end if
        end if

    end do
    write(*,*) "##########################################"
    write(*,*) "######### Write solution binary ##########"
    write(*,*) "######### Solution in OUTPUT/   ##########"
    write(*,*) "##########################################"

    inquire(iolength=reclsnaps) Uout
    open(15,file=outname,access="direct",recl=reclsnaps)
    write(15,rec=1) Uout
    close(15)

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




    !deallocate(xi,wi,v1D,v1Dgll,rho1D,rho1Dgll,mu1Dgll,xgll,Cij,M,Me,Minv)
    deallocate(u,uold,unew,lprime,Kg,Ke)
    deallocate(Uout)
end program

