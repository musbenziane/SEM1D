# SEM1D 
Basic 1D Spectral Elements Method Implementation Solver for the wave equation.
<br>Parameters input file should be as follow:<br>

Compile "create1Dmodel_files.f90" seperately to generate velocity and density 1D models (double precision floats)
<br> 
If OpenMP API is not installed, comment the OMP lines in the CMAKE file.
<br>
Check Cmake version before compiling 

<pre>
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
    
    </pre>


