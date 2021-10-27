# SEM1D 
Basic 1D Spectral Elements Method Implementation Solver for the wave equation.
<br>Parameters input file should be as follow:<br>

4                    ! Polynomial order <br>
1000                 ! Number of elements<br>
3                    ! Element size<br>
100.                 ! Wavelet's peak frequency<br>
0.00007215           ! Time step<br>
1000                 ! Number of time steps<br>
1500                 ! Source location<br>
50                   ! Snapshot interval<br>

<h3> Build with gfortran </h3>
git clone or just download it<br> 
mkdir src <br>
unzip inside src/<br>
mkdir build && cd build/<br>
cmake ../src<br>
Make <br>

Compile "create1Dmodel_files.f90" seperately to generate velocity and density 1D models (double precision floats)
<br> 
OpenMP lines in the CMAKE file are for learning OMP purposes, omitting my prodcuce errors.
<br>
Check Cmake version before compiling 
