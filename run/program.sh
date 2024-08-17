#!/bin/bash

# Input Output scripts --------------------------------------------------------------------------------

ioFunctions="../src/ioFunctions/ioFunctions.cpp"

# Acquisition geometry scripts ------------------------------------------------------------------------

geometry="../src/geometry/geometry.cpp"

# Seismic modeling scripts ----------------------------------------------------------------------------

eikonal="../src/modeling/eikonal.cu"

modeling="../src/modeling/modeling.cpp"

modeling_main="../src/modeling_main.cpp"

modeling_all="$modeling $eikonal"

# Seismic inversion scripts ---------------------------------------------------------------------------

inversion="../src/inversion/inversion.cpp"

least_squares="../src/inversion/least_squares_tomography.cu"
adjoint_state="../src/inversion/adjoint_state_tomography.cu"

inversion_main="../src/main/inversion_main.cpp"

inversion_all="$inversion $least_squares $adjoint_state"

# Seismic migration scripts ---------------------------------------------------------------------------

migration="../src/migration/migration.cpp"

kirchhoff="../src/migration/kirchhoff/kirchhoff.cu"

migration_main="../src/main/migration_main.cpp"

migration_all="$migration $kirchhoff"

# Compiler flags --------------------------------------------------------------------------------------

flags="--std=c++11 -lm -w -g -O3"

# Main dialogue ---------------------------------------------------------------------------------------

USER_MESSAGE="-------------------------------------------------------------------------------
                                 \033[34mSeisFAT2D\033[0;0m
-------------------------------------------------------------------------------
\nUsage:\n
    $ $0 -compile             # Generate executables 
    $ $0 -modeling            # Run eikonal solver          
    $ $0 -inversion           # Run first arrival tomography
    $ $0 -migration           # Run kirchhoff depth migration   

Tests:\n
    $ $0 -test_modeling       # Perform a small modeling experiment          
    $ $0 -test_inversion      # Perform a small inversion experiment
    $ $0 -test_migration      # Perform a small migration experiment          
"

[ -z "$1" ] && 
{
	echo -e "\nYou didn't provide any parameter!" 
	echo -e "Type $0 -help for more info\n"
    exit 1 
}

case "$1" in

-h) 

	echo -e "$USER_MESSAGE"
	exit 0
;;

-compile) 

    echo -e "Compiling stand-alone executables!\n"

    echo -e "../bin/\033[31mmodeling.exe\033[m" 
    nvcc $ioFunctions $geometry $modeling_all $modeling_main $flags -o ../bin/modeling.exe

    # echo -e "../bin/\033[31minversion.exe\033[m" 
    # nvcc $ioFunctions $geometry $modeling_all $inversion_all $inversion_main $flags -o ../bin/inversion.exe

    # echo -e "../bin/\033[31mmigration.exe\033[m"
    # nvcc $ioFunctions $geometry $modeling_all $migration_all $migration_main $flags -o ../bin/migration.exe

	exit 0
;;

-modeling) 

    ./../bin/modeling.exe parameters.txt
	
    exit 0
;;

-inversion) 
    
    ./../bin/inversion.exe parameters.txt
	
    exit 0
;;

-migration) 
    
    ./../bin/migration.exe parameters.txt
	
    exit 0
;;

-test_modeling)

    python3 ../tests/modeling/eikonal_equation/generate_models.py

    spacings=(100 50 25)
    methods=("pod" "fim" "fsm")

    for method in ${methods[@]}; do 
        for spacing in ${spacings[@]}; do 
            ./../bin/modeling.exe ../tests/modeling/eikonal_equation/parFiles/parameters_"$method"_"$spacing"m.txt; 
        done    
    done 

    python3 ../tests/modeling/eikonal_equation/generate_figures.py

    python3 ../tests/modeling/wave_equation/generate_models.py 

    ./../bin/modeling.exe ../tests/modeling/wave_equation/parFiles/parameters_scalar_homogeneous.txt
    ./../bin/modeling.exe ../tests/modeling/wave_equation/parFiles/parameters_acoustic_homogeneous.txt
    ./../bin/modeling.exe ../tests/modeling/wave_equation/parFiles/parameters_elastic_homogeneous.txt

    python3 ../tests/modeling/wave_equation/generate_figures.py 

	exit 0
;;

-test_inversion) 

    python3 ../tests/inversion/generate_models.py

    ./../bin/modeling.exe ../tests/inversion/parFiles/parameters_obsData.txt

    ./../bin/inversion.exe ../tests/inversion/parFiles/parameters_leastSquares.txt
    ./../bin/inversion.exe ../tests/inversion/parFiles/parameters_adjointState.txt

    ./../bin/modeling.exe ../tests/inversion/parFiles/parameters_lsFinalModeling.txt
    ./../bin/modeling.exe ../tests/inversion/parFiles/parameters_adjFinalModeling.txt

    python3 ../tests/inversion/generate_figures.py
	
    exit 0
;;

-test_migration)

    python3 ../tests/migration/generate_models.py

    ./../bin/modeling.exe ../tests/migration/parFiles/parameters_mod_diffraction.txt
    ./../bin/modeling.exe ../tests/migration/parFiles/parameters_mod_homogeneous.txt

    python3 ../tests/migration/prepare_data.py

    ./../bin/migration.exe ../tests/migration/parFiles/parameters_migration.txt

    python3 ../tests/migration/generate_figures.py

	exit 0
;;

-check_geometry)

    ./../bin/geometry.exe parameters.txt

    python3 ../src/visualization/check_geometry.py parameters.txt

	exit 0
;;

* ) 

	echo -e "\033[31mERRO: Option $1 unknown!\033[m"
	echo -e "\033[31mType $0 -h for help \033[m"
	
    exit 3
;;

esac