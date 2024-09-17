#!/bin/bash

# Input Output scripts --------------------------------------------------------

ioFunctions="../src/ioFunctions/ioFunctions.cpp"

# Acquisition geometry scripts ------------------------------------------------

geometry="../src/geometry/geometry.cpp"

# Seismic modeling scripts ----------------------------------------------------

eikonal="../src/modeling/hfreq/eikonal.cpp"
serial_aFSM="../src/modeling/hfreq/isotropic/serial_aFSM.cpp"
parallel_aFSM="../src/modeling/hfreq/isotropic/parallel_aFSM.cu"
podvin_lecomte="../src/modeling/hfreq/isotropic/podvin_lecomte.cu"

scalar="../src/modeling/lfreq/isotropic/scalar.cu"
elastic="../src/modeling/lfreq/isotropic/elastic.cu"
acoustic="../src/modeling/lfreq/isotropic/acoustic.cu"
wavefield="../src/modeling/lfreq/wavefield.cu"

modeling="../src/modeling/modeling.cpp"

modeling_all="$modeling $serial_aFSM $parallel_aFSM $podvin_lecomte 
              $eikonal $scalar $elastic $acoustic $wavefield"

modeling_main="../src/modeling_main.cpp"

# Seismic inversion scripts ---------------------------------------------------

# least_squares="../src/inversion/least_squares.cpp"
# adjoint_state="../src/inversion/adjoint_state.cu"

# tomography="../src/inversion/tomography.cpp"

# inversion_all="$tomography $least_squares $adjoint_state"

# inversion_main="../src/inversion_main.cpp"

# Seismic migration scripts ---------------------------------------------------

# kirchhoff="../src/migration/kirchhoff.cpp"

# migration_main="../src/main/migration_main.cpp"

# migration_all="$kirchhoff"

# Compiler flags --------------------------------------------------------------

flags="--std=c++11 -lm -w -g -O3"

# Main dialogue ---------------------------------------------------------------

USER_MESSAGE="
-------------------------------------------------------------------------------
                                 \033[34mSeisFAT2D\033[0;0m
-------------------------------------------------------------------------------
\nUsage:\n
    $ $0 -compile             # Generate executables 
    $ $0 -modeling            # Run eikonal equation solver          
    $ $0 -inversion           # Run first arrival tomography
    $ $0 -migration           # Run kirchhoff depth migration   
-------------------------------------------------------------------------------
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

    echo -e "../bin/\033[34mmodeling.exe\033[m" 
    nvcc $ioFunctions $geometry $modeling_all $modeling_main $flags -o ../bin/modeling.exe

    # echo -e "../bin/\033[34minversion.exe\033[m" 
    # nvcc $ioFunctions $geometry $modeling_all $inversion_all $inversion_main $flags -o ../bin/inversion.exe

    # echo -e "../bin/\033[34mmigration.exe\033[m"
    # nvcc $ioFunctions $geometry $modeling_all $migration_all $migration_main $flags -o ../bin/migration.exe

	exit 0
;;

-modeling) 

    ./../bin/modeling.exe parameters.txt
	
    exit 0
;;

-inversion) 
    
    # ./../bin/inversion.exe parameters.txt
	
    exit 0
;;

-migration) 
    
    # ./../bin/migration.exe parameters.txt
	
    exit 0
;;

* ) 

	echo -e "\033[31mERRO: Option $1 unknown!\033[m"
	echo -e "\033[31mType $0 -h for help \033[m"
	
    exit 3
;;

esac