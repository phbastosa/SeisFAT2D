___

## Seismic First-Arrival Toolkit : SeisFAT2D

Modeling, inversion and migration using massive computational parallelism in object-oriented programming.

### Requirements:

- A linux distribution or the Windows Subsystem for linux
- A Nvidia Graphic Processing Unit (GPU)
- Cuda 12 or higher with compatible drivers
- C++ and Python 3 programming languages
- FFTW3 libary ---> https://www.fftw.org/
____

### Usage:

```console
SeisFAT2D$ cd run/
SeisFAT2D/run/$ ./program.sh -h

Usage:

    $ ./program.sh -compile              
    $ ./program.sh -modeling                      
    $ ./program.sh -inversion           
    $ ./program.sh -migration           

Tests:

    $ ./program.sh -test_modeling                 
    $ ./program.sh -test_inversion      
    $ ./program.sh -test_migration      
```

# Test results

### Modeling test

![modeling_test_models](https://github.com/user-attachments/assets/86667335-d7f0-4b29-84ad-3118901ce334)

![modeling_test_data](https://github.com/user-attachments/assets/1fbdbbe0-4a7a-48ee-9540-fced501ed6e2)

### Inversion test

![inversion_test_configuration](https://github.com/user-attachments/assets/f62e36cd-96c8-41ea-941b-e0a628e69d25)

![inversion_test_results](https://github.com/user-attachments/assets/11be0be4-e8a8-4564-99fe-59f0716103fb)

![inversion_test_convergence](https://github.com/user-attachments/assets/5782be67-6289-4dd7-8730-ac400b0f1f33)

### Migration test

![migration_test_results](https://github.com/user-attachments/assets/804a8791-a949-4117-8442-2fd966a37d61)