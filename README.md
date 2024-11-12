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

![modeling_test_models](https://github.com/user-attachments/assets/28a8ba3e-b844-4617-b6a7-c0f7faa20a90)

![modeling_test_data](https://github.com/user-attachments/assets/1fbdbbe0-4a7a-48ee-9540-fced501ed6e2)

### Inversion test

![inversion_test_configuration](https://github.com/user-attachments/assets/7ab8b1cc-2f42-4c1d-898a-e723e698abf2)

![inversion_test_results](https://github.com/user-attachments/assets/394695bf-fca6-474b-8168-93d49c0cfdb6)

![inversion_test_convergence](https://github.com/user-attachments/assets/b8da8946-fe18-4776-ae09-53052ccc7dcb)

### Migration test

![migration_test_results](https://github.com/user-attachments/assets/c62ac4e1-c4b8-423a-8872-8b46273b32ab)