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
___

![modeling_test_models](https://github.com/user-attachments/assets/28a8ba3e-b844-4617-b6a7-c0f7faa20a90)

![modeling_test_data](https://github.com/user-attachments/assets/1fbdbbe0-4a7a-48ee-9540-fced501ed6e2)

___
### Inversion test
___

![inversion_test_configuration](https://github.com/user-attachments/assets/e5f4e736-5d92-4010-8fa5-f007aabec6fa)

![inversion_test_results](https://github.com/user-attachments/assets/a420a288-e7dd-4b19-a291-6026423b7cf4)

![inversion_test_convergence](https://github.com/user-attachments/assets/c27ee1ca-ecb8-45a1-b3cc-68f3f1c40ba9)
___
### Migration test
___

![migration_test_results](https://github.com/user-attachments/assets/c62ac4e1-c4b8-423a-8872-8b46273b32ab)