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
### Citation:

```console
@software{paulo_h_b_alves_2024_14170687,
  author       = {Paulo H. B. Alves},
  title        = {phbastosa/SeisFAT2D: Initial},
  month        = November,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {1.0},
  doi          = {10.5281/zenodo.14170687},
  url          = {https://doi.org/10.5281/zenodo.14170687}
}
```
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

![modeling_test_accuracy](https://github.com/user-attachments/assets/5c02d496-b7e0-4e07-9ab9-7710a89497f5)

___
### Inversion test
___

![inversion_test_configuration](https://github.com/user-attachments/assets/e5f4e736-5d92-4010-8fa5-f007aabec6fa)

![inversion_test_results](https://github.com/user-attachments/assets/46bb6452-5d6e-4c53-838c-8d8b3bf2d156)

![inversion_test_convergence](https://github.com/user-attachments/assets/f1a3a417-3e9d-46cd-be25-a87462f8c9a4)

___
### Migration test
___

![migration_test_results](https://github.com/user-attachments/assets/c62ac4e1-c4b8-423a-8872-8b46273b32ab)
