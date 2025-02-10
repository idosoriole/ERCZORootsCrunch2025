## Supporting resources for 'Deep roots supply reactivity and enhance silicate weathering in the bedrock vadose zone'

Osorio-Leon, I.D., Rempe, D.M., Golla, J.K., Bouchez, and Druhan, J.L. (2025) 'Deep roots supply reactivity and enhance silicate weathering in the bedrock vadose zone', AGU Advances.

### Batch_Experiments_Model
Directory that contains input and selected output files for the batch simulations of the batch experiments. The models were run with a version of CrunchTope used in Golla et al. (2021; [https://doi.org/10.1016/j.epsl.2021.116988](https://doi.org/10.1016/j.epsl.2021.116988 "https://doi.org/10.1016/j.epsl.2021.116988")) and with PETSc 3.7.7 ([https://petsc.org/release/](https://petsc.org/release/ "https://petsc.org/release/")).

### 1D_model
Directory that contains the input and selected output files for the 1D RTM simulations used to build figures 3 and 4. The models were run with a version of CrunchTope used in Golla et al. (2021; [https://doi.org/10.1016/j.epsl.2021.116988](https://doi.org/10.1016/j.epsl.2021.116988 "https://doi.org/10.1016/j.epsl.2021.116988")) and with PETSc 3.7.7 ([https://petsc.org/release/](https://petsc.org/release/ "https://petsc.org/release/")).

### ReadRTM.ipynb
Jupyter notebook to treat the output files in the Batch_Experiments_Model and 1D_model directories and reproduce main computations reported in the manuscript. 

### Krunchmanday.py
Set of functions to process Crunchflow output files and manipulate data. Imported by ReadRTM.ipynb
