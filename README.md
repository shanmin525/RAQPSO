This repository contains the implementation of RAQPSOfor optimization on Riemannian manifolds.
The algorithm has been encapsulated as a MATLAB function and can be directly applied to any test problem.


---

## Requirements

- MATLAB R2020 or later  
- Manopt toolbox  
  - Website: https://www.manopt.org  

Add Manopt to your MATLAB path:

```matlab
addpath(genpath('path/to/manopt'));


Install Manopt
1.Either git-clone or unzip the whole manopt directory you just downloaded in a location of your choice, say, in /my/directory/.
2.Go to /my/directory/manopt/ at the Matlab prompt and execute importmanopt.
3.You will be asked whether you want to save this path for your next Matlab sessions. If you reply Y (yes) and you have the rights to run savepath, then you won’t need to go through these steps next time you open Matlab.
