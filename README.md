This repository contains the implementation of RAQPSO for optimization on Riemannian manifolds.
The algorithm has been encapsulated as a MATLAB function and can be directly applied to any test problem.

---

## Requirements

- MATLAB R2020 or later  
- Manopt toolbox  
  - Website: https://www.manopt.org  

Add Manopt to your MATLAB path:

```matlab
addpath(genpath('path/to/manopt'));

## Usage

- Define a Test Problem
  - **`problem.M`** — the manifold
  - **`problem.cost`** — the cost function
  - **`problem.egrad`** — (Optional) Euclidean gradient
  - **Initialization point** — (Optional) starting point



