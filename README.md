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

## Usage

### 1. Define a Test Problem

A valid test problem should include:

- **`problem.M`** — the manifold  
- **`problem.cost`** — the cost function  
- **(Optional) `problem.egrad`** — Euclidean gradient  
- **(Optional) initialization point**


