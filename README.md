This repository contains the implementation of RAQPSO for optimization on Riemannian manifolds.
The algorithm has been encapsulated as a MATLAB function and can be directly applied to any test problem.

---
# Requirements

- MATLAB R2020 or later  
- Manopt toolbox  
  - Website: https://www.manopt.org  

Add Manopt to your MATLAB path:

```matlab
addpath(genpath('path/to/manopt'));
```

# Usage

- Define a Test Problem
  - **`problem.M`** — the manifold
  - **`problem.cost`** — the cost function
  - **`problem.egrad`** — (Optional) Euclidean gradient
Example: test_problem.m
```matlab
function problem = test_problem()
    % Stiefel manifold example
    n = 50; 
    k = 5;

    M = stiefelfactory(n, k);

    A = randn(n); 
    A = 0.5 * (A + A');  % Symmetric

    problem.M = M;
    problem.cost = @(X) -trace(X' * A * X);
    problem.egrad = @(X) -2 * A * X;
end
```

