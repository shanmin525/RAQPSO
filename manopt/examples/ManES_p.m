function [Ybest, fbest, info] = ManES_p(problem, Y, options)

if ~exist('Y','var') || isempty(Y)
    Y = problem.M.rand();   % Y is struct: Y.A, Y.R
end

fbest = problem.cost(Y);
Ybest = Y;

lambda = options.lambda;
mu = ceil(lambda/2);
weights = log(mu+1/2)-log(1:mu)';
weights = weights/sum(weights);

for iter = 1:options.maxiter
    
    %%% === 1. Sample perturbations ===
    for i = 1:lambda
        U = struct();
        U.A = randn(size(Y.A));
        U.R = randn(size(Y.R));
        U = problem.M.proj(Y, U);
        Z{i} = U;
    end
    
    %%% === 2. Generate candidates ===
    for i = 1:lambda
        Yi{i} = problem.M.retr(Y, Z{i}, options.sigma0);
        fit(i) = problem.cost(Yi{i});
    end
    
    %%% === 3. Selection & recombination ===
    [~,idx] = sort(fit);
    MeanZ = weighted_sum(Z(idx(1:mu)), weights, problem.M, Y);
    Y = problem.M.retr(Y, MeanZ, options.sigma0);
    
    %%% === 4. Elitism ===
    if fit(idx(1)) < fbest
        fbest = fit(idx(1));
        Ybest = Yi{idx(1)};
    end
end

info = [];
end
