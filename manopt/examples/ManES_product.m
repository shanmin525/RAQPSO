function [Ybest, fbest, info, S] = ManES_product(problem, Y, options)
% ManES_product - Manifold Evolution Strategy adapted to product/struct-type points.
% Usage:
%   [Ybest, fbest, info, S] = ManES_product(problem, Y, options)
%
% This version supports problem.M being a product manifold where points Y are
% structures with multiple fields (e.g., Y.A, Y.R). The code flattens such
% structures to vectors for sampling and recombination, and converts back when
% calling manifold operations (proj, retr, transp, cost, etc.).
%
% Author: adapted from original ManES by Xiaoyu He / 74006
% Date: adapted for product manifold.

timetic = tic;

if ~exist('Y','var') || isempty(Y)
    Y = problem.M.rand();
end

% Field info and flatten/unflatten helpers
[fieldNames, fieldSizes, fieldLengths, D] = extract_fieldinfo(Y);

% Evaluate dimension and shapes
dim = problem.M.dim(); % manifold dimension
% D is total ambient elements count (sum of prod of field sizes)
% For compatibility with original naming, treat D analogous to n*p
fbest = problem.cost(Y);
Ybest = Y;

% default options (same as original)
localdefaults = struct('minsigma',1e-6,'maxiter',inf,...
    'maxcostevals',D*10000,...
    'sigma0',1,...
    'lambda',4+floor(3*log(dim)),...
    'm',['4+floor(3*log(dim))'],...
    'damping',1,...
    'verbosity',0);
options = mergeOptions(mergeOptions(getGlobalDefaults(), localdefaults), options);

lambda = options.lambda;
mu = ceil(lambda/2);
weights = log(mu+1/2)-log(1:mu)'; weights = weights/sum(weights);
mueff = 1/sum(weights.^2);

len = 1:2*lambda; scum = 0;
cs = 0.3;
damping = options.damping;
qtarget = 0.05;
sigma_ = options.sigma0;

prevfits = fbest*ones(1, lambda);

ccov = 1/(3*sqrt(dim)+5);
m = eval(options.m);
PCs = zeros(D, m); % principal components stored in flattened form
cc = lambda/dim./4.^(0:m-1);

% Preallocate containers
Z = zeros(D, 2*lambda); % will fill columns
Yi = cell(1, lambda);
fit_ = zeros(1, lambda);

FEs = 1; iter = 0;
recordGap = 1; recordIter = 1;
recordLength = ceil(options.maxcostevals/lambda/recordGap)+10;
record = zeros(recordLength,5);
record(recordIter,:) = [iter, FEs, fbest, sigma_, toc(timetic)];
storedb = StoreDB(); % keep minimal

% main loop
while iter < options.maxiter && FEs < options.maxcostevals && sigma_ > options.minsigma
    iter = iter + 1;
    %% sample
    % rank-1 line-Gaussian perturbations in D-dimensional tangent space
    % sample half set
    half = ceil(lambda/2);
    Z_half = PCs * ( sqrt(ccov*(1-ccov).^(0:m-1))' .* randn(m, half) );
    for i = 1:half
        % isotropic ambient Gaussian of same structure, then project
        ambient_rand_vec = randn(D,1);
        ambient_rand_struct = vec2struct(ambient_rand_vec, fieldNames, fieldSizes, fieldLengths);
        % project to tangent at Y
        t_struct = problem.M.proj(Y, ambient_rand_struct);
        t_vec = struct2vec(t_struct, fieldNames, fieldSizes, fieldLengths);
        Z_half(:,i) = Z_half(:,i) + sqrt((1-ccov)^m) * t_vec;
    end
    % mirror to get symmetric set as in original
    Z(:,1:half) = Z_half;
    Z(:,half+1:lambda) = -Z_half(:,1:floor(lambda/2));
    % If lambda is odd, we may have one extra column; but original used symmetric
    
    %% evaluate and sort
    for i = 1:lambda
        % reshape to struct tangent
        tangent_struct = vec2struct(Z(:,i), fieldNames, fieldSizes, fieldLengths);
        % retract (with step sigma_)
        % Some M.retr may not accept sigma; handle both cases
        try
            Yi{i} = problem.M.retr(Y, tangent_struct, sigma_);
        catch
            Yi{i} = problem.M.retr(Y, tangent_struct); % fallback
        end
        fit_(i) = problem.cost(Yi{i});
    end
    FEs = FEs + lambda;
    [~, sortedIdx] = sort(fit_);
    
    %% recombine and move
        %% recombine and move
    MeanZ = Z(:, sortedIdx(1:mu)) * weights;
    MeanZ_struct = vec2struct(MeanZ, fieldNames, fieldSizes, fieldLengths);
    % --- store previous point before retraction ---
    Y_prev = Y;
    try
        Y = problem.M.retr(Y_prev, MeanZ_struct, sigma_);
    catch
        Y = problem.M.retr(Y_prev, MeanZ_struct);
    end

    %% adapt distribution (update PCs)
    PCs = (1-cc) .* PCs + real(sqrt(cc.*(2-cc)*mueff)) .* MeanZ;
    % move PCs to new tangent space by transporting each PC column
    for i = 1:m
        pc_vec = PCs(:,i);
        pc_struct = vec2struct(pc_vec, fieldNames, fieldSizes, fieldLengths);
        % --- use Y_prev and Y as source and target in transp ---
        pc_transp = problem.M.transp(Y_prev, Y, pc_struct); 
        pc_vec2 = struct2vec(pc_transp, fieldNames, fieldSizes, fieldLengths);
        PCs(:,i) = pc_vec2;
    end

    
    % adapt step size (same as original)
    [~, imax] = sort([prevfits, fit_]);
    R1 = sum(len(imax<=lambda));
    U1 = R1 - lambda*(lambda+1)/2;
    W = U1 - lambda^2/2;
    W = W / sqrt(lambda^2*(2*lambda+1)/12);
    scum = (1-cs)*scum + sqrt(cs*(2-cs))*W;
    sigma_ = sigma_ * exp((normcdf(scum)/(1-qtarget)-1)/damping);
    prevfits = fit_;
    
    %% trace elite
    if fit_(sortedIdx(1)) < fbest
        fbest = fit_(sortedIdx(1));
        Ybest = Yi{sortedIdx(1)};
    end
    
    %% save records
    record(recordIter,:) = [iter, FEs, fbest, sigma_, toc(timetic)];
    recordIter = recordIter + 1;
end

% finalize record
record(recordIter:end,:) = [];
info = array2table(record, 'VariableNames', {'iter','costevals','cost','sigma','time'});
info = table2struct(info);
S = []; % placeholder for backward compatibility

% compute grad/hess at Ybest (optional)
storedb = StoreDB();
key = storedb.getNewKey();
% getGradient/getHessian accept structured Ybest
try
    grad = getGradient(problem, Ybest, storedb, key);
    hess = getHessian(problem, Ybest, grad, storedb, key);
catch
    grad = []; hess = [];
end

end

%%%%%% Helper functions %%%%%%

function [fieldNames, fieldSizes, fieldLengths, total] = extract_fieldinfo(S)
% returns ordered list of field names and sizes, lengths and total length
fieldNames = fieldnames(S);
numf = numel(fieldNames);
fieldSizes = cell(numf,1);
fieldLengths = zeros(numf,1);
total = 0;
for k = 1:numf
    val = S.(fieldNames{k});
    sz = size(val);
    fieldSizes{k} = sz;
    fieldLengths(k) = prod(sz);
    total = total + fieldLengths(k);
end
end

function v = struct2vec(S, fieldNames, fieldSizes, fieldLengths)
% flatten struct fields into column vector in order of fieldNames
numf = numel(fieldNames);
v = zeros(sum(fieldLengths),1);
idx = 1;
for k = 1:numf
    f = fieldNames{k};
    fld = S.(f);
    len = fieldLengths(k);
    if len == 0
        continue;
    end
    v(idx:idx+len-1) = fld(:);
    idx = idx + len;
end
end

function S = vec2struct(v, fieldNames, fieldSizes, fieldLengths)
% convert vector back to struct with same field shapes
numf = numel(fieldNames);
S = struct();
idx = 1;
for k = 1:numf
    len = fieldLengths(k);
    sz = fieldSizes{k};
    if len == 0
        S.(fieldNames{k}) = zeros(sz);
    else
        block = v(idx:idx+len-1);
        S.(fieldNames{k}) = reshape(block, sz);
        idx = idx + len;
    end
end
end
