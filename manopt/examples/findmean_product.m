function mbest = findmean_product(problem, id, pbest_pop, pbest, popsize, rtt)
% FINDMEAN_PRODUCT: weighted mean computation on product/single manifold
% relative to particle pbest_pop{id}.
%
% Supports both:
%   - Product manifolds (cell structure)
%   - Submanifolds whose tangent vectors are structs (e.g., Stiefel)
%
% Inputs:
%   problem.M : manifold object supporting proj() and retr()
%   id        : current particle index
%   pbest_pop : cell array of all particles {1:popsize}
%   pbest     : fitness array of all particles
%   popsize   : number of particles
%   rtt       : retraction type or step parameter
%
% Output:
%   mbest     : weighted mean point on manifold

k = popsize;
index = randperm(popsize);

% --- Compute normalized inverse-fitness weights ---
fit = pbest(index);
ww = max(fit) - fit + 1e-4;
ww = ww / sum(ww);

x0 = pbest_pop{id};

% --- Determine tangent structure type ---
proj_sample = problem.M.proj(x0, pbest_pop{index(1)});

% --- Initialize weighted direction container ---
if iscell(proj_sample)
    numBlocks = numel(proj_sample);
    meanw = cell(1, numBlocks);
    for j = 1:numBlocks
        meanw{j} = initialize_like(proj_sample{j});
    end
else
    meanw = initialize_like(proj_sample);
end

% --- Weighted tangent accumulation ---
for i = 1:k
    proj_i = problem.M.proj(x0, pbest_pop{index(i)});
    if iscell(proj_i)
        for j = 1:numel(proj_i)
            meanw{j} = add_weighted(meanw{j}, proj_i{j}, ww(i));
        end
    else
        meanw = add_weighted(meanw, proj_i, ww(i));
    end
end

% --- Retraction back to manifold ---
mbest = problem.M.retr(x0, meanw, rtt);
end

%% --- helper: initialize a zero-like object (matrix, struct, or cell)
function out = initialize_like(template)
    if isstruct(template)
        fields = fieldnames(template);
        for f = 1:numel(fields)
            val = template.(fields{f});
            if isnumeric(val)
                out.(fields{f}) = zeros(size(val));
            else
                out.(fields{f}) = val; % non-numeric fields kept unchanged
            end
        end
    elseif iscell(template)
        out = cell(size(template));
        for k = 1:numel(template)
            out{k} = initialize_like(template{k});
        end
    elseif isnumeric(template)
        out = zeros(size(template));
    else
        error('Unsupported tangent vector type.');
    end
end

%% --- helper: weighted sum of tangent vectors
function out = add_weighted(accum, term, weight)
    if isstruct(term)
        out = accum;
        fields = fieldnames(term);
        for f = 1:numel(fields)
            if isnumeric(term.(fields{f}))
                out.(fields{f}) = accum.(fields{f}) + weight * term.(fields{f});
            end
        end
    elseif iscell(term)
        out = accum;
        for k = 1:numel(term)
            out{k} = add_weighted(accum{k}, term{k}, weight);
        end
    elseif isnumeric(term)
        out = accum + weight * term;
    else
        error('Unsupported data type in weighted addition.');
    end
end
