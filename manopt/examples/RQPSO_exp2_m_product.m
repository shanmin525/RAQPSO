function [Ybest, fbest, info, gbestlog] = RQPSO_exp2_m_product(problem, options, beta, e, ps)
% RQPSO_exp2_m_product_fast
% 高效版 Riemannian QPSO 算法 (Product Manifold)
% 减少 struct2vec / vec2struct 调用，显著提升速度
%
% 主要优化：
%   - 避免在每次更新中重复 struct2vec/vec2struct；
%   - 保持 retraction / projection / transport 完整调用；
%   - 向量化噪声和缩放操作。

% -------------------- 初始化 --------------------
dim = problem.M.dim();
Y0 = problem.M.rand();
[fieldNames, fieldSizes, fieldLengths, D] = extract_fieldinfo(Y0);

localdefaults = struct('minsigma',1e-6, 'maxiter', D*200, 'damping',1);
options = mergeOptions(mergeOptions(getGlobalDefaults(), localdefaults), options);

switch ps
    case 1, popsize = 4*floor(log(dim));
    case 2, popsize = 10;
    case 3, popsize = 30;
    case 4, popsize = 50;
    otherwise, popsize = max(10, round(4*log(dim)));
end

maxcostevals = D * 10000;
iter = 1;

% 种群初始化
pop = cell(popsize,1);
for i = 1:popsize
    pop{i} = problem.M.rand();
end
pbest_pop = pop;
pbestvalue = inf(1,popsize);
fit = zeros(1,popsize);

% 初始评估
for i = 1:popsize
    fit(i) = problem.cost(pop{i});
    pbestvalue(i) = fit(i);
end

[~, index] = sort(pbestvalue);
Ybest = pbest_pop{index(1)};
fbest = pbestvalue(index(1));
FES(1) = popsize;
gbestlog = fbest;
maxiter = ceil(D * 10000 / popsize);
rtt = 1 - (1:maxiter)/maxiter;

% -------------------- 主循环 --------------------
while FES(iter) <= maxcostevals
    pop_pre = pop;
    prevfits = fit;

    % ----------- 计算群体平均中心 mbest -----------
    if FES(iter)<e*maxcostevals
        tt = randi(popsize);
        mbest = findmean_product(problem, tt, pbest_pop, pbestvalue, popsize, rtt(iter));
    else
        mbest = findmean1_product(problem, Ybest, pbest_pop, pbestvalue, popsize, rtt(iter));
    end

    % ----------- 预提取 Ybest 的向量形式（仅一次）-----------
    y_vec = struct2vec(Ybest, fieldNames, fieldLengths);

    % ----------- 更新每个粒子 -----------
    for i = 1:popsize
        % Step 1: pbest → Ybest 的方向
%         diff_struct = problem.M.log(pbest_pop{i}, Ybest);
        pv_struct = problem.M.proj(pbest_pop{i}, Ybest);

        % Step 2: pbest 上回撤小步
        pop_i_base = problem.M.retr(pbest_pop{i}, pv_struct, 0.1*rand());

        % Step 3: 在欧氏空间中构造扰动（只转一次向量）
        pv_vec = struct2vec(pv_struct, fieldNames, fieldLengths);
        scale = e .* abs(pv_vec) + eps;
        tmpvec_vec = y_vec + scale .* randn(D,1);
        tmpvec_struct = vec2struct(tmpvec_vec, fieldNames, fieldSizes, fieldLengths);

        % Step 4: 投影 & 回撤
        tmp_proj = problem.M.proj(pop_i_base, tmpvec_struct);
        pop_temp = problem.M.retr(pop_i_base, tmp_proj, 0.1*rand());

        % Step 5: pop{i} → mbest 向量
%         diff2_struct = problem.M.log(pop{i}, mbest);
        proj_popi_struct = problem.M.proj(pop{i}, mbest);

        % Step 6: 向量传输到 pop_temp
        p_temp2_struct = problem.M.transp(pop{i}, pop_temp, proj_popi_struct);

        % Step 7: 随机缩放与符号
        step_vec = beta .* log(1./rand(D,1)) .* struct2vec(p_temp2_struct, fieldNames, fieldLengths);
        if rand()<0.5
            step_vec = -step_vec;
        end
        step_struct = vec2struct(step_vec, fieldNames, fieldSizes, fieldLengths);

        % Step 8: 最终回撤更新
        pop{i} = problem.M.retr(pop_temp, step_struct, rtt(min(iter,numel(rtt))));
    end

    % ----------- 更新适应度与最优解 -----------
    for i = 1:popsize
        fit(i) = problem.cost(pop{i});
        if fit(i) < pbestvalue(i)
            pbestvalue(i) = fit(i);
            pbest_pop{i} = pop{i};
        end
    end
    [minp, idxmin] = min(pbestvalue);
    if minp < fbest
        fbest = minp;
        Ybest = pbest_pop{idxmin};
    end
    gbestlog(iter) = fbest;

    % ----------- 精英保留合并 -----------
    combined_fits = [prevfits, fit];
    pop_sum = [pop_pre, pop];
    [~, imax] = sort(combined_fits);
    for k = 1:popsize
        pop{k} = pop_sum{imax(k)};
    end

    % ----------- 计数与终止条件 -----------
    FES(iter+1) = FES(iter) + popsize;
    iter = iter + 1;
    if iter > maxiter + 5
        break;
    end
end

% 输出信息
info = struct();
info.FES = FES;
info.best = gbestlog;
info.Y = Ybest;
info.fit_ = fit;
info.prevfits = prevfits;

end

%% ---------------- Helper functions ----------------
function [fieldNames, fieldSizes, fieldLengths, total] = extract_fieldinfo(S)
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

function v = struct2vec(S, fieldNames, fieldLengths)
numf = numel(fieldNames);
v = zeros(sum(fieldLengths),1);
idx = 1;
for k = 1:numf
    f = fieldNames{k};
    fld = S.(f);
    len = fieldLengths(k);
    if len == 0, continue; end
    v(idx:idx+len-1) = fld(:);
    idx = idx + len;
end
end

function S = vec2struct(v, fieldNames, fieldSizes, fieldLengths)
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
