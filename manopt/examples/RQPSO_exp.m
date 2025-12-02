function [Ybest, fbest, info,time] = RQPSO_exp(problem, Y, options,beta,e)
dim = problem.M.dim();
[n,p] = size(problem.M.rand());
localdefaults = struct('minsigma',1e-6,...
    'maxiter',n*p*500,...
    'damping',1,...
    'popsize',40);
options = mergeOptions(mergeOptions(getGlobalDefaults(), localdefaults),options);
maxiter = n*p*500;
% popsize=30;
% sigma_ = options.sigma;
popsize = options.popsize;
maxcostevals = n*p*500*popsize;
iter = 1;
pop = cell(popsize,1);

for i = 1:popsize
    pop{i} = problem.M.rand();%随机生成粒子位置
end
pbest_pop = cell(1,popsize);%pbest
% pvalue = zeros(popsize,maxiter);
% pvalue = zeros(popsize,1);
pbestvalue = zeros(1,popsize);
gbest = ones(maxiter,1)*(inf);
fit=zeros(1,popsize);
for i = 1:popsize
    fit(i) = problem.cost(pop{i});
end
pbest_pop =pop;
pbestvalue= fit;
[~,index] = sort(pbestvalue);
Ybest = pop{index(1)};
FES(1) = popsize;
fbest=pbestvalue(index(1));
gbest(iter) =fbest;
time=0;
for i=1:maxiter
     rtt(i)=0.5*exp(-15*(i/maxiter));
%           rtt(i)=exp(-15*(i/maxiter));%SC
end
%       rtt(i)=0.1*exp(-10*(i/(maxiter)));
% %         rtt(i)=rt_start-(rt_start-rt_end)/(maxiter/10)*i;
% %       rtt(i)=rt_start*(rt_end/rt_start)^(i/(maxiter/10));
%       
%     
%%
% while iter <= maxiter
 while iter <= maxiter && FES(iter) <=maxcostevals %&& sigma_ > options.minsigma
    pop_pre = pop;
    prevfits = fit;
%% mbest
%     tic
    tt=randi(popsize);
    mbest = pbest_pop{tt};
%     toc
%     time(iter)= toc;
%     for i=1:popsize
%         rtt(i,iter)=rt_end+(rt_start-rt_end)*exp((fit(i)-fmin)/(fmax-fmin))^0.35;
% %         rtt(i,iter)=rt_end+(rt_start-rt_end)*(1/(1+exp(-abs((fit(i)-fmin))/fmax)));
% % %         rtt(i,iter)=rt_end+(rt_start-rt_end)*fit(i)/sum(fit);
% %          rtt(i,iter)=rt_start-(rt_start-rt_end)*(1/(1+exp(-abs((fit(i)-fmin)/(fmax-fmin)))));
%    end
%     rtt(iter)=rt_start-(rt_start-rt_end)/maxiter*iter;
%% update pos
    for i = 1:popsize
        pv=problem.M.proj(pbest_pop{i},Ybest);
        tmp1 = problem.M.retr(pbest_pop{i},pv,0.1*rand());
%         tmpvec=Ybest+randn(n,p).*e.*sqrt(pv);
%         tmpvec=Ybest+randn(n,p).*e.*pv;
        tmpvec=Ybest+randn(n,p).*e.*pbest_pop{i};
        tmpvec=real(tmpvec);
        tmpvec1=problem.M.proj(tmp1,tmpvec);
        p_temp= problem.M.retr(tmp1,tmpvec1,0.1*rand());
        p_temp2 = problem.M.transp(pop{i},p_temp,(problem.M.proj(pop{i},mbest)));
        if rand() < 0.5
             pop{i} = problem.M.retr(p_temp,beta.*log(1./rand(n,p)).*abs(p_temp2), rtt(iter));
        else
            pop{i} = problem.M.retr(p_temp,beta.*log(1./rand(n,p)).*abs(p_temp2), rtt(iter));
        end       
    end
 %% update pbest and gbest
    for i = 1:popsize
        fit(i) = problem.cost(pop{i});
        gbest(iter)=fbest;
        if fit(i) < pbestvalue(i)
            pbestvalue(i) = fit(i);
            pbest_pop{i} = pop{i};
        end
        if pbestvalue(i)<fbest
            Ybest=pbest_pop{i};
            fbest=pbestvalue(i);
            gbest(iter)=fbest;
        end  

    end
     pop_sum = [pop_pre,pop];
    [~,imax] = sort([prevfits,fit]);
    for i = 1:popsize
        pop{i} = pop_sum{imax(i)};
    end

%    fprintf('[#%d FEs:%d] fit = %+.16e rt = %+.16e\n',iter,FES(iter),fbest,rtt(iter));
    FES(iter+1) = FES(iter)+ popsize;
    iter = iter + 1;
     
 end
% gbest(iter+1:end) = [];
% FES(iter+1:end) = [];
info = struct();
info.FES = FES;
info.best = gbest;
% fbest = min(gbest);
info.Y = Ybest;
info.fit_ = fit;
info.prevfits = prevfits;
