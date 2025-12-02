function [Ybest, fbest, info,gbestlog] = RAQPSO(problem, options)
dim = problem.M.dim();
[n,p] = size(problem.M.rand());
localdefaults = struct('popsize',4*floor(log(dim)),...
    'maxcostevals',nn*p*10000);
options = mergeOptions(mergeOptions(getGlobalDefaults(), localdefaults),options);
iter = 1;
%initialtion
pop = cell(popsize,1);
fit=zeros(1,popsize);
pbest_pop = cell(1,popsize);%pbest
pbestvalue = zeros(1,popsize);
for i = 1:popsize
    pop{i} = problem.M.rand();
    fit(i) = problem.cost(pop{i});
end
pbest_pop =pop;
pbestvalue= fit;
[~,index] = sort(pbestvalue);
Ybest = pop{index(1)};
FES(1) = popsize;
fbest=pbestvalue(index(1));
gbest(iter) =fbest;
maxiter = n*p*10000/popsize;
for   i=1:maxiter
    rtt(i)=1-(i/maxiter);
end
gbestlog(iter)=fbest;

while FES(iter) <=maxcostevals
    pop_pre = pop;
    prevfits = fit;
    if iter<e*maxcostevals
        if FES(iter)<e*maxcostevals
            tt=randi(popsize);
            mbest=findmean(problem,tt,pbest_pop,pbestvalue,popsize,rtt(iter));
        else
            mbest=findmean1(problem,Ybest,pbest_pop,pbestvalue,popsize,rtt(iter));
        end
    end
        
        
        %% update pos
        for i = 1:popsize
            
            pv=problem.M.proj(pbest_pop{i},Ybest);
            tmp1 = problem.M.retr(pbest_pop{i},pv,0.1*rand());
            tmpvec=normrnd(Ybest,10.*abs(pv));
            tmpvec1=problem.M.proj(tmp1,tmpvec);
            p_temp= problem.M.retr(tmp1,tmpvec1,0.1*rand());
            p_temp2 = problem.M.transp(pop{i},p_temp,(problem.M.proj(pop{i},mbest)));
            if rand() < 0.5
                pop{i} = problem.M.retr(p_temp,beta.*log(1./rand(n,p)).*p_temp2, rtt(iter));
            else
                pop{i} = problem.M.retr(p_temp,-beta.*log(1./rand(n,p)).*p_temp2, rtt(iter));
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
        
        fprintf('[#%d FEs:%d] fit = %+.16e rt = %+.16e\n',iter,FES(iter),fbest,rtt(iter));
        FES(iter+1) = FES(iter)+ popsize;
        iter = iter + 1;
        gbestlog(iter)=fbest;
end
info = struct();
info.FES = FES;
info.best = gbest;
info.Y = Ybest;
info.fit_ = fit;
info.prevfits = prevfits;
