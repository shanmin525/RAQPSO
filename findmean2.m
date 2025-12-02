function mbest=findmean2(problem,pbest_pop,pbest,popsize)
[~,index] = sort(pbest);
m1=pbest_pop{index(1)};
for i=2:popsize
    v=problem.M.proj(m1,pbest_pop{index(i)});
    m1=problem.M.exp(m1,v);
end
mbest=m1;
         
     
   