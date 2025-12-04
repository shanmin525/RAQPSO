function mbest=findmean(problem,id,pbest_pop,pbest,popsize,rtt)
k=popsize;
w=cell(1,k);  
fit=zeros(1,k);
for i=1:k
    fit(i)=pbest(i);
    w{i}=problem.M.proj(pbest_pop{id},pbest_pop{i});
end
ww = max(fit) - fit + 0.0001;
ww = ww / sum(ww);
meanw = zeros(size(w{1}));
for i = 1:k
    meanw = meanw + ww(i) * w{i};
end
mbest=problem.M.retr(pbest_pop{id},meanw,rtt);
         
     
   
