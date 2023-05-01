n=80;
iter2=500;
iter1=50;
startsize2=16;
globit=100;
lyaps=zeros(4,globit);
deads=zeros(1,globit);
solvs=zeros(1,globit);
times=zeros(1,globit);
dists=zeros(1,globit);
halls=zeros(4,globit);
rng(1);

for l=1:globit
cornerseed=randi([0,1],startsize2,startsize2);
N=max(2*n,n+100);
Domain=zeros(N,N);

startpos=[(N)/2-floor(startsize2/2),(N)/2-floor(startsize2/2)];

Domain(startpos(1):startpos(1)+startsize2-1,startpos(2):startpos(2)+startsize2-1)=cornerseed;

Domain2=Domain;
Domain2(startpos(1),startpos(2))=mod(Domain2(startpos(1),startpos(2))+1,2);
Domain5=Domain;
Domain5(startpos(1)+startsize2-1,startpos(2))=mod(Domain5(startpos(1)+startsize2-1,startpos(2))+1,2);
Domain4=Domain;
Domain4(startpos(1)+startsize2-1,startpos(2)+startsize2-1)=mod(Domain4(startpos(1)+startsize2-1,startpos(2)+startsize2-1)+1,2);
Domain3=Domain;
Domain3(startpos(1),startpos(2)+startsize2-1)=mod(Domain3(startpos(1),startpos(2)+startsize2-1)+1,2);

diff=zeros(4,iter1+1);
diff(:,1)=ones(4,1);

for i = 1:iter1
    newDomain=zeros(N,N);
    newDomain2=zeros(N,N);
    newDomain3=zeros(N,N);
    newDomain4=zeros(N,N);
    newDomain5=zeros(N,N);
for j=2:N-1
    for k=2:N-1
        if Domain(j,k)==1
            if sum(Domain(j-1:j+1,k-1:k+1),'all')>=2 && sum(Domain(j-1:j+1,k-1:k+1),'all')<=6
            newDomain(j,k)=1;
            else 
            newDomain(j,k)=0;
            end
        else
            if sum(Domain(j-1:j+1,k-1:k+1),'all')==3
            newDomain(j,k)=1;
            else
                newDomain(j,k)=0;
            end
        end
        if diff(1,i)>0
        if Domain2(j,k)==1
            if sum(Domain2(j-1:j+1,k-1:k+1),'all')>=2 && sum(Domain2(j-1:j+1,k-1:k+1),'all')<=6
            newDomain2(j,k)=1;
            else 
            newDomain2(j,k)=0;
            end
        else
            if sum(Domain2(j-1:j+1,k-1:k+1),'all')==3
            newDomain2(j,k)=1;
            else
                newDomain2(j,k)=0;
            end
        end
        end
        if diff(2,i)>0
        if Domain3(j,k)==1
            if sum(Domain3(j-1:j+1,k-1:k+1),'all')>=2 && sum(Domain3(j-1:j+1,k-1:k+1),'all')<=6
            newDomain3(j,k)=1;
            else 
            newDomain3(j,k)=0;
            end
        else
            if sum(Domain3(j-1:j+1,k-1:k+1),'all')==3
            newDomain3(j,k)=1;
            else
                newDomain3(j,k)=0;
            end
        end
        end
        if diff(3,i)>0
        if Domain4(j,k)==1
            if sum(Domain4(j-1:j+1,k-1:k+1),'all')>=2 && sum(Domain4(j-1:j+1,k-1:k+1),'all')<=6
            newDomain4(j,k)=1;
            else 
            newDomain4(j,k)=0;
            end
        else
            if sum(Domain4(j-1:j+1,k-1:k+1),'all')==3
            newDomain4(j,k)=1;
            else
                newDomain4(j,k)=0;
            end
        end
        end
                if diff(4,i)>0
        if Domain5(j,k)==1
            if sum(Domain5(j-1:j+1,k-1:k+1),'all')>=2 && sum(Domain5(j-1:j+1,k-1:k+1),'all')<=6
            newDomain5(j,k)=1;
            else 
            newDomain5(j,k)=0;
            end
        else
            if sum(Domain5(j-1:j+1,k-1:k+1),'all')==3
            newDomain5(j,k)=1;
            else
                newDomain5(j,k)=0;
            end
        end
        end
    end
end
    Domain=newDomain;
    Domain2=newDomain;
    Domain3=rot90(newDomain);
    Domain4=rot90(newDomain,2);
    Domain5=rot90(newDomain,3);
    if diff(1,i)>0&&sum(abs(newDomain-newDomain2),'all')>0
    diff(1,i+1)=sum(abs(newDomain-newDomain2),'all');
    [row,col]=find(~newDomain2==newDomain);
    Domain2(row(1),col(1))=mod(Domain2(row(1),col(1))+1,2);
    end
    if diff(2,i)>0&&sum(abs(newDomain-newDomain3),'all')>0
    diff(2,i+1)=sum(abs(newDomain-newDomain3),'all');
    [row,col]=find(~rot90(newDomain3)==rot90(newDomain));
    Domain3(row(1),col(1))=mod(Domain3(row(1),col(1))+1,2);
    Domain3=rot90(Domain3,-1);
    end
    if diff(3,i)>0&&sum(abs(newDomain-newDomain4),'all')>0
         diff(3,i+1)=sum(abs(newDomain-newDomain4),'all');
    [row,col]=find(~rot90(newDomain4,2)==rot90(newDomain,2));
    Domain4(row(1),col(1))=mod(Domain4(row(1),col(1))+1,2);
    Domain4=rot90(Domain4,-2);
    end
    if diff(4,i)>0&&sum(abs(newDomain-newDomain5),'all')>0
         diff(4,i+1)=sum(abs(newDomain-newDomain5),'all');
    [row,col]=find(~rot90(newDomain5,3)==rot90(newDomain,3));
    Domain5(row(1),col(1))=mod(Domain5(row(1),col(1))+1,2);
    Domain5=rot90(Domain5,-3);
    end
end

for d=1:4
if min(diff(d,:))==0
    lyaps(d,l)=0;
else
    lyaps(d,l)=sum(log2(diff(d,:)))/length(diff(d,:));
end
end

domain=zeros(n,n);
domain(2:2+startsize2-1,2:2+startsize2-1)=cornerseed;
domain(end-startsize2:end-1,end-startsize2:end-1)=cornerseed;
domain(2:2+startsize2-1,end-startsize2:end-1)=cornerseed;
domain(end-startsize2:end-1,2:2+startsize2-1)=cornerseed;

olddomain=domain;
for i = 1:iter2
newdomain=zeros(n,n);
for j=2:n-1
    for k=2:n-1
        if domain(j,k)==1
            if sum(domain(j-1:j+1,k-1:k+1),'all')>=2 && sum(domain(j-1:j+1,k-1:k+1),'all')<=6
            newdomain(j,k)=1;
            else 
            newdomain(j,k)=0;
            end
        else
            if sum(domain(j-1:j+1,k-1:k+1),'all')==3
            newdomain(j,k)=1;
            else
                newdomain(j,k)=0;
            end
        end
    end
end
if isequal(domain,newdomain) || isequal(olddomain,newdomain)
    break
else
    olddomain=domain;
    domain=newdomain;
end
end
%% 
dists(l)=sum(domain,'all')/n^2;
for i=2:n-1
    for j=2:n-1
        if domain(i,j)==1
            if sum(domain(i-1,j)+domain(i+1,j)+domain(i,j-1)+domain(i,j+1))==1
                halls(1,l)=halls(1,l)+1;
            elseif sum(domain(i-1,j)+domain(i+1,j)+domain(i,j-1)+domain(i,j+1))==2
                halls(2,l)=halls(2,l)+1;
            elseif sum(domain(i-1,j)+domain(i+1,j)+domain(i,j-1)+domain(i,j+1))==3
                halls(3,l)=halls(3,l)+1; 
            elseif sum(domain(i-1,j)+domain(i+1,j)+domain(i,j-1)+domain(i,j+1))==4
                halls(4,l)=halls(4,l)+1;
            end
        end
    end
end

domain2=domain;
domain2(2:3,2:3)=ones(2,2);
domain2(end-2:end-1,end-2:end-1)=ones(2,2);

Domain3=domain2;
olddeads=zeros(size(Domain3));
deadends=zeros(size(Domain3));
for k=1:100
for i=2:n-1
    for j=2:n-1
        if Domain3(i-1,j)+Domain3(i+1,j)+Domain3(i,j-1)+Domain3(i,j+1)<=1 && Domain3(i,j)==1
            deadends(i,j)=1;
            Domain3(i,j)=0;
        end
    end
end
if olddeads==deadends
    break
end
olddeads=deadends;
end
deads(l)=sum(olddeads,'all');

Domain5=Domain3;

paths=[2,2];
orients=[1,0];
solution=zeros(size(Domain5));
solution(2,2)=1;

solved=-1;

while (solved==-1)
    if isequal(orients,[1,0])
        if Domain5(paths(end,1)-1,paths(end,2))==1
            solution(paths(end,1)-1,paths(end,2))=1;
            paths=[paths;paths(end,1)-1,paths(end,2)];
            orients=[0,1];
        elseif Domain5(paths(end,1),paths(end,2)+1)==1
             solution(paths(end,1),paths(end,2)+1)=1;
            paths=[paths;paths(end,1),paths(end,2)+1];
            orients=[1,0];       
        else 
            orients=[0,-1];
        end
    elseif isequal(orients,[0,1])
        if Domain5(paths(end,1),paths(end,2)-1)==1
            solution(paths(end,1),paths(end,2)-1)=1;
            paths=[paths;paths(end,1),paths(end,2)-1];
            orients=[-1,0];
        elseif Domain5(paths(end,1)-1,paths(end,2))==1
            solution(paths(end,1)-1,paths(end,2))=1;
            paths=[paths;paths(end,1)-1,paths(end,2)];
            orients=[0,1];
        else 
            orients=[1,0];
        end
     elseif isequal(orients,[-1,0])
        if Domain5(paths(end,1)+1,paths(end,2))==1
            solution(paths(end,1)+1,paths(end,2))=1;
            paths=[paths;paths(end,1)+1,paths(end,2)];
            orients=[0,-1];
        elseif Domain5(paths(end,1),paths(end,2)-1)==1
            solution(paths(end,1),paths(end,2)-1)=1;
            paths=[paths;paths(end,1),paths(end,2)-1];
            orients=[-1,0];
        else 
            orients=[0,1];
        end 
     elseif isequal(orients,[0,-1])
        if Domain5(paths(end,1),paths(end,2)+1)==1
             solution(paths(end,1),paths(end,2)+1)=1;
            paths=[paths;paths(end,1),paths(end,2)+1];
            orients=[1,0];
        elseif Domain5(paths(end,1)+1,paths(end,2))==1
             solution(paths(end,1)+1,paths(end,2))=1;
            paths=[paths;paths(end,1)+1,paths(end,2)];
            orients=[0,-1];
        else 
            orients=[-1,0];
        end
    end
    if isequal(paths(end,:),[2,2]) && isequal(orients,[1,0])
        solved=0;
    end
    if isequal(paths(end,:),[(n-1),(n-1)])
        solved=1;
    end
end
    solvs(l)=solved;
    times(l)=length(paths);
end

figure();
Corr=[dists',sum(lyaps)',max(lyaps)',solvs',times',deads',halls(3,:)'+halls(4,:)'];
imagesc(corrcoef(Corr));
xticks([1 2 3 4 5 6 7]);
xticklabels({'Dists','Sum L','Max L','Solved','Path','Dead','Branch'});
yticks([1 2 3 4 5 6 7]);
yticklabels({'Dists','Sum L','Max L','Solved','Path','Dead','Branch'});
title('100 seeds, 16 vs 80')
colorbar();

figure();
CorrSolv1=[dists(solvs==1)',sum(lyaps(:,solvs==1))',max(lyaps(:,solvs==1))',times(solvs==1)',deads(solvs==1)',halls(3,solvs==1)'+halls(4,solvs==1)'];
imagesc(corrcoef(CorrSolv1));
xticks([1 2 3 4 5 6]);
xticklabels({'Dists','Sum L','Max L','Path','Dead','Branch'});
yticks([1 2 3 4 5 6]);
yticklabels({'Dists','Sum L','Max L','Path','Dead','Branch'});
title('Solved, 16 vs 80')
colorbar();

figure();
CorrSolv0=[dists(solvs==0)',sum(lyaps(:,solvs==0))',max(lyaps(:,solvs==0))',times(solvs==0)',deads(solvs==0)',halls(3,solvs==0)'+halls(4,solvs==0)'];
imagesc(corrcoef(CorrSolv0));
xticks([1 2 3 4 5 6]);
xticklabels({'Dists','Sum L','Max L','Path','Dead','Branch'});
yticks([1 2 3 4 5 6]);
yticklabels({'Dists','Sum L','Max L','Path','Dead','Branch'});
title('Not Solved, 16 vs 80')
colorbar();
