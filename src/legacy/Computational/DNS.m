function FAB = DNS(gsource,gtarget,varargin)
if nargin==4
    param = varargin{2};
else
    param = 2;
end
if nargin==2
    sox = gsource.xyz(:,1);
    soy = gsource.xyz(:,2);
    soz = gsource.xyz(:,3);
    snx = gsource.n(:,1);
    sny = gsource.n(:,2);
    snz = gsource.n(:,3);
    sA  = gsource.A;
    
    tox = gtarget.xyz(:,1)';
    toy = gtarget.xyz(:,2)';
    toz = gtarget.xyz(:,3)';
    tnx = gtarget.n(:,1)';
    tny = gtarget.n(:,2)';
    tnz = gtarget.n(:,3)';
    tA  = gtarget.A';
    
    Sx  = tox-sox;
    Sy  = toy-soy;
    Sz  = toz-soz;
    Sm2 = Sx.^2+Sy.^2+Sz.^2;
    p1  = Sx.*snx+Sy.*sny+Sz.*snz;
    p2  =-Sx.*tnx-Sy.*tny-Sz.*tnz;
    dF  = p1.*p2.*sA.*tA./Sm2.^2;
    FAB = sum(sum(dF))/sum(sA)/pi;
else
    sox = gsource.xyz(:,1);
    soy = gsource.xyz(:,2);
    soz = gsource.xyz(:,3);
    snx = gsource.n(:,1);
    sny = gsource.n(:,2);
    snz = gsource.n(:,3);
    sA  = gsource.A;
    ns  = length(sA);
    
    tox = gtarget.xyz(:,1);
    toy = gtarget.xyz(:,2);
    toz = gtarget.xyz(:,3);
    tnx = gtarget.n(:,1);
    tny = gtarget.n(:,2);
    tnz = gtarget.n(:,3);
    tA  = gtarget.A;
    nt  = length(tA);
    
    pclude=varargin{1};
    
    FAB = 0;
    for jj=1:length(gsource.A)        
        vis  = HPR([tox,toy,toz;pclude],[sox(jj),soy(jj),soz(jj)],param);
        vis  = vis<nt+1;
        Sx  = tox(vis)-sox(jj);
        Sy  = toy(vis)-soy(jj);
        Sz  = toz(vis)-soz(jj);
        Sm2 = Sx.^2+Sy.^2+Sz.^2;
        p1  = Sx.*snx(jj)+Sy.*sny(jj)+Sz.*snz(jj);
        p2  =-Sx.*tnx(vis)-Sy.*tny(vis)-Sz.*tnz(vis);
        dF  = p1.*p2.*sA(jj).*tA(vis)./Sm2.^2;
        FAB = FAB+sum(dF);
    end
    FAB = FAB/sum(sA)/pi;
end
end

function bugplot(pclude,tox,toy,toz,sox,soy,soz,jj,vis)
scatter3(pclude(:,1),pclude(:,2),pclude(:,3),1,[0,0,0])
hold on
scatter3(tox,toy,toz,1,[0,0,0])
scatter3(tox(vis),toy(vis),toz(vis),1,[1,0,0])
scatter3(sox(jj),soy(jj),soz(jj),1,[0,1,0])
end