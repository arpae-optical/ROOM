function FAB = DNS(gsource,gtarget)
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


% FAB = 0;
% for jj=1:length(gsource.A)
%     dF = 0;
%     sxyz = gsource.xyz(jj,:);
%     sn   = gsource.n(jj,:);
%     for ii=1:length(gtarget.A)
%         txyz = gtarget.xyz(ii,:);
%         tn   = gtarget.n(ii,:);
%         S  = txyz - sxyz;
%         dF = dF + dFdA(S,sn,tn)*gtarget.A(ii);
%     end
%     FAB = FAB + gsource.A(jj)*dF/pi;
% end
% FAB = FAB/sum(gsource.A);
% end
% 
% function result=dFdA(S,sn,tn)
% p1=dot(S,sn);
% p2=dot(-S,tn);
% result=p1*p2/norm(S)^4;
% end