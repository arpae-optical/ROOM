% for triangular walls with side length l, and height h:
%  2/\2
% 2/BT\2 W
%   11
% let us find all view factors.
function epa=triwall(l_h,eps)
epa=eps.*0;
h  =1;
l  =l_h*h;
s  =sqrt(3/4)*l;
AW =3*s*h;
AB =1/2*s*l;

F12=D6plus(l,l,h,pi/3);

FWW=2*F12;
FWB=(1-FWW)/2;
FBW=AW*FWB/AB;
for ii=1:length(eps)
    if eps(ii)==0
        epa(ii)=0;
    else
        epa(ii)=(1-FBW/2)/(1/eps(ii)-FBW/2*(1/eps(ii)-1));
    end
end