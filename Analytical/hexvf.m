% for a hex with center to center spacing l, and height h, side length s:
%  2/\3
% 1|  |4
%  2\/3
% let us find all view factors. We will start by declaring necessary
% variables.
function epa=hexwall(l_h,eps)
h  =1;
l  =l_h*h;
s  =l/sqrt(3);
p12=2*pi/3;
% p13=pi/3;
A=s/h;
B=s/h;
C=A^2+B^2-2*A*B*cos(p12);
D=sqrt(1+A^2*sin(p12)^2);
syms xsi

% then, breaking up the giant term into several smaller ones:
t21=A*B*sin(p12)+(pi/2-p12)*(A^2+B^2)+...
    B^2*atan((A-B*cos(p12))/(B*sin(p12)))+...
    A^2*atan((B-A*cos(p12))/(A*sin(p12)));
t22=(2/sin(p12)^2-1)*log((1+A^2)*(1+B^2)/(1+C))+...
    B^2*log((B^2*(1+C))/(C*(1+B^2)))+...
    A^2*log((A^2*(1+A^2)^(cos(2*p12)))/(C*(1+C)^(cos(2*p12))));
t23=1/pi*atan(1/B)+A/pi/B*atan(1/A)-sqrt(C)/pi/B*atan(1/sqrt(C));
t24=atan(A*cos(p12)/D)+atan((B-A*cos(p12))/D);
t2h=sqrt(1+xsi^2*sin(p12)^2);
t25=t2h*(atan(xsi*cos(p12)/t2h)+atan((A-xsi*cos(p12))/t2h));
F12=-sin(2*p12)/4/pi/B*t21+...
    sin(p12)^2/4/pi/B*t22+...
    t23+...
    sin(p12)*sin(2*p12)/2/pi/B*A*D*t24+...
    cos(p12)/pi/B*int(t25,xsi,0,B);
F12=double(F12);
for ii=1:length(eps)
        if eps(ii)==0
            epa(ii,jj)=0;
        else
            epa(ii,jj)=(1-f12(jj)/2)/(1/eps(ii)-f12(jj)/2*(1/eps(ii)-1));
        end
    end