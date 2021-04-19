# C16.jl returns the view factor between two adjacent rectangles, of
# shared length h, width of surface 2 a, width of surface 1 b, and angle
# between them phi. From UT Austin, C-16, thermalradiation.net

using SymPy
function C16(a,b,h,phi)
A=a/h;
B=b/h;
C=A^2+B^2-2*A*B*cos(phi);
D=sqrt(1+A^2*sin(phi)^2);
@vars x

# then, breaking up the giant term into several smaller ones:
t21=A*B*sin(phi)+(pi/2-phi)*(A^2+B^2)+
    B^2*atan((A-B*cos(phi))/(B*sin(phi)))+
    A^2*atan((B-A*cos(phi))/(A*sin(phi)));
t22=(2/sin(phi)^2-1)*log((1+A^2)*(1+B^2)/(1+C))+
    B^2*log((B^2*(1+C))/(C*(1+B^2)))+
    A^2*log((A^2*(1+A^2)^(cos(2*phi)))/(C*(1+C)^(cos(2*phi))));
t23=1/pi*atan(1/B)+A/pi/B*atan(1/A)-sqrt(C)/pi/B*atan(1/sqrt(C));
t24=atan(A*cos(phi)/D)+atan((B-A*cos(phi))/D);
t2h(x)=sqrt(1+x^2*sin(phi)^2);
t25(x)=t2h(x)*(atan(x*cos(phi)/t2h(x))+atan((A-x*cos(phi))/t2h(x)));
F12=-sin(2*phi)/4/pi/B*t21+
     sin(phi)^2/4/pi/B*t22+
                       t23+
     sin(phi)*sin(2*phi)/2/pi/B*A*D*t24+
     cos(phi)/pi/B*N(integrate(t25(x),(x,0,B)));
end