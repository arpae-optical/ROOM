% C11.m returns the view factor between two directly opposed, parallel
% rectangles, of width a, length b, and distance between them c. 
% From UT Austin, C-11, thermalradiation.net
function f12 = C11(a,b,c)

x = a/c;
y = b/c;

f12 = 2/pi/x/y*(1/2*log((1+x^2)*(1+y^2)/(1+x^2+y^2))+...
                x*sqrt(1+y^2)*atan(x/sqrt(1+y^2))+...
                y*sqrt(1+x^2)*atan(y/sqrt(1+x^2))-...
                x*atan(x)-y*atan(y));

end