% C121.m returns the view factor between a sphere distance d away from a 
% rectangle of width w and length l, centered on and perpendicular to the 
% diameter of the sphere. From UT Austin, C-121, thermalradiation.net
function F12=C121(d,l,w)
B=((l/2)/d);
C=((w/2)/d);
n1=(2*B^2-(1-B^2)*(B^2+C^2))/((1+B^2)*(B^2+C^2));
n2=(2*C^2-(1-C^2)*(B^2+C^2))/((1+C^2)*(B^2+C^2));
F12 = 1/2/pi*(asin(n1)+asin(n2));