% Analytical View Factor from Modest, M. "Radiative Heat Transfer", 3rd Ed.
% Appendix D. 39: Two finite rectangles of same length, with one common 
% edge, and at an angle of 90* to each other.
% l: length of shared edge
% w: width of area 1
% h: height of area 2

function F12=D39(l,w,h)
H=h/l;
W=w/l;
F12=1/pi/W*(W*atan(1/W)+...
            H*atan(1/H)-...
            sqrt(H^2+W^2)*atan(1/sqrt(H^2+W^2))+...
            1/4*log((1+W^2)*(1+H^2)/(1+W^2+H^2)*...
           ((W^2*(1+W^2+H^2)/(1+W^2)/(W^2+H^2)))^(W^2)*...
            (H^2*(1+H^2+W^2)/(1+H^2)/(H^2+W^2))^(H^2)));
end