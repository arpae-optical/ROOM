function [FWW,t]=FhexWW(L_H)
tic
H   = 1;
L   = L_H*H;
s   = L/sqrt(3);

F12 = C16(1*s,1*s,H,2*pi/3);

A1=s*H;
A1m=s*H;
A1p=2*s*H;
% A3=A1;
A3m=A1m;
A3p=A1p;

F3p1p=C16(2*s,2*s,H,pi/3);
F3p1m=C16(1*s,2*s,H,pi/3);
F3p1=F3p1p-F3p1m;
F13p=A3p*F3p1/A1;

F3m1m=C16(1*s,1*s,H,pi/3);
F3m1p=C16(2*s,1*s,H,pi/3);
F3m1=F3m1p-F3m1m;
F13m=A3m*F3m1/A1;
F13 = F13p - F13m;

F14 = D38(H,s,L);

FWW = 2*F12+2*F13+F14;
t=toc;
end