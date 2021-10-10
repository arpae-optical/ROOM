% This function takes an emissivity spectra for the base and the sphere, as
% well as a L/D ratio and a "polling" set of wavelengths at which to make
% the calculations. It outputs an array of apparent emissivities at the
% wavelengths from l_out.

function ea=hexsphere2mat(e_base,l_base,e_sphe,l_sphe,L_D,l_out)
%% geometry definitions
D=1;          % We'll take D to be equal to one.
L=L_D*D;      % L, then, is L/D * D
s=sqrt(L^2/3);% The length of one side of the hexagon is s

A_sphe=pi*D^2; % area of the sphere
A_base=3/2*s*L; % area of both of the hexagons

%% view factor definitions
a=D/2;        % Distance from center of sphere to plane of triangles
b1=L/2;       % Distance from center of sphere to base of triangle
b2=b1;        % Length of edge of scalene triangle (same as b1)
b3=s;         % Length of hypotenuse = length of side of hexagon

B1=b1/a;      % Nondim param for calculating F
B2=b2/b1;
B3=b3/b1;

F_st=1/(4*pi)*(acos(1/B3)-acos(1/B2))-...
     1/(8*pi)*(asin(((1-B1^2)*B3^2-2)/((1+B1^2)*B3^2))-...
               asin(((1-B1^2)*B2^2-2)/((1+B1^2)*B2^2)));

%% emissivity stuff
e1= interp1(l_sphe,e_sphe,l_out,'makima',min(e_sphe));
e2= interp1(l_base,e_base,l_out,'makima',min(e_base));
ea= 0.*e1;

%% 1 2 3 4 definitions
A1=A_sphe;
A2=A_base;

F11=0;
F12=F_st*12;
F13=F12;
F14=1-F11-F12-F13;

F21=A1*F12/A2;
F22=0;
F23=(1-F12)/2;
F24=1-F21-F22-F23;

for ii=1:length(e1)
    num=e1(ii)*F14*(1+F23*(1-F21))+...
        e2(ii)*(F24-F21*F14+F23*F24-F23*F21*F14)+...
        e1(ii)*e2(ii)*(F21*F14+F12*F24-F14*F23+F23*F24+2*F23*F21*F14);
    dom=F24-F21*F14-F23*F24+F23*F21*F14+...
        e1(ii)*(F14+F21*F14-F12*F24-F23*F21*F14)+...
        e2(ii)*(F24+F23*F24-F21*F14-F23*F21*F14)+...
        e1(ii)*e2(ii)*(F12*F24+F21*F14+F23*F21*F14);
    ea(ii)=num/dom;
end