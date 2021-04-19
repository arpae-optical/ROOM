function emissout=hexsphere(emissin,l_d)
%% hexagonal stack, sphere on surface
% view factors from 
% Feingold, A. and Gupta, K.G., 1970,
% "New analytical approach to the evaluation of configuration factors in radiation 
% from spheres and infinitely long cylinders,"
% J. Heat Transfer, vol. 92, no. 1, pp. 69-76, February.
D  =1; %diameter of sphere
R  =D/2; %radius of sphere    /\
L  =l_d*D; %width of hexagon |  | distance between |  |
%                             \/
s=sqrt(L^2/3); %length of one side of hexagon |
% 
% a=R+eps; %sphere resting on surface
% b1=L/2; %base of triangle half of hexagon width
% b2=b1+eps; %right triangle, not scalene
% b3=s; %hypotenuse of triangle is the length of one side |
% 
% B1=b1/a;
% B2=b2/b1;
% B3=b3/b1;

As=4*pi*R^2; %surface area of sphere
Ah=3/2*s*L; %surface area of hexagon (side times base over two times six)
A1=As+Ah; %surface area of solid
% Av=Ah+s*D*6; %virtual surface: hexagon (A2) plus six faces (s*D*6)

% Fst=1/(4*pi)*(acos(1/B3)-acos(1/B2))-1/(8*pi)*...
% (asin(((1-B1^2)*B3^2-2)/((1+B1^2)*B3^2))-...
%  asin(((1-B1^2)*B2^2-2)/((1+B1^2)*B2^2)));
%Fst is amount of radiation leaving sphere, 
%and arriving at one right triangle of hexagon-
%the factor we want is that value times twelve
% Fsh=Fst*12;
%Fsh is amount of radiation leaving sphere
%and arriving at hexagonal surface
% Fsv=1-Fsh;
%Fsv is the radiation leaving sphere and arriving
%all virtual surfaces
% Fvs=Fsv*As/Av;
%Fvs is the radiation leaving virtual surface
%and arriving at sphere
% Fhs=Fsh*As/Ah;
%Fhs is radiation leaving hexagon and arriving
%at sphere
% Fhv=1-Fhs;
%hex->virtual
% Fvh=Fhv*Ah/Av;
%virtual->hex
% Fv1=Fvh+Fvs;
%virtual->hex+sphere
% F1v=Fv1*Av/A1;
%hex+sphere->virtual
% F11=1-F1v;
%hex+sphere->hex+sphere
emissout=0.*emissin;
A_1=A1; %emitting surface
A_2=Ah; %apparent area
% K=1-F11-A_2/A_1;
for iii=1:length(emissin)
    if emissin(iii)==0
        emissout(iii)=0;
    else
        epp=emissin(iii);
        ep1=1/epp;
        emissout(iii)=1/(1+A_2/A_1*(ep1-1));
    end
end

end