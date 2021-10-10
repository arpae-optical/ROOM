function emissout=hexcolumn(emissin,h_d,l_d)
%% hexagonal stack, cylindrical column
D  =1; %diameter of cylinder
R  =D/2; %radius
L  =l_d*D; %horizontal spacing between centers of cylinder
H  =h_d*D; %height of cylinders
s=L/sqrt(3); %length of one side of hex
A_H=s*L*3/2; %area of hexagon
A_C=pi*R^2; %area of circular cylinder cap
A_V=pi*D*H; %area of vertical wall of cylinder
emissout=0.*emissin;
A_1=A_H-A_C+A_V; %emitting surface
A_2=A_H-A_C; %apparent area
for iii=1:length(emissin)
    epp=emissin(iii);
    ep1=1/epp;
    emissout(iii)=(A_2/A_H)*1/(1+A_2/A_1*(ep1-1))+epp*A_C/A_H;
end

end