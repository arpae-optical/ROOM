function emissout=hexcone(emissin,h_d,l_d)
%% hexagonal stack, cylindrical column
D  =1; %diameter of cylinder
R  =D/2; %radius
L  =l_d*D; %horizontal spacing between centers of cylinder
H  =h_d*D; %height of cylinders
l  =sqrt(R^2+H^2); %slant height of cone
s=L/sqrt(3); %length of one side of hex
A_H=s*L*3/2; %area of hexagon
A_C=pi*R^2; %area of circular base
A_V=pi*R*l; %area of outside of cone
emissout=0.*emissin;
A_1=A_H-A_C+A_V; %emitting surface
A_2=A_H; %apparent area
F11=1-A_2/A_1; %emitting surface self view factor
for iii=1:length(emissin)
    epp=emissin(iii);
    ep1=1/epp;
    emissout(iii)=(A_1/A_2)*(1-F11)/(ep1-F11*(ep1-1));
end

end