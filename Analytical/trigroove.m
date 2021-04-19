function emissout=trigroove(emissin,H,L)
%% hexagonal stack, cylindrical cavity
D=1;
emissout=0.*emissin;
A_1=sqrt(4*H^2+D^2);
A_2=D;
A23=L;
A_3=A23-A_2;
F11=1-A_2/A_1;
for iii=1:length(emissin)
    epp=emissin(iii);
    ep1=1/epp;
    emissout(iii)=(A_1/A23)*(1-F11)/(ep1-F11*(ep1-1))+epp*A_3/A23;
end
end