function emissout=hexcavity(epsc,lamc,epss,lams,h_d,l_d)
% inputs cavity emissivity and
%% hexagonal stack, cylindrical cavity
D  =1;
R  =D/2;
A_2=pi*R^2;
A_1=A_2+2*pi*R*h_d;
F11=1-A_2/A_1;
A23=3*l_d^2/2/sqrt(3);
A_3=A23-A_2;
emissout=epss+NaN;
epsc=interp1(epsc,lamc,lams,'makima');
for ii=1:length(epss)
    epp=epsc(ii);
    if epp==0
        emissout(ii)=0;
    else
        ep1=1/epp;
        emissout(ii)=(A_1/A23)*(1-F11)/(ep1-F11*(ep1-1))+epss(ii)*A_3/A23;
    end
end

end