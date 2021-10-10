function emissout=hexcav_2mat(epsc,lamc,epss,lams,h_d,l_d)
% inputs cavity emissivity and spectra, surface emissivity and spectra,
% depth of cavity, and center to center spacing, and returns a spectra with
% the same sampling as the surface spectra.
%% hexagonal stack, cylindrical cavity
D  =1;
R  =D/2;
A_2=pi*R^2;
A_1=A_2+2*pi*R*h_d;
A23=3*l_d^2/2/sqrt(3);
A_3=A23-A_2;
emissout=epss+NaN;
for ii=1:length(epss)
    pp=find(lamc>lams(ii),1);
    if pp<2
        epc=epsc(1)*lams(ii)/lamc(1);
    elseif isempty(pp)
        epc=min(epsc);
    else
        epc=epsc(pp-1)+(epsc(pp)-epsc(pp-1))*...
            (lams(ii)-lamc(pp-1))/(lamc(pp)-lamc(pp-1));
    end
    if epc==0
        ep1=1e10;
    else
        ep1=1/epc;
    end
    emissout(ii)=(A_2/A23)*1/(1+A_2/A_1*(ep1-1))+epss(ii)*A_3/A23;
end

end