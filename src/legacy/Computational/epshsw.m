function ea = epshsw(e1,e2,dl,hl)
% This function returns the apparent emissivity for a hex-sphere-adwall
% configuration where the substrate and the sphere are two dissimilar
% materials. e1 is the emissivity of the substrate, e2 is the emissivity of
% the sphere. THIS FUNCTION ASSUMES A SINGULAR VALUE OF e1 AND e2, OR THAT
% THEY ARE THE SAME LENGTH. Any spectra interpolation/extrapolation to
% satisfy this assumption must be done outside of the function. 

if dl==0&&hl==0
    ea = e1;
elseif dl==0
% first call the HSW view factor function- it will return analytic view
% factors for the hex-to-wall configuration, in an array, vfr
    vfr=vfrhsw(dl,hl); 
    F14=vfr(4);

    ea = (1-F14/2)./(1./e1-F14/2*(1./e1-1));
elseif hl==0
% first call the HSW view factor function- it will return analytic view
% factors as necessary, and monte carlo simulation for view factors that
% require computational solutions, in an array, vfr
    vfr = vfrhsw(dl,hl); 
    F12 = vfr(2);
    F21 = vfr(5);
    F22 = vfr(6);
    F31 = vfr(9);
    F32 = vfr(10);
    n1 = -e1*(F31*(1-F22)+F32*F21);
    n2 = -e2*(F32+F31*F12);
    n12= e1.*e2*(F32*F21+F31*F12-F31*F22);

    d12 = (1-e1).*(1-e2)*F12*F21;
    d20 = (1-e2).*F22-1;

    ea=(n1+n2+n12)./(d12+d20);
else

% first call the HSW view factor function- it will return analytic view
% factors as necessary, and monte carlo simulation for view factors that
% require computational solutions, in an array, vfr
    vfr=vfrhsw(dl,hl);
    F12=vfr(2);
    F14=vfr(4);
    F21=vfr(5);
    F22=vfr(6);
    F24=vfr(8);
    F31=vfr(9);
    F32=vfr(10);
    F34=vfr(12);
    F41=vfr(13);
    F42=vfr(14);
    F44=vfr(16);

n0  = -e1.*e2*(F31*(1-F44)+F34*F41)-e2*(F32*(1-F44)+F34*F42);

n12 = (e1.*e2-e1)*(F31*(1-F22-F44-F24*F42+F22*F44)+...
                  F32*F21*(1-F44)+...
                  F34*F41*(1-F22)+...
                  F32*F24*F41+F34*F42*F21);

n21 = (e1.*e2-e2)*(F31*F12*(1-F44)+...
                  F31*F14*F42+F34*F41*F12-F14*F41*F32);

d0  =    (e2-e1.*e2)*F14*F41-(1-F44);

d2  =        (1-e2)*(F24*F42 + F22*(1-F44));

dsq = (1-e1).*(1-e2)*(F12*F21*(1-F44) + F14*F41*(1-F22) + F12*F24*F41 + F14*F42*F21);

na = n0+n12+n21;
da = d0+d2 +dsq;
ea = na./da;

end
end