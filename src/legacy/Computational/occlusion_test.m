%%% This script runs MCHexCyl.m and DNSHexCyl.m, comparing view factors %%%
%%% and speed using monte carlo and direct numerical solution of the    %%%
%%% approximate geometry, for a hexagonal prism enclosure with a        %%%
%%% cylindrical occlusion at the origin.                                %%%

% First, let's pick a single geometry to compare the view factors:
L_H = sqrt(3); 

%% MC Test
NMC = [10^3,ceil(10^3.5),10^4,ceil(10^4.5),10^5,ceil(10^5.5),10^6]; % Test points for Monte Carlo
WMC = zeros(size(NMC));  % Numeric wall-to-wall view factor results (MC)
TMC = WMC;               % Wall time for MC view factor test
for ii=1:length(NMC)
    tic
    WMC(ii)=MCHexCyl(NMC(ii),L_H,1);
    TMC(ii)=toc;
end
%% DNS Test
NDNS= [2^2,ceil(2^2.5),2^3,ceil(2^3.5),2^4,ceil(2^4.5),2^5];
WDNS= zeros(size(NDNS));
TDNS= WDNS;
for ii=1:length(NDNS)
    tic
    WDNS(ii)=DNSHexCyl(NDNS(ii),L_H,1);
    TDNS(ii)=toc;
end
%% Plot
figure
loglog(NMC,WMC,NMC,TMC,10.*NDNS.^2,WDNS,10.*NDNS.^2,TDNS)
legend('MC VF','MC Time','DNI VF','DNI Time','Location','SE')
xlabel('Number of Bundles/Cells, N (-)')
ylabel('View Factor (-), Wall Time (s)')
title({'Comparison of Direct Numeric Integration and Monte Carlo';...
    'for the Calculation of Occluded View Factors,';...
    'Cylinder in Hexagonal Enclosure, L/H=1.7, D/H=1'})
axis([10^2 10^6 10^-4 10^1])