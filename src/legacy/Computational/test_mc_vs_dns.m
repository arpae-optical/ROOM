%%% This script runs MCTest.m and DNSTest.m, comparing accuracy and     %%%
%%% speed using monte carlo and direct numerical solution of the        %%%
%%% approximate geometry, for a hexagonal prism enclosure.              %%%
%%% Future work: Test with cylinder in middle of enclosure for          %%%
%%%              occlusion performance.                                 %%%

% First, let's pick a single geometry and calculate the exact view factors:
L_H = sqrt(3); 
tic
FWW = FhexWW(L_H); % Exact solution for view factors wall-to-walls of hex
Tex = toc;

%% MC Test
NMC = [10^3,floor(10^3.5),10^4,floor(10^4.5),10^5,floor(10^5.5),10^6]; % Test points for Monte Carlo
WMC = zeros(size(NMC));  % Numeric wall-to-wall view factor results (MC)
TMC = WMC;               % Wall time for MC view factor test
for ii=1:length(NMC)
    tic
    WMC(ii)=MCTest(NMC(ii),L_H);
    TMC(ii)=toc;
end
EMC = abs(WMC-FWW)/FWW;
figure
loglog(NMC,EMC,NMC,TMC)
legend('Error','Wall Time','Location','E')
xlabel('Number of Bundles, N (-)')
ylabel('Relative Error (-), Wall Time (s)')
OEM=log(EMC(2)/EMC(end))/log(NMC(end)/NMC(2));
OTM=log(TMC(end)/TMC(2))/log(NMC(end)/NMC(2));
OETM=OEM/OTM;
disp(['MC Error Order (1/N)^(' num2str(OEM) ')'])
disp(['MC Wall Time (N)^(' num2str(OTM) ')'])
disp(['MC Error order (1/t)^(' num2str(OETM) ')'])

%% DNS Test
NDNS= [4,8,16,32,64];
WDNS= zeros(size(NDNS));
TDNS= WDNS;
for ii=1:length(NDNS)
    tic
    WDNS(ii)=DNSTest(NDNS(ii),L_H);
    TDNS(ii)=toc;
end
EDNS = abs(WDNS-FWW)/FWW;
figure
loglog(NDNS.^2,EDNS,NDNS.^2,TDNS);
legend('Error','Wall Time','Location','E')
xlabel('Number of Cells, N (-)')
ylabel('Relative Error (-), Wall Time (s)')
OED=log(EDNS(2)/EDNS(end))/log(NDNS(end)/NDNS(2));
OTD=log(TDNS(end)/TDNS(2))/log(NDNS(end)/NDNS(2));
OETD=OED/OTD;
disp(['DNS Error Order (1/N)^(' num2str(OED) ')'])
disp(['DNS Wall Time (N)^(' num2str(OTD) ')'])
disp(['DNS Error order (1/t)^(' num2str(OETD) ')'])
