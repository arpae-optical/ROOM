%%% This script runs MCTest.m and DNSTest.m, comparing accuracy and     %%%
%%% speed using monte carlo and direct numerical solution of the        %%%
%%% approximate geometry, for a hexagonal prism enclosure.              %%%
%%% Future work: Test with cylinder in middle of enclosure for          %%%
%%%              occlusion performance.                                 %%%

% First, let's pick a single geometry and calculate the exact view factors:
L_H = sqrt(3); 

[FWW,tex] = FhexWW(L_H); % Exact solution for view factors wall-to-walls of hex

% FWB = (1-FWW)/2;   % Exact for wall-to-base of hex
% FWT = FWB;         % Exact for wall-to-top

NMC = [10^3,floor(10^3.5),10^4,floor(10^4.5),10^5,floor(10^5.5),10^6]; % Test points for Monte Carlo
WMC = zeros(size(NMC));  % Numeric wall-to-wall view factor results (MC)
% BMC = WMC;               % Numeric wall-to-base view factor results (MC)
TMC = WMC;               % Wall time for MC view factor test
for ii=1:length(NMC)
    [WMC(ii),~,~,TMC(ii)]=MCTest(NMC(ii),L_H);
end
EMC = abs(WMC-FWW)/FWW;
figure
loglog(NMC,EMC,NMC,TMC)
legend('Error','Wall Time','Location','E')
xlabel('Number of Bundles, N (-)')
ylabel('Relative Error (-), Wall Time (s)')
OEM=log(EMC(1)/EMC(end))/log(NMC(end)/NMC(1));
OTM=log(TMC(end)/TMC(1))/log(NMC(end)/NMC(1));
OETM=OEM/OTM;
disp(['MC Error Order (1/N)^(' num2str(OEM) ')'])
disp(['MC Wall Time (N)^(' num2str(OTM) ')'])
disp(['MC Error order (1/t)^(' num2str(OETM) ')'])

%% DNS Test
NDNS= [10,14,20,28,40];
WDNS= zeros(size(NDNS));
TDNS= WDNS;
for ii=1:length(NDNS)
    [WDNS(ii),~,~,TDNS(ii)]=DNSTest(NDNS(ii),L_H);
end
EDNS = abs(WDNS-FWW)/FWW;
figure
loglog(NDNS.^2,EDNS,NDNS.^2,TDNS);
legend('Error','Wall Time','Location','E')
xlabel('Number of Cells, N (-)')
ylabel('Relative Error (-), Wall Time (s)')
OED=log(EDNS(1)/EDNS(end))/log(NDNS(end)/NDNS(1));
OTD=log(TDNS(end)/TDNS(1))/log(NDNS(end)/NDNS(1));
OETD=OED/OTD;
disp(['DNS Error Order (1/N)^(' num2str(OED) ')'])
disp(['DNS Wall Time (N)^(' num2str(OTD) ')'])
disp(['DNS Error order (1/t)^(' num2str(OETD) ')'])
