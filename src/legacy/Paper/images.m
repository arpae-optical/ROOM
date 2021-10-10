%% Matlab setup:
set(groot,'defaultLineLineWidth',2.0)
set(groot,'defaultAxesFontName','Times')
set(groot,'defaultAxesFontSize',12)

%% check it out
cd ../Paper/

pullal
pullalumina
pullcsic
pullsi
pullti
pullinc
pullni
pullcu

plot(l_al,e_al,l_al2o3,e_al2o3,l_csic,e_csic,l_si,e_si,l_ti,e_ti,l_inc,e_inc,l_ni,e_ni,l_cu,e_cu)
legend('al','alu','csic','si','ti','inc','ni','cu')

%% solar cooling: low emissivity near solar max (0.5um) high emissivity 
% near blackbody max at ~300-400K (7-10um)

cd ../Paper/
pullal
figure('Position',[60 60 600 600])
pullalumina

l_ipsc=[0,6,6,10,11];
e_ipsc=[.1,.1,.9,.9,.1];

e_alu=interp1(l_al2o3,e_al2o3,l_al,'makima',min(e_al2o3));
cd ../Computational/
e_opt24=epshsw(e_al,e_alu,0.25,4);
e_opt54=epshsw(e_al,e_alu,0.5,4);
e_opt74=epshsw(e_al,e_alu,0.75,4);
cd ../Paper/
plot(l_al,e_al,'--k',l_al,e_opt24,'-b',l_al,e_opt54,'-c',l_al,e_opt74,'-g');

legend('Aluminum Emissivity','Walls + Alumina Spheres: D/L=0.25','Walls + Alumina Spheres: D/L=0.50','Walls + Alumina Spheres: D/L=0.75')
xlabel('Wavelength, \lambda (\mum)')
ylabel('Spectral Emissivity, \epsilon_\lambda (-)')
axis([0 12 0 1])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','FontName'),'FontName','Times')
%% If happy, print this picture to file
print -dtiff SolarCooling
%% analyze
i_ir=find(abs(l_al-7)==min(abs(l_al-7)));
i_sm=find(abs(l_al-0.5)==min(abs(l_al-0.5)));

e_sm_base=e_al(i_sm);
e_ir_base=e_al(i_ir);

e_sm_24=e_opt24(i_sm)/e_sm_base;
e_ir_24=e_opt24(i_ir)/e_ir_base;

e_sm_54=e_opt54(i_sm)/e_sm_base;
e_ir_54=e_opt54(i_ir)/e_ir_base;

e_sm_74=e_opt74(i_sm)/e_sm_base;
e_ir_74=e_opt74(i_ir)/e_ir_base;

%% solar absorber ideal spectra: high in solar max (~0.5um), low in IR (~8-10um)

pullsi
% pullinc
pullni
% pullcu



e_ni=interp1(l_ni,e_ni,l_si,'makima');
cd ./../../
e_opt25=epshsw(e_si,e_ni,.25,2);
e_opt55=epshsw(e_si,e_ni,.5,2);
e_opt15=epshsw(e_si,e_ni,1,2);

cd ./Legacy/Paper/
figure('Position',[60 60 600 600])
plot(l_si,e_si,'--k',l_si,e_opt25,'-b',l_si,e_opt55,'-c',l_si,e_opt15,'-g')
legend('Silicon Emissivity','Walls + Ni Spheres: D/L=0.25',...
       'Walls + Ni Spheres: D/L=0.5','Walls + Ni Spheres: D/L=1.0',...
       'Location','NE')
xlabel('Wavelength, \lambda (\mum)')
ylabel('Spectral Emissivity, \epsilon_\lambda (-)')
axis([0 12 0 1])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','FontName'),'FontName','Times');

%% If happy, print this picture to file
print -dtiff SolarCollect

%% hypersonic skin ideal spectra (high for heat rejection around 1.5-3um, low for parasitic heating (7-10)
pullti
pullalumina

e_ti=interp1(l_ti,e_ti,l_al2o3,'makima',min(e_ti));
cd ./../../
e_opt25=epshsw(e_al2o3,e_al2o3,.5,0);
e_opt55=epshsw(e_al2o3,e_al2o3,.75,0);
e_opt15=epshsw(e_al2o3,e_al2o3,1,0);
cd ./Legacy/Paper/
figure('Position',[60 60 600 600])
plot(l_al2o3,e_al2o3,'--k',l_al2o3,e_opt25,'-b',l_al2o3,e_opt55,'-c',l_al2o3,e_opt15,'-g');

legend('Alumina Emissivity','Alumina Spheres: D/L=0.5',...
       'Alumina Spheres: D/L=0.75','Alumina Spheres: D/L=1.0',...
       'Location','SE')
xlabel('Wavelength, \lambda (\mum)')
ylabel('Spectral Emissivity, \epsilon_\lambda (-)')
axis([0 12 0 1])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','FontName'),'FontName','Times')

%% If happy, print this picture to file
print -dtiff Hypersonic

%% tpv emitter ideal spectra, high around 1000k (3um) low around 500K (6um)
pullni
pullti

e_al=interp1(l_ni,e_ni,l_ti,'makima');
cd ./../../
e_opt25=epshsw(e_ti,e_al,.25,2);
e_opt55=epshsw(e_ti,e_al,.5,2);
e_opt15=epshsw(e_ti,e_al,1,2);
cd ./Legacy/Paper/

figure('Position',[60 60 600 600])
plot(l_ti,e_ti,'--k',l_ti,e_opt25,'-b',l_ti,e_opt55,'-c',l_ti,e_opt15,'-g');

legend('Titanium Emissivity','Walls + Ni Spheres: D/L=0.25',...
       'Walls + Ni Spheres: D/L=0.5','Walls + Ni Spheres: D/L=1.0',...
       'Location','N')
xlabel('Wavelength, \lambda (\mum)')
ylabel('Spectral Emissivity, \epsilon_\lambda (-)')
axis([0 12 0 1])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','FontName'),'FontName','Times')

%% If happy, print this picture to file
print -dtiff TPV