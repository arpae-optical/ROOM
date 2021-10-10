% pullcsic.m
% loads spectral emissivity values into workspace
%
% data taken from:
%"EFFECTS OF TEMPERATURE ON THE SPECTRAL EMISSIVITY OF C/SiC COMPOSITES"
%YUFENG ZHANG
%Ceramics-Silik�ty
%journal article
%15 September 2015
%https://www.irsm.cas.cz/materialy/cs_content/2016_doi/Zhang_CS_2016_0023.pdf

%% c/sic 1000k
k1000=[3.235555555555557, 0.5527421236872812;...
3.5288888888888894, 0.6024504084014002;...
3.662222222222223, 0.6416569428238039;...
3.813333333333333, 0.6794632438739789;...
3.955555555555556, 0.6948658109684948;...
4.106666666666667, 0.7116686114352393;...
4.462222222222223, 0.7571761960326722;...
4.675555555555556, 0.7662777129521587;...
4.897777777777779, 0.7837806301050175;...
5.13777777777778, 0.798483080513419;...
5.395555555555556, 0.8075845974329057;...
5.697777777777779, 0.8313885647607935;...
6.026666666666668, 0.8383897316219371;...
6.40888888888889, 0.8635939323220537;...
6.826666666666668, 0.8502917152858811;...
7.315555555555557, 0.8460910151691949;...
7.866666666666667, 0.8551925320886815;...
8.515555555555558, 0.8684947491248542;...
9.27111111111111, 0.8992998833138858;...
10.195555555555556, 0.9476079346557762;...
11.297777777777776, 0.7151691948658111;...
12.666666666666666, 0.675262543757293;...
14.435555555555554, 0.806184364060677;...
16.773333333333333, 0.827887981330222;...
20.04444444444444, 0.887397899649942];
l_csic=k1000(:,1);
e_csic=k1000(:,2);
clear k1000
% 
% %% csic 1200k
% k1200=[3.1377777777777793, 0.512835472578763;...
% 3.315555555555559, 0.5912485414235706;...
% 3.5288888888888894, 0.5884480746791131;...
% 3.662222222222223, 0.61365227537923;...
% 3.813333333333333, 0.6409568261376896;...
% 3.955555555555556, 0.6710618436406067;...
% 4.106666666666667, 0.6857642940490083;...
% 4.462222222222223, 0.7214702450408401;...
% 4.675555555555556, 0.7396732788798135;...
% 4.897777777777779, 0.7571761960326722;...
% 5.13777777777778, 0.7634772462077013;...
% 5.386666666666667, 0.7781796966161028;...
% 5.697777777777779, 0.7907817969661611;...
% 6.026666666666668, 0.8033838973162195;...
% 6.40888888888889, 0.8131855309218204;...
% 6.826666666666668, 0.8152858809801635;...
% 7.315555555555557, 0.8187864644107352;...
% 7.866666666666667, 0.8285880980163363;...
% 8.506666666666668, 0.8411901983663945;...
% 9.262222222222222, 0.8684947491248542;...
% 10.186666666666667, 0.9280046674445742;...
% 11.297777777777776, 0.6815635939323222;...
% 12.666666666666666, 0.661260210035006;...
% 14.435555555555554, 0.7970828471411902;...
% 16.773333333333333, 0.7774795799299887;...
% 20.000000000000004, 0.9553092182030342];
% lam1200kcsic=k1200(:,1);
% emm1200kcsic=k1200(:,2);
% clear k1200
% 
% %% 1400k
% k1400=[3.1377777777777793, 0.4358226371061844;...
% 3.315555555555559, 0.5142357059509919;...
% 3.5288888888888894, 0.5765460910151692;...
% 3.662222222222223, 0.5968494749124854;...
% 3.813333333333333, 0.624154025670945;...
% 3.955555555555556, 0.6577596266044341;...
% 4.106666666666667, 0.6717619603267213;...
% 4.462222222222223, 0.7046674445740957;...
% 4.675555555555556, 0.724970828471412;...
% 4.897777777777779, 0.7396732788798133;...
% 5.13777777777778, 0.752975495915986;...
% 5.395555555555556, 0.7648774795799302;...
% 5.697777777777779, 0.7893815635939324;...
% 6.026666666666668, 0.7963827304550759;...
% 6.40888888888889, 0.8236872812135356;...
% 6.826666666666668, 0.8096849474912488;...
% 7.315555555555557, 0.8089848308051342;...
% 7.866666666666667, 0.8173862310385066;...
% 8.506666666666668, 0.8257876312718788;...
% 9.262222222222222, 0.8523920653442242;...
% 10.186666666666667, 0.905600933488915;...
% 11.297777777777776, 0.7004667444574098;...
% 12.666666666666666, 0.655659276546091;...
% 14.435555555555554, 0.7900816802800465;...
% 16.773333333333333, 0.8173862310385068;...
% 20.000000000000004, 0.9308051341890318];
% lam1400kcsic=k1400(:,1);
% emm1400kcsic=k1400(:,2);
% clear k1400
% 
% %% 1600k
% k1600=[3.1377777777777793, 0.4582263710618437;...
% 3.315555555555559, 0.5233372228704785;...
% 3.4311111111111114, 0.5653442240373396;...
% 3.5288888888888894, 0.5968494749124853;...
% 3.662222222222223, 0.6248541423570595;...
% 3.813333333333333, 0.6563593932322052;...
% 3.955555555555556, 0.6738623103850641;...
% 4.106666666666667, 0.6941656942823805;...
% 4.284444444444445, 0.6864644107351227;...
% 4.462222222222223, 0.7256709451575263;...
% 4.675555555555556, 0.7396732788798135;...
% 4.897777777777779, 0.7543757292882147;...
% 5.13777777777778, 0.7697782963827303;...
% 5.395555555555556, 0.7788798133022172;...
% 5.697777777777779, 0.8019836639439907;...
% 6.026666666666668, 0.8096849474912486;...
% 6.40888888888889, 0.8341890315052509;...
% 6.826666666666668, 0.8166861143523922;...
% 7.315555555555557, 0.8138856476079346;...
% 7.866666666666667, 0.8215869311551927;...
% 8.506666666666668, 0.8334889148191367;...
% 9.262222222222222, 0.857992998833139;...
% 10.186666666666667, 0.9021003500583433;...
% 11.297777777777776, 0.7221703617269546;...
% 12.666666666666666, 0.6563593932322054;...
% 14.435555555555554, 0.7963827304550757;...
% 16.773333333333333, 0.8299883313885652;...
% 20.000000000000004, 0.8866977829638276];
% lam1600kcsic=k1600(:,1);
% emm1600kcsic=k1600(:,2);
% clear k1600