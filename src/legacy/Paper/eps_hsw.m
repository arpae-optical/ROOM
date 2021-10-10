% epsa.m
syms q1 q2 q3 e1 e2 e4 EE EA F12 F13 F14 F21 F22 F23 F24
syms F31 F32 F34 F41 F42 F43 F44

q4=0;
e3=1;
F11=0;
F33=0;

nre1 = q1/e1-F11*(1/e1-1)*q1-F12*(1/e2-1)*q2-...
             F13*(1/e3-1)*q3-F14*(1/e4-1)*q4==...
          EE-F11*EE-F12*EE-F13*0-F14*EA;
nre2 = q2/e2-F21*(1/e1-1)*q1-F22*(1/e2-1)*q2-...
             F23*(1/e3-1)*q3-F24*(1/e4-1)*q4==...
          EE-F21*EE-F22*EE-F23*0-F24*EA;
nre3 = q3/e3-F31*(1/e1-1)*q1-F32*(1/e2-1)*q2-...
             F33*(1/e3-1)*q3-F34*(1/e4-1)*q4==...
           0-F31*EE-F32*EE-F33*0-F34*EA;
nre4 = q4/e4-F41*(1/e1-1)*q1-F42*(1/e2-1)*q2-...
             F43*(1/e3-1)*q3-F44*(1/e4-1)*q4==...
          EA-F41*EE-F42*EE-F43*0-F44*EA;
[~,~,q3s,~] = solve([nre1 nre2 nre3 nre4],[q1 q2 q3 EA]);

e_a = simplify(-q3s/EE);

% pretty(e_a)

%% Section 2: actual equation
n0  = -e1*e2*(F31*(1-F44)+F34*F41)-e2*(F32*(1-F44)+F34*F42);

n12 = (e1*e2-e1)*(F31*(1-F22-F44-F24*F42+F22*F44)+...
                  F32*F21*(1-F44)+...
                  F34*F41*(1-F22)+...
                  F32*F24*F41+F34*F42*F21);

n21 = (e1*e2-e2)*(F31*F12*(1-F44)+...
                  F31*F14*F42+F34*F41*F12-F14*F41*F32);

d0  =    (e2-e1*e2)*F14*F41-(1-F44);

d2  =        (1-e2)*(F24*F42 + F22*(1-F44));

dsq = (1-e1)*(1-e2)*(F12*F21*(1-F44) + F14*F41*(1-F22) + F12*F24*F41 + F14*F42*F21);

na = n0+n12+n21;
da = d0+d2 +dsq;
ea = na/da;

pretty(simplify((ea-e_a)*da))
%% All sections

n1_12=(e1*e2-e1)*(F31*(F22*F44-F24*F42)-F32*F21*F44-...
                  F34*F41*F22+F32*F24*F41+F34*F42*F21);
n2_12=(e1*e2-e2)*(-F32*F14*F41-F31*F12*F44-F31*F14*(1-F42)+F34*F41*F12);
n12  =(1-e1)*(1-e2)*(F31*F14*(1-F22)-F34*F12*F21+F31*F12*F24+F32*F21*F14);
nr   = F34*(1-F22)+F32*F24+e1*(F31*F44-F34*F41)+...
       e2*(F32*F44-F34*F42+F34*F22-F32*F24);

d2   =(1-e2)*(F24*F42-F22*F44);
d12  =(1-e1)*(1-e2)*(F14*F41*(1-F22)-F12*F21*F44+F12*F24*F41+F14*F42*F21);
dr = (e1*e2-e2)*(-F14*F41);
ea=(n1_12+n2_12+n12+nr)/(d2+d12+F44+dr);
pretty(simplify(ea-e_a))


%% Check for 0

pretty(simplify(subs(ea,[e1,e2],[0,0])))

%% Wall only
% epsa.m
syms q1 q3 e1 e4 EE EA F13 F14 
syms F31 F34 F41 F42 F43 F44
F43=F41;
F34=F14;
F31=F13;
q4=0;
e3=1;
F11=0;
F33=0;

nre1 = q1/e1==EE-F14*EA;
nre3 = q3-F13*(1/e1-1)*q1==-F13*EE-F14*EA;
nre4 = -F41*(1/e1-1)*q1==EE-F41*EE-F44*EA; %-F42*EE
[~,q3s,~] = solve([nre1 nre3 nre4],[q1 q3 EA]);

e_a = simplify(-q3s/EE);

pretty(e_a)

%% ball only
%epsam?

syms q1 q2 q3 EE F12 F13 F21 F22 F23 F31 F32


n1 = -e1*(F31*(1-F22)+F32*F21);
n2 = -e2*(F32+F31*F12);
n12= e1*e2*(F32*F21+F31*F12-F31*F22);

d12 = (1-e1)*(1-e2)*F12*F21;
d20 = (1-e2)*F22-1;

ea=(n1+n2+n12)/(d12+d20);
