% Onur Kucuktas

% We have relaxed the black substrate assumption. May God help me with
% geometry.

clc; clear all; format short;
tic
    
% Given/chosen variables
N = 10000; % Number of bundles
n1 = 1; % Index of refraction of 1st medium
n2 = 1.5; % Index of refraction of 2nd medium, chose Pyrex arbitrarily
n3 = 1; % Index of refraction of 3rd medium
r2 = .00015; % outer radius in millimeters
r1 = r2*.92; % inner radius in millimeters
lambda0 = 550E-9;% wavelength of incident light in vacuum in meters
k2 = 2.918E-5; % absorptive index of medium 2 for use in transmisivity calculation
            % much smaller than one in order to be weakly absorbing
            % material
k = (4*pi*k2)/lambda0; % for use in transmisivity
%k = .1;
thetaCrit = asind(n3/n2); % any angle higher than this results in TIR

% Initialize incident angle vector
theta1 = [];

% Loop N times and fill the theta1 vector
for i = 1:N
  theta1(i) = asind((i-1)/N);
end

% Because of the perfectly specular surface, angle of reflection is also
% equal to theta 1

% Define theta 2, this is the angle of the transmitted beam into medium 2
theta2 = asind((n1/n2*sind(theta1)));

% Define theta3, this is the angle of incidence of the transmitted beam
% within medium 2 onto medium 3
theta3 = asind((r2/r1)*sind(theta2));

% NOTE: theta3 has complex solutions towards the end. 

% Reflection angle is the same because of Snell's law and perfectly
% specular surface

% Define theta4, this is the angle of the transmitted beam into medium 3
theta4 = asind((n2/n3)*sind(real(theta3)));

% Define theta5, angle of the transmitted beam out of the shell into medium
% 3
theta5 = asind((n2/n1)*sind(real(theta2)));

% Identify how far each bundle travels within the shell.
for i = 1:N
    distance(i) = r2*cosd(theta2(i)) - r1*cosd(theta3(i));
end

TIRBundleIndexes = [];
j = 1;

for i = 1:N

    if theta3(i) > thetaCrit
        TIRBundleIndexes(j) = i;
        j = j + 1 ;
    end
end

% Now for every theta1, theta6(i) needs to be calculated. Each theta 1 will
% produce many theta6 values. For now, we'll truncate at 5 theta6 values.
% Analyzing the Transmittance series will tell us how many terms we can
% truncate at.

truncation = 7;
theta6=[];%zeros(truncation, N-length(TIRBundleIndexes));
theta7=[];%zeros(truncation, N-length(TIRBundleIndexes));
beta=[];
alpha = [];
returnPos = {};
% Loop thru all theta1's which don't correspond to bundles undergoing
% TIR
for i = 1:TIRBundleIndexes(1)-1
      j = 1;
      alpha = 90-theta1(i);
      while j <= truncation
          % short position vector finding, yields a theta6           
          Pold = [r2*cosd(alpha), r2*sind(alpha)];
          beam = [distance(i)*cosd(180+alpha+theta2(i)), distance(i) * sind(180+alpha+theta2(i))];
          Pnew = Pold + beam;
          alpha = atan2d(Pnew(2),Pnew(1));
          outsideAngle = theta4(i) + alpha;%180-theta4(i)-alpha;
          
          % beam moving left
          if outsideAngle < 90 
              theta6(j,i) = 90-outsideAngle;
              beta(j,i) = 90 - theta6(j,i);
              l = Pnew(1) - Pnew(2)*tand(theta6(j,i));
              theta7(j,i) = asind((l*sind(beta(j,i)))/r1);
              returnPos{i,j} = [r1*cosd(180 - beta(j,i) - theta7(j,i)),r1*sind(180 - beta(j,i) - theta7(j,i))];%[r1*cosd(180 - (beta+2*theta6(j,i)) - theta7(j,i)),r1*sind(180 - (beta+2*theta6(j,i)) - theta7(j,i))];
          
          % beam moving right
          else                   
              theta6(j,i) = outsideAngle-90;
              beta(j,i) = 90 + theta6(j,i);
              l = Pnew(1) + Pnew(2)*tand(theta6(j,i));
              theta7(j,i) = asind((l*sind(90+theta6(j,i))/r1));
              returnPos{i,j} = [r1*cosd(180 - beta(j,i) - theta7(j,i)),r1*sind(180 - beta(j,i) - theta7(j,i))];
          end
          
          %long position vector finding
          Pold = Pnew;
          beam = [distance(i) * cosd(alpha-theta3(i)), distance(i) * sind(alpha-theta3(i))];
          Pnew = Pold + beam;
          alpha = atan2d(Pnew(2),Pnew(1))
          %Pold = Pnew;
          
          %increase iterator variable
          j = j + 1;
          
      end
        
end
    

% %short
% alpha = 90 - theta1(50);
% Pold = [r2*cosd(alpha), r2*sind(alpha)];
% beam = [distance(50)*cosd(180+alpha+theta2(50)), distance(50) * sind(180+alpha+theta2(50))];
% Pnew = Pold + beam;
% 
% %long
% Pold = Pnew;
% alpha = atan2d(Pold(2),Pold(1));
% beam = [distance(50) * cosd(alpha-theta3(50)), distance(50) * sind(alpha-theta3(50))];
% Pnew = Pold + beam;
% 
% %short
% Pold = Pnew;
% alpha = atan2d(Pold(2),Pold(1));
% beam = [distance(50)*cosd(180+alpha+theta2(50)), distance(50) * sind(180+alpha+theta2(50))];
% Pnew = Pold + beam;
% 
% %long
% Pold = Pnew;
% alpha = atan2d(Pold(2),Pold(1));
% beam = [distance(50) * cosd(alpha-theta3(50)), distance(50) * sind(alpha-theta3(50))];
% Pnew = Pold + beam;





% Identify by index which of the incident bundles on the inner interface
% are going to undergo total internal reflection (theta3 > thetaCrit). When
% TIR happens, theta4 becomes meaningless and will be complex.

TIRBundleIndexes = [];
j = 1;

for i = 1:N

    if theta3(i) > thetaCrit
        TIRBundleIndexes(j) = i;
        j = j + 1 ;
    end
end

% Identify how far each bundle travels within the shell.
for i = 1:N
    distance(i) = r2*cosd(theta2(i)) - r1*cosd(theta3(i));
end

% Compute the reflectivities and transmisivity for each bundle
for i = 1:N
    rho12(i) = 0.5 * ( ((tand(theta1(i) - theta2(i)))^2)/((tand(theta1(i) + theta2(i)))^2)...
    + ((sind(theta1(i) - theta2(i)))^2)/((sind(theta1(i) + theta2(i)))^2));

    rho21(i) = 0.5 * ( ((tand(theta2(i) - theta5(i)))^2)/((tand(theta2(i) + theta5(i)))^2)...
    + ((sind(theta2(i) - theta5(i)))^2)/((sind(theta2(i) + theta5(i)))^2));

    if theta3(i) < thetaCrit

        rho23(i) = 0.5 * ( ((tand(theta3(i) - theta4(i)))^2)/((tand(theta3(i) + theta4(i)))^2)...
        + ((sind(theta3(i) - theta4(i)))^2)/((sind(theta3(i) + theta4(i)))^2));

    else
        rho23(i) = 1;
    end

    tau(i) = exp(-k*r2*(cosd(theta2(i)) - (r1/r2)*cosd(theta3(i))));
end

% Will compute total apparent reflectance, transmittance, and absorptance

% Initialize

R = [];T = []; A=[];

% Compute R, T, A for each bundle
for i = 1:N

    R(i) = (rho12(i) + (rho23(i) * (tau(i))^2)*(1 - rho21(i) - rho12(i)))...
        /(1 - rho21(i)*rho23(i)*(tau(i))^2);

    T(i) = ((1 - rho12(i))*(1-rho23(i))*tau(i))/ (1 - rho21(i)*rho23(i)*(tau(i))^2);

    A(i) = 1 - R(i) - T(i);   

    %A(i) = ((1-tau(i))*(1-rho12(i))*(1+rho23(i)*tau(i)))/ (1 - rho21(i)*rho23(i)*(tau(i))^2);
end
% truncation = [];
% for i = 1:N
%    a = (1-rho12(i))*(1-rho23(i))*tau(i);
%    truncation(i) = (log(.000001/a)/log(rho12(i)*rho23(i)*tau(i)*tau(i)))+1;
% end
% Now for every theta1, theta6(i) needs to be calculated. Each theta 1 will
% produce many theta6 values. For now, we'll truncate at 5 theta6 values.
% Analyzing the Transmittance series will tell us how many terms we can
% truncate at.

% truncation = 5;
% theta6=zeros(truncation, N-length(TIRBundleIndexes));
% theta7=zeros(truncation, N-length(TIRBundleIndexes));
% returnPos = {};
% % Loop thru all theta1's which don't correspond to bundles undergoing
% % TIR
% for i = 1:TIRBundleIndexes(1)
%       j = 1;
%       alpha = 90-theta1(i);
%       while j <= truncation
%           % short position vector finding, yields a theta6
%           Pold = [r2*cosd(alpha), r2*sind(alpha)];
%           beam = [distance(i)*cosd(180+alpha+theta2(i)), distance(i) * sind(180+alpha+theta2(i))];
%           Pnew = Pold + beam;
%           alpha = atan2d(Pnew(2),Pnew(1));
%           outsideAngle = 180-theta4(i)-alpha;
%           
%           % beam moving left
%           if outsideAngle < 90 
%               theta6(j,i) = 90-outsideAngle;
%               beta = 90 - theta6(j,i);
%               l = Pnew(1) - Pnew(2)*tand(theta6(j,i));
%               theta7(j,i) = asind((l*sind(beta))/r1);
%               returnPos{i,j} = [r1*cosd(180 - (beta+2*theta6(j,i)) - theta7(j,i)),r1*sind(180 - (beta+2*theta6(j,i)) - theta7(j,i))];
%           
%           % beam moving right
%           else                   
%               theta6(j,i) = outsideAngle-90;
%               beta = 90 - theta6(j,i);
%               l = Pnew(1) + Pnew(2)*tand(theta6(j,i));
%               theta7(j,i) = asind((l*sind(90+theta6(j,i))/r1));
%               returnPos{i,j} = [r1*cosd(180 - beta - theta7(j,i)),r1*sind(180 - beta - theta7(j,i))];
%           end
%           
%           %long position vector finding
%           Pold = Pnew;
%           beam = [distance(i) * cosd(alpha-theta3(i)), distance(i) * sind(alpha-theta3(i))];
%           Pnew = Pold + beam;
%           alpha = atan2d(Pnew(2),Pnew(1));
%           Pold = Pnew;
%           
%           %increase iterator variable
%           j = j + 1;
%           
%       end
%         
% end

% % Compute the total apparent properties
% Rtot = sum(R(2:N)) / (N-1);
% 
% Ttot = sum(T(2:N))/ (N-1);
% 
% Atot = sum(A(2:N))/ (N-1);
% 
% % Results Matrix
% Results = [rho12; rho21; rho23; tau; R; T; A];
% 
% % Number of terms to convergence for transmittance
% 
% ternNumber50 = log(0.0000001)/log(rho12(50)*rho23(50)*tau(50)^2)+1
% 
% ternNumber50 = log(0.0000001)/log(rho12(60)*rho23(60)*tau(60)^2)+1
% 
% ternNumber50 = log(0.0000001)/log(rho12(90)*rho23(90)*tau(90)^2)+1

%Print out the total apparent properties

% fprintf('Total Apparent Reflectance (R): %f\n', Rtot);
% fprintf('Total Apparent Transmittance (T): %f\n', Ttot);
% fprintf('Total Apparent Absorbtance (A): %f\n\n', Atot);

% % Plot Total Absorptance vs Normalized Shell Thickness
% plot(normalizedShellThickness,Adata); xlabel ('r1/r2','fontweight','bold'); ylabel ('Total Absorptance'...
%     ,'fontweight', 'bold');
% title('Total Apparent Absorptance vs. Normalized Shell Thickness');

toc







