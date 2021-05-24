
# Onur Kucuktas

# Hemispherical Shell on Non-Black Substrate capable of reflection

clc; clear all; format short
tic()

#NST = range(.93,.99,length=10)
# ref = range(1.1,1.8,length=50)

# for x = 1:length(ref)


# Given/chosen variables
N = 10000; # Number of bundles
n1 = 1; # Index of refraction of 1st medium
n2 = 1.1;#1.463;#1.5; # Index of refraction of 2nd medium
n3 = 1; # Index of refraction of 3rd medium
r2 = .1;     # outer radius in meters
r1 = r2*.93; # inner radius in meters
lambda0 = 9.5E-7;#3.1E-6;#4E-6;#550E-9;# wavelength of incident light in vacuum in meters
k2 =0;# 4.17E-6;#1.06E-4;#2.918E-5; # absorptive index of medium 2 for use in transmisivity calculation
            # much smaller than one in order to be weakly absorbing
            # material()
k = (4*pi*k2)/lambda0; # for use in transmisivity
#k = .1
thetaCrit = asind[n3/n2]; # any angle higher than this results in TIR

# **********************************************************************
#              INITIAL SHELL INTERACTIONS
# **********************************************************************

# Initialize incident angle vector
theta1 = []

# Loop N times & fill the theta1 vector
for i = 1:N
  theta1[i] = asind[(i-1)/N]
end

# Because of the perfectly specular surface; angle of reflection is also
# equal to theta 1

# Define theta 2; this is the angle of the transmitted beam into medium 2
theta2 = asind[(n1/n2*sind[theta1])]

# Define theta3; this is the angle of incidence of the transmitted beam
# within medium 2 onto medium 3
theta3 = asind[(r2/r1)*sind[theta2]]

# NOTE: theta3 has complex solutions towards the end. 

# Reflection angle is the same because of Snell's law & perfectly
# specular surface()

# Define theta4; this is the angle of the transmitted beam into medium 3
theta4 = asind[(n2/n3)*sind[real(theta3)]]

# Define theta5; angle of the transmitted beam out of the shell into medium
# 3
theta5 = asind[(n2/n1)*sind[real(theta2)]]

# Identify how far each bundle travels within the shell.
for i = 1:N
    distance[i] = r2*cosd[theta2[i]] - r1*cosd[theta3[i]]
end

# Identify by index which of the incident bundles on the inner interface
# are going to undergo total internal reflection [theta3 .> thetaCrit]. When
# TIR happens; theta4 becomes meaningless & will be complex.

TIRBundleIndexes = []
j = 1

for i = 1:N

    if theta3[i] .> thetaCrit
        TIRBundleIndexes[j] = i
        j = j + 1 
    end
end

# Identify how far each bundle travels within the shell.
for i = 1:N
    distance[i] = r2*cosd[theta2[i]] - r1*cosd[theta3[i]]
end

# Compute the reflectivities & transmisivity for each bundle
for i = 1:N
    rho12[i] = 0.5 * ( ((tand[theta1[i] - theta2[i]])^2)/((tand[theta1[i] + theta2[i]])^2)...
    + ((sind[theta1[i] - theta2[i]])^2)/((sind[theta1[i] + theta2[i]])^2))

    rho21[i] = 0.5 * ( ((tand[theta2[i] - theta5[i]])^2)/((tand[theta2[i] + theta5[i]])^2)...
    + ((sind[theta2[i] - theta5[i]])^2)/((sind[theta2[i] + theta5[i]])^2))

    if theta3[i] .< thetaCrit

        rho23[i] = 0.5 * ( ((tand[theta3[i] - theta4[i]])^2)/((tand[theta3[i] + theta4[i]])^2)...
        + ((sind[theta3[i] - theta4[i]])^2)/((sind[theta3[i] + theta4[i]])^2))

    elseif theta3[i] .> thetaCrit && isreal(theta3[i])
        rho23[i] = 1
    else()
        rho23[i] = 0
    end
    tau[i] = exp(-k * distance[i])
    #tau[i] = exp(-k*r2*(cosd[theta2[i]] - (r1/r2)*cosd[theta3[i]]))
end

# Will compute total apparent reflectance; transmittance; & absorptance

# Initialize

R = [] T = []' A=[]
for i = 1:N
    DENOM[i] = ((i/N)*r2)
end

sumDENOM = sum(DENOM)

# Compute R; T; A for each bundle
for i = 1:N

    R[i] = ((rho12[i] + (rho23[i] * (tau[i])^2)*(1 - rho21[i] - rho12[i]))...
        /(1 - rho21[i]*rho23[i]*(tau[i])^2))* ((i/N)*r2) / sumDENOM

    T[i] = (((1 - rho12[i])*(1-rho23[i])*tau[i])/ (1 - rho21[i]*rho23[i]*(tau[i])^2))*((i/N)*r2) / sumDENOM

    A[i] = ((((i/N)*r2) / sumDENOM )- R[i] - T[i]); 
    
#     DENOM[i] = ((i/N)*r2)

    #A[i] = ((1-tau[i])*(1-rho12[i])*(1+rho23[i]*tau[i]))/ (1 - rho21[i]*rho23[i]*(tau[i])^2)
end


Rtot = sum(R[2:N]) #/ #(DENOM*(N-1))

Ttot = sum(T[2:N])#/ #(N-1)#(DENOM*(N-1))

Atot = sum(A[2:N])#/#(N-1)#(DENOM *(N-1))

j = 1
for i = 1:N
    
    alpha = 90-theta1[i]
    Pold= [r2*cosd[alpha], r2*sind[alpha]]
    beam = [distance[i]*cosd[180+alpha+theta2[i]], distance[i] * sind[180+alpha+theta2[i]]]
    Pnew = Pold + beam
    
    if Pnew[2] .< 0
       stop = "yes"
       break()
    end
    j = j+1
end

#Rtot+Ttot
# cutoff[x] = j-1

# end
# 
# plot(ref, ((N - cutoff)/N)*100)
# title("Percentage of Direct Impact Beams vs. n2")
# xlabel("n2")
# ylabel("Direct Impact Beams Percentage [#]")


# **********************************************************************
#              SUBSTRATE INTERACTIONS
# **********************************************************************

# Complex index of refraction for substrate
n4 = 3.4302;#3.4302;#4.087
k4 = 0;#0;#5.395

truncation = []

if isempty(TIRBundleIndexes)
    y = N
else()
    y = (TIRBundleIndexes[1]-1)
end
for i = 1:y
   a = (1-rho12[i])*(1-rho23[i])*tau[i]
   truncation[i] = (log(.000001/a)/log(rho12[i]*rho23[i]*tau[i]*tau[i]))+1
end

#truncation = ceil(maximum(truncation))
truncation = 3
intensities = []
theta6=[];#zeros(truncation, N-length(TIRBundleIndexes))
theta7=[];#zeros(truncation, N-length(TIRBundleIndexes))
theta8 = []
theta9 = []
theta10 = []
theta11 = []
phi=[]
L = []
alpha = []
returnPos = {}
rho34 = []
rho32 = []
Rsub_bundle = zeros(1,N);#TIRBundleIndexes[1]-1)
Asub_bundle = zeros(1,N);#TIRBundleIndexes[1]-1)
Aglass_bundle = zeros(1,N);#TIRBundleIndexes[1]-1)
# Loop thru all theta1"s which don"t correspond to bundles undergoing
# TIR. This loop includes a single bounce out of the substrate

for i = 1:y
#       j = 1
      alpha = 90-theta1[i]
      for  j = 1:truncation
          intensities[j,i] = (1-rho12[i])*(1-rho23[i])*(tau[i])*(rho21[i] * rho23[i] * (tau[i]^2))^(j-1)

          # short position vector finding; yields a theta6           
          Pold = [r2*cosd[alpha], r2*sind[alpha]]
          beam = [distance[i]*cosd[180+alpha+theta2[i]], distance[i] * sind[180+alpha+theta2[i]]]
          Pnew = Pold + beam
          alpha = atan2d[Pnew[2],Pnew[1]]
          outsideAngle = theta4[i] + alpha;#180-theta4[i]-alpha()
          
          # beam moving left
          if outsideAngle .< 90 
              theta6[j,i] = 90-outsideAngle
              phi[j,i] = 90 - theta6[j,i]
              L[j,i] = Pnew[1] - Pnew[2]*tand[theta6[j,i]]
              theta7[j,i] = asind[(L[j,i]*sind[phi[j,i]])/r1]
              returnPos(i,j) = [r1*cosd[180 - phi[j,i] - theta7[j,i]],r1*sind[180 - phi[j,i] - theta7[j,i]]] #[r1*cosd[180 - (beta+2*theta6[j,i]) - theta7[j,i]],r1*sind[180 - (beta+2*theta6[j,i]) - theta7[j,i]]]
          
          # beam moving right
          else()                   
              theta6[j,i] = outsideAngle-90
              phi[j,i] = 90 + theta6[j,i]
              L[j,i] = Pnew[1] + Pnew[2]*tand[theta6[j,i]]
              theta7[j,i] = asind[(L[j,i]*sind[phi[j,i]])/r1]
              returnPos(i,j) = [r1*cosd[180 - phi[j,i] - theta7[j,i]],r1*sind[180 - phi[j,i] - theta7[j,i]]]
          end
          
          #long position vector finding
          Pold = Pnew
          beam = [distance[i] * cosd[alpha-theta3[i]], distance[i] * sind[alpha-theta3[i]]]
          Pnew = Pold + beam
          alpha = atan2d[Pnew[2],Pnew[1]]
          #Pold = Pnew
          
          #Calculate theta8
          theta8[j,i] = asind[(n3/n2)*sind[theta7[j,i]]]
          
          #Calculate theta9
          theta9[j,i] = asind[(r1/r2)*sind[theta8[j,i]]]
          
          #Calculate theta10
          theta10[j,i] = asind[(n2/n1)*sind[theta9[j,i]]]
          
          #Calculate theta11
          theta11[j,i] = asind[(n2/n3)*sind[theta8[j,i]]]
          
          #FINDING REFLECTIVITY RHO34
          
          term1 = ((n4^2) - (k4^2)-(n3^2) * (sind[theta6[j,i]])^2)^2
          term2 = 4 * (n4^2) * (k4^2)
          term3 = sqrt(term1)

          psquared = 0.5 * ((sqrt(term1 + term2)) + term3)
          qsquared = 0.5 * ((sqrt(term1 + term2)) - term3)

          p = sqrt(psquared)
          q = sqrt(qsquared)

          # Calculate reflectivity

          numterm1 = (n3*cosd[theta6[j,i]]-p)^2
          denomterm1 = (n3*cosd[theta6[j,i]]+p)^2

          numterm2 = (p - n3*sind[theta6[j,i]]*tand[theta6[j,i]])^2
          denomterm2 = (p + n3*sind[theta6[j,i]]*tand[theta6[j,i]])^2

          leftPortion = 0.5 * ((numterm1 + q^2)/(denomterm1 + q^2))
          rightPortion = 1 + ((numterm2 + q^2)/(denomterm2 + q^2))

          rho34[j,i] =  leftPortion * rightPortion
          
          #FINDING REFLECTIVITY RHO32
          rho32[j,i] = 0.5 * ( ((tand[theta7[j,i] - theta8[j,i]])^2)/((tand[theta7[j,i] + theta8[j,i]])^2)...
          + ((sind[theta7[j,i] - theta8[j,i]])^2)/((sind[theta7[j,i] + theta8[j,i]])^2))
      
          #FINDING TAU1 FOR BEAMS ON WAY OUT OF SHELL
          distance1[j,i] = r2*cosd[theta9[j,i]] - r1*cosd[theta8[j,i]]
          
          tau1[j,i] = exp(-k * distance1[j,i])
          
          #Finding rho21_o (reflectivity for beam off substrate onto top
          #surface of shell:
          rho21_o[j,i] = 0.5 * ( ((tand[theta9[j,i] - theta10[j,i]])^2)/((tand[theta9[j,i] + theta10[j,i]])^2)...
          + ((sind[theta9[j,i] - theta10[j,i]])^2)/((sind[theta9[j,i] + theta10[j,i]])^2))
      
      
          #Finding rho23_o [reflectivity for beams bouncing on inner
          #interface of shell after substrate bounce]
          rho23_o[j,i] = 0.5 * ( ((tand[theta8[j,i] - theta11[j,i]])^2)/((tand[theta8[j,i] + theta11[j,i]])^2)...
          + ((sind[theta8[j,i] - theta11[j,i]])^2)/((sind[theta8[j,i] + theta11[j,i]])^2))
          
          #Intensities after bouncing off of the substrate & absorbed by
          #the substrate
          
          
          while intensities[j,i] .> .00001
              Asub[j,i] = intensities[j,i]*(1-rho34[j,i])
              
              Asub_bundle[i] = Asub_bundle[i] + Asub[j,i]

              intensities[j,i] = intensities[j,i] * rho34[j,i]

              #Calculate the beams sent back out [reflectance] for every i,j
              #beam that enters the shell & bounces off of the substrate

              Rsub[j,i] = intensities[j,i]*((1 - rho21_o[j,i])*(1-rho32[j,i])*tau1[j,i])/ (1 - rho21_o[j,i]*rho23_o[j,i]*(tau1[j,i])^2)
              
              Rsub_bundle[i] = Rsub_bundle[i] + Rsub[j,i]

              #Calculate the beams sent back to the substrate after bouncing
              #off of the interior of the shell

              BackToSub[j,i] = intensities[j,i]* (rho32[j,i] + (rho21_o[j,i] * (tau1[j,i])^2)*(1 - rho23_o[j,i] - rho32[j,i]))...
              /(1 - rho21_o[j,i]*rho23_o[j,i]*(tau1[j,i])^2)


              #Calculate amount absorbed by glass on way back out after
              #substrate bounce

              Aglass[j,i] = intensities[j,i]-Rsub[j,i]-BackToSub[j,i]
              
              Aglass_bundle[i] = Aglass_bundle[i] + Aglass[j,i]


              # Theses are the intensities sent back to the substrate for each
              # transmission [j,i]. Reset & continue the while loop until
              # attenuated.
              intensities[j,i] = BackToSub[j,i]
          
          end
          
           
       end
        
end

#Compute total apparent properties & check to see if conservation of
#energy is satisfied
for i = 1:y
    
    a[i] = Rsub_bundle[i] *((i/N)*r2) / sumDENOM
    
    b[i] = Asub_bundle[i] *((i/N)*r2) / sumDENOM
    
    c[i] = Aglass_bundle[i] *((i/N)*r2) / sumDENOM
    
end

Rsubtot = sum(a)#length(Rsub_bundle)

Asubtot = sum(b)#length(Asub_bundle)

Aglasstot = sum(c)#length(Aglass_bundle)

# Rsubtot = sum(Rsub_bundle)/(N-1)#length(Rsub_bundle)
# 
# Asubtot = sum(Asub_bundle)/(N-1)#length(Asub_bundle)
# 
# Aglasstot = sum(Aglass_bundle)/(N-1)#length(Aglass_bundle)



check1 = Rtot + Atot + Rsubtot + Asubtot + Aglasstot

check2 = Rsubtot+Asubtot +Aglasstot

check3 = Ttot - check2


# # Plot Total Reflectance vs Wavelength
#  PLOT2 = plot(lambda*(1e6),Rtotdata+RsubtotData)
#  #plot(lambda*(10e6),Rdata)
# xlabel ("lambda [wavelength]in microns','fontweight','bold'); ylabel ('Total Reflectance"...
#     ,"fontweight', 'bold")
# ax = ancestor[PLOT2, "axes"]
# ax.YAxis.Exponent = 0
# ytickformat["#.4f"]
# title("Total Apparent Reflectance vs. Incident Wavelength")
toc()

# # Plot Total Substrate Absorptance vs Wavelength
#  PLOT2 = plot(lambda*(1e6),AsubtotData)
#  #plot(lambda*(10e6),Rdata)
# xlabel ("lambda [wavelength]in microns','fontweight','bold'); ylabel ('Total Absorbed By Substrate"...
#     ,"fontweight', 'bold")
# ax = ancestor[PLOT2, "axes"]
# ax.YAxis.Exponent = 0
# ytickformat["#.4f"]
# title("Total Absorbed by Substrate vs. Incident Wavelength")

#plot(theta6[1,:],rho34[1,:])
# 
# 
# 
