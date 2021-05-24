## MCCylCyl.m
# Given geometry, number of bundles, and number of reflections to track,
# returns the view factor of an array of cylindrical structures to
# themselves, given perfect reflection at the boundary. 
# arguments: N, L/H, D/H
# N: number of bundles
# L/H: edge to edge distance of hex divided by the height
# D/H: diameter of occluding cylinder divided by the height
# n: number of reflections (default: 5)
# Geometry is of form:
#   |Y             |Z
#  ___            ___
# / O \  _X     || O || _X
# \___/         ||___||
#
function MCCylCyl(N::Int64,L_H::Float64,D_H::Float64,nr::Int64=5)

# First, let's initialize three booleans to keep track of our bundles:

ray = ones(Bool,(N,1))  # Currently "active" rays
out = zeros(Bool,(N,1)) # Terminated rays out of bounds (top/bottom planes)
faa = out               # Terminated rays on the target (occlusion)

# We will also initialize some more arrays for memory happiness later:
L2  = zeros(N) # Occlusion helper variable
tc  = L2
d2  = tc
t1c = d2
t1  = t1c
cob = faa

# Geometry of our unit cell:

H   = 1         # Height of walls
L   = H*L_H     # Edge-to-edge distance of unit cell borders
s   = L/sqrt(3) # Edge length of unit cell borders
R   = H*D_H/2   # Radius of cylinder in middle of unit cell
R2  = R^2       # Radius of cylinder in middle of unit cell, squared

# Now let's establish our monte carlo RNG:

MCW = rand(N) # Monte-Carlo RNG for psi distribution of bundles
MCZ = rand(N) # Monte-Carlo RNG for z distribution of bundles
MCP = rand(N) # Monte-Carlo RNG for phi distribution of bundles
MCQ = rand(N) # Monte-Carlo RNG for theta distribution of bundles

# Converting our MC origin variables to XYZ coordinates

OW  = 2*pi*MCW      # Origin of bundles, psi (angle from x-axis)
OX  = R*cos.(OW)    # Origin of bundles
OY  = R*sin.(OW)    # Origin of bundles
OZ  = (-1/2 .+ MCZ)*H  # Origin of bundles

# Converting our MC direction variables to XYZ directions

NP  = 2*pi*MCP                                        # Phi of bundles
NQ  = asin.(sqrt.(MCQ))                               # Theta of bundles
NX  = cos.(OW).*cos.(NQ)+sin.(OW).*sin.(NQ).*sin.(NP) # X direction of bundles
NY  = sin.(OW).*cos.(NQ)-cos.(OW).*sin.(NQ).*sin.(NP) # Y direction of bundles
NZ  = sin.(NQ).*cos.(NP)                              # Z direction of bundles

# Now that we have established the location and direction of all bundles
# from the central occlusion, we wish to positively identify which wall the
# bundle will intercept, rather than testing all six until one is found. We
# know that until the first reflection, the bundle has two possibilities:
# intercept a perfectly reflecting wall, or intercept the top/bottom. We
# are only interested in the first case. Inscribing a circle around the
# unit cell, we can identify the wall it intercepts first by the angle from
# the origin to the location on the bounding circle.
# The circle represented by x^2+y^2=s^2, and our bundle's location given by
# x=xo+nx*t,y=yo+ny*t,z=zo+nz*t, we can solve for t:

a   = NX.^2+NY.^2
b   = 2*NX.*OX+2*NY.*OY
c   = OX.^2+OY.^2 .- s^2

ti  = (-b+sqrt.(b.^2-4*a.*c))./(2*a)

# The location on that circle, then, is:

DX  = OX+ti.*NX
DY  = OY+ti.*NY

# We can find the angle from the x axis to this point on the circle by:

the = floor.(((DY .> 0).*2*pi .+ sign.(DY).*acos.(DX/s))/pi/3) .+ pi/6

# This lets us solve for the normal vector of the plane intercepted:

nx  = -cos.(the)
ny  = -sin.(the)

# Solving for the new interception distance (to the plane, not the bounding
# circle)

ti  = (L/2 .- (OX.*nx+OY.*ny))./(NX.*nx+NY.*ny)

# We can use this to find the reflection location and direction.

OX  = OX+NX.*ti
OY  = OY+NY.*ti
OZ  = OZ+NZ.*ti

NX  = NX+2*NX.*nx
NY  = NY+2*NY.*ny

# If the reflection location is outside of the bounds of our unit cell, it
# is out of bounds, otherwise we continue.

out = Bool(abs.(OZ) .> H/2)

### This line is giving me trouble- do you know how I would compare two boolean arrays?
### I want ray::Bool to be equal to .false. anywhere out::Bool OR faa::Bool are .true.
### and .true. anywhere neither out nor faa are .true. My error I'm getting is:
## 'no method matching !(::BitArray{1})

ray = !out .| !faa

# Until we terminate based on number of reflections,
for ii=1:nr
    # Check for interception with occlusion:
    L2(ray)  = OX(ray).^2+OY(ray).^2        # Distance from center of occlusion
    tc(ray)  = OX(ray).*NX(ray)+OY(ray).*NY(ray)      # Distance to closest approach
    d2(ray)  = L2(ray)-tc(ray).^2 # Distance of closest approach 
    cob(ray) = d2(ray)<R2         # Check Occlusion Boolean: If closest
                                  # approach is less than radius, occlude:
                                   
    t1c(cob) = sqrt.(R2-d2(cob))  # Distance from ca to radius
    t1(cob)  = tc(cob)-t1c(cob)   # Distance from origin to intercept
    if dbb # Only need to update occluded x&y location if debugging
        OX(cob) = OX(cob)+NX(cob).*t1(cob)
        OY(cob) = OY(cob)+NY(cob).*t1(cob)
    end
    OZ(cob)  = OZ(cob)+NZ(cob).*t1(cob)# Z coord of interception
    out(cob) = abs(OZ(cob))>H/2        # If out of bounds, mark
    faa(cob) = true-out(cob)           # Otherwise, occluded
    ray      = ray&~cob       # Update active bundles, and cont.
    
    # If our ray has not been occluded by the cylinder, we will check for
    # interception and do all calculations same as the original bundle:
    a(ray)   = NX(ray).^2+NY(ray).^2
    b(ray)   = 2*NX(ray).*OX(ray)+2*NY(ray).*OY(ray)
    c(ray)   = OX(ray).^2+OY(ray).^2-s^2

    ti(ray)  = (-b(ray)+sqrt.(b(ray).^2-4*a(ray).*c(ray)))./(2*a(ray))

    DX(ray)  = OX(ray)+ti(ray).*NX(ray)
    DY(ray)  = OY(ray)+ti(ray).*NY(ray)

    the(ray) = (floor(((DY(ray)>0).*2*pi+
                 sign.(DY(ray)).*acos.(DX(ray)/s))/pi/3)+pi/6)

    nx(ray)  = -cos.(the(ray))
    ny(ray)  = -sin.(the(ray))

    ti(ray)  = ((L/2-(OX(ray).*nx(ray)+OY(ray).*ny(ray)))/
                 (NX(ray).*nx(ray)+NY(ray).*ny(ray)))
             
    OX(ray)  = OX(ray)+NX(ray).*ti(ray)
    OY(ray)  = OY(ray)+NY(ray).*ti(ray)
    OZ(ray)  = OZ(ray)+NZ(ray).*ti(ray)

    NX(ray)  = NX(ray)+2*NX(ray).*nx(ray)
    NY(ray)  = NY(ray)+2*NY(ray).*ny(ray)

    out(ray) = abs(OZ(ray))>H/2
    ray(ray) = ~out(ray)|~faa(ray)

end

if any(faa)
    FCC = sum(faa)/N
    if sum(ray)~=0
        print(string(sum(ray)," unterminated rays!"))
    end
elseif sum(ray)==N
    FCC = 0
    print("No terminated rays!")
else
    FCC = 0
    print("No rays on target!")
end

return FCC

end