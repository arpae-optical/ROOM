# for rectangular walls with side lengths l and w, and height h:
#  1__
# 2| T|2
# 2|B |2 W
#    1
# let us find all view factors, and then calculate apparent emissivity:

include("D39.jl")

function recwall(l,w,h,eps::Vector)
epa=Vector{Float64}(undef,length(eps));

FT1=D39(l,w,h);
FT2=D39(w,l,h);
FBW=2*FT1+2*FT2;
for ii=1:length(eps)
    if eps[ii]==0
        epa[ii]=0;
    else
        epa[ii]=(1-FBW/2)/(1/eps[ii]-FBW/2*(1/eps[ii]-1));
    end
end
end

function recwall(l,w,h,eps::Float64)
FT1=D39(l,w,h);
FT2=D39(w,l,h);
FBW=2*FT1+2*FT2;
    if eps==0
        epa=0;
    else
        epa=(1-FBW/2)/(1/eps-FBW/2*(1/eps-1));
    end
end