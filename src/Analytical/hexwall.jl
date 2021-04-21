# for a hex with center to center spacing l, and height h, side length s:
#  2/\3
# 1|BT|4 W
#  2\/3
# let us find all view factors, and then calculate apparent emissivity.

include("C16.jl")
include("D38.jl")
function fhex(l_h)
    h = 1
    l = l_h * h
    s = l / sqrt(3)
    AW = 6 * s * h
    AB = 3 / 2 * s * l

    F12 = C16(s, s, h, 2 * pi / 3)

    F1p3p = C16(2 * s, 2 * s, h, pi / 3)
    F1p3m = C16(1 * s, 2 * s, h, pi / 3)
    F1p3 = F1p3p - F1p3m

    F1m3p = C16(2 * s, 1 * s, h, pi / 3)
    F1m3m = C16(1 * s, 1 * s, h, pi / 3)
    F1m3 = F1m3p - F1m3m

    F13p = 2 * F1p3
    F13m = F1m3

    F13 = F13p - F13m

    F14 = D38(s, h, l)

    FWW = 2 * F12 + 2 * F13 + F14
    FWB = (1 - FWW) / 2
    fhex = AW * FWB / AB
end

function hexwall(l_h, eps::AbstractVector)
    epa = eps .* 0
    FBW = fhex(l_h)

    for ii in 1:length(eps)
        if eps(ii) == 0
            epa(ii) = 0
        else
            epa(ii) = (1 - FBW / 2) / (1 / eps(ii) - FBW / 2 * (1 / eps(ii) - 1))
        end
    end
end
