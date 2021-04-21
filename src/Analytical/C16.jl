# C16.jl returns the view factor between two adjacent rectangles, of
# shared length h, width of surface 2 a, width of surface 1 b, and angle
# between them phi. From UT Austin, C-16, thermalradiation.net

# TODO switch to symbolics.jl. should be as simple as changing @vars -> @variables and replacing sympy with `using Symbolics`
using SymPy
function C16(a, b, h, ϕ)
    A = a / h
    B = b / h
    C = A^2 + B^2 - 2A * B * cos(ϕ)
    D = sqrt(1 + A^2 * sin(ϕ)^2)
    @vars x

    # then, breaking up the giant term into several smaller ones:
    t21 =
        A * B * sin(ϕ) +
        (π / 2 - ϕ) * (A^2 + B^2) +
        B^2 * atan((A - B * cos(ϕ)) / (B * sin(ϕ))) +
        A^2 * atan((B - A * cos(ϕ)) / (A * sin(ϕ)))
    t22 =
        (2 / sin(ϕ)^2 - 1) * log((1 + A^2) * (1 + B^2) / (1 + C)) +
        B^2 * log((B^2 * (1 + C)) / (C * (1 + B^2))) +
        A^2 * log((A^2 * (1 + A^2)^(cos(2ϕ))) / (C * (1 + C)^(cos(2ϕ))))
    t23 = 1 / π * atan(1 / B) + A / π / B * atan(1 / A) - sqrt(C) / π / B * atan(1 / sqrt(C))
    t24 = atan(A * cos(ϕ) / D) + atan((B - A * cos(ϕ)) / D)
    t2h(x) = sqrt(1 + x^2 * sin(ϕ)^2)
    t25(x) = t2h(x) * (atan(x * cos(ϕ) / t2h(x)) + atan((A - x * cos(ϕ)) / t2h(x)))
    F12 =
        -sin(2ϕ) / 4π / B * t21 +
        sin(ϕ)^2 / 4π / B * t22 +
        t23 +
        sin(ϕ) * sin(2ϕ) / 2π / B * A * D * t24 +
        cos(ϕ) / π / B * N(integrate(t25(x), (x, 0, B)))
end

