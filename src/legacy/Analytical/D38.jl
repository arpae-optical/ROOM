"""
Analytical View Factor from Modest, M. "Radiative Heat Transfer", 3rd Ed.

Appendix D. 38: Identical, parallel, directly opposed rectangles.

- a: width of rectangle
- b: depth of rectangle
- c: distance between rectangles
"""
function D38(a, b, c)
    X = a / c
    Y = b / c
    T1 = 1 / 2 * log((1 + X^2) * (1 + Y^2) / (1 + X^2 + Y^2))
    T2 = X * sqrt(1 + Y^2) * atan(X / sqrt(1 + Y^2))
    T3 = Y * sqrt(1 + X^2) * atan(Y / sqrt(1 + X^2))
    T4 = X * atan(X) + Y * atan(Y)
    F12 = 2 / pi / X / Y * (T1 + T2 + T3 - T4)
end
