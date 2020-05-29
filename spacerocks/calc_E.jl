#function calc_E(e, M)
#
#    output = Array{Float64}(undef, 1, length(e))
#
#    for idx in 1:length(e)
#
#        E = M[idx]
#
#        for idx in 1:10
#            E = E - (E - M[idx] - e[idx] * sin(E)) / (1 - e[idx] * cos(E))
#        end
#
#        output[idx] = E
#
#    end
#
#    return output
#
#end

function calc_E(e, M)
    output = Array{Float64}(undef, 1, length(e))
    for idx in 1:length(e)

        e_i = e[idx]
        M_i = rem2pi(M[idx], RoundNearest)

        if iszero(M_i) || iszero(e_i)
            output[idx] = mod2pi(M_i)

        else
            α = (3 * pi^2 + 1.6 * (pi^2 - pi * abs(M_i))/(1 + e_i))/(pi^2 - 6)
            d = 3 * (1 - e_i) + α * e_i
            q = 2 * α * d * (1 - e_i) - M_i^2
            r = 3 * α * d * (d - 1 + e_i) * M_i + M_i^3
            w = cbrt(abs2(abs(r) + sqrt(q^3 + r^2)))
            E1 = (2 * r * w / @evalpoly(w, q^2, q, 1) + M_i)/d
            f2 = e_i * sin(E1)
            f3 = e_i * cos(E1)
            f0 = E1 - f2 - M_i
            f1 = 1 - f3
            δ3 = -f0 / (f1 - f0 * f2 / (2 * f1))
            δ4 = -f0 / @evalpoly(δ3, f1, f2 / 2, f3 / 6)
            δ5 = -f0 / @evalpoly(δ4, f1, f2 / 2, f3 / 6, - f2 / 24)
            output[idx] = mod2pi(E1 + δ5)

        end

    end
    return output
end
