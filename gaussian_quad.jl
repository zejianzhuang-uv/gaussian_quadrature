


import FastGaussQuadrature

function gaussian_quad(fcn::Function, a::Float64, b::Float64; npt::Int64=10)
    a *= 1.0
    b *= 1.0
    EPS = 5e-2
    Δ = 1.
    last_inte = 220e0
    i::Int64 = npt

    while abs(Δ) >= EPS
        inte = 0.
        xw = FastGaussQuadrature.gausslegendre(i)#xw[i]
        x, w = xw
        xp, wp = (b - a)/2 .* x .+ (a+b)/2, (b-a)/2 .* w

        for (ii, xxp) in enumerate(xp)
            inte += wp[ii] * fcn(xxp)
        end
        Δ = (inte - last_inte) / inte
        last_inte = inte
        i += 10
        if i == 100
            break
        end
        # if length(x) == 50
        #     #println("Reaching the maximum order of Gaussian quadrature")
        #     break
        # end
    end
    return last_inte
end