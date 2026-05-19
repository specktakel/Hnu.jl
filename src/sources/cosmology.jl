module Cosmology

using Integrals: IntegralProblem, solve, QuadGKJL
include("./cosmo_units.jl")
using Unitful
using .CosmoUnits
using Roots

Om = 0.3
Ol = 0.7
H0 = 70u"km/s/Mpc"
c = 3e5u"km/s"
DH = c / H0

function E(z::Real)
    Omp = Om * (1 + z)^3
    return sqrt(Omp + Ol)
end

function hubble_factor(z::Real)
    return H0 * E(z)
end

function comoving_distance(z::Real)
    function scale(z::Real)
        return 1 / E(z)
    end

    sol = solve(IntegralProblem(scale, (0, z)), QuadGKJL())
    return sol.u * DH
end

function luminosity_distance(z::Real)
    return (1 + z) * comoving_distance(z)
end

function redshift(dL::Unitful.Length)
    find_zero(x -> ustrip(u"Mpc", dL) - ustrip(u"Mpc", luminosity_distance(x)), 0.2)
end

end