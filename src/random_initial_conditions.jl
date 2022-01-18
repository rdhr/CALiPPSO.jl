#= 
BE SURE TO INCLUDE/CALL THE SCRIPT `packing_type.jl` BEFORE INCLUDING THIS FILE, BECAUSE IN SUCH FILE THERE ARE SOME FUNCTIONS AND LIBRARIES
THAT ARE USED HERE
=#
using Distributions, StaticArrays   
#################################################
# Functions to obtain particle radii from packing fraction and other simple geometric properties
################################################

"Compute the radius from a given packing fraction 'ϕ', in arbitrary dimensions, 'd'."
function radius_from_phi(d::Int64, ϕ::Real, N::Int64, L::Real=1.0)::Float64
    #d: dimension; ϕ: packing fraction; N: number of particles; L: size of box
    Vd = volume_d_ball(d)
    return L*(ϕ/(N*Vd))^(1/d)
end
radius_from_phi(3, 0.64, 1000)

#################################################
# Functions to generate a random configuration that is NOT jammed, but can be used as initial condition for iLP
################################################

"Asses whether two particles of radius 'r' and centers 'X1' and 'X2' are overlapping. "
function are_overlapping(X1::SVector{d, PeriodicNumber{T}}, X2::SVector{d, PeriodicNumber{T}}, r::T)::Bool where {d, T<:Real}
    if norm(X1-X2)< 2*r
        return true
    else
        return false
    end
end
are_overlapping(PeriodicVector(rand(4)), PeriodicVector(rand(4)), 0.0)


"""
    generate_random_configuration(d::Int64, N::Int64, ϕ::T, L::T=1.0; max_tries::Int64=5000 )

Generate 'N' random centers of monodisperse particles of radius 'r' and in 'd' dimensions, 
without any overlaps.

The output is the radius (that corresponds to the packing fraction 'ϕ' used as input), and 
the vector containing the centers (each as a SVector{d, PeriodicNumber} type). 
Each center is placed at uniformly random in all space, and when an overlap is detected a 
new random position is drawn. 'max_tries' (default 5000) attempts are tried for each particle 
and when this bound surpassed an error is thrown, with the index of the center that was not
created.
"""
function generate_random_configuration(d::Int64, N::Int64, ϕ::T, L::T=1.0; max_tries::Int64=5000 ) where {T<:Real}
    r = radius_from_phi(d, ϕ, N, L)
    Xs_distr = Uniform(0, L)
    centers = Vector{SVector{d, PeriodicNumber{T}}}(undef, N)

    centers[1] = PeriodicVector(rand(Xs_distr, d), L)

    for i in 2:N
        new_center = PeriodicVector(rand(Xs_distr, d), L)
        overlaps = are_overlapping.(centers[1:i-1], [new_center], r)
        c=0
        
        while any(overlaps)
            new_center = PeriodicVector(rand(Xs_distr, d), L)
            c+=1
            overlaps = are_overlapping.(centers[1:i-1], [new_center], r)

            (c>max_tries) && error("Could not assign a new center to particle ", i, "    after ", max_tries, " tries")
        end
        centers[i] = new_center
    end
    return r, centers
end
generate_random_configuration(3, 10, 0.2)