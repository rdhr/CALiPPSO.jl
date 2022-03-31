using Distributions

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

"Compute the radii of a *bi*disperse configuration of (N1, N2) particles with radii (R1, R2), from a given packing fraction 'ϕ', in arbitrary dimensions, 'd'. It is assumed that R2/R1=c."
function radii_from_phi(d::Int64, ϕ::Real, N1::Int64, N2::Int64=N1, c::T=1.4, L::T=1.0) where T<:AbstractFloat
    #d: dimension; ϕ: packing fraction; N: number of particles; L: size of box
    Vd = volume_d_ball(d)
    den = Vd*(N1 + N2*c^d)
    R1 = L*(ϕ/den)^(1/d)
    return R1, c*R1
end
radii_from_phi(2, 0.84, 1000)



#################################################
# Functions to generate a random configuration that is NOT jammed, but can be used as initial condition for CALiPPSO
################################################

"Asses whether two particles of radius 'r' and centers 'X1' and 'X2' are overlapping. "
function are_overlapping(X1::SVector{d, PeriodicNumber{T}}, X2::SVector{d, PeriodicNumber{T}}, r::T)::Bool where {d, T<:AbstractFloat}
    if norm(X1-X2)< 2*r
        return true
    else
        return false
    end
end
are_overlapping(PeriodicVector(rand(4)), PeriodicVector(rand(4)), 0.0)

"Asses whether two particles whose centers are 'X1' and 'X2', and with radii 'r1' and 'r2' are overlapping. "
function are_overlapping(X1::SVector{d, PeriodicNumber{T}}, X2::SVector{d, PeriodicNumber{T}}, r1::T, r2::T)::Bool where {d, T<:AbstractFloat}
    if norm(X1-X2)< r1+r2
        return true
    else
        return false
    end
end
are_overlapping(PeriodicVector(rand(4)), PeriodicVector(rand(4)), 0.0, 0.0)


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
function generate_random_configuration(d::Int64, N::Int64, ϕ::T, L::T=1.0; max_tries::Int64=5000 ) where {T<:AbstractFloat}
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


function generate_random_configuration(d::Int64, N1::Int64, R1::T, N2::Int64, R2::T, L::T=1.0; max_tries::Int64=5000, verbose::Bool=false ) where {T<:AbstractFloat}
    N = N1+N2; Rs = [R1*ones(N1); R2*ones(N2)]

    verbose && println("Trying to generate a d-", d, " dimensional configuration of N = ", N, " spheres with packing fraction φ= ", packing_fraction(d, R1, N1, R2, N2, L))

    Xs_distr = Uniform(0, L)
    centers = Vector{SVector{d, PeriodicNumber{T}}}(undef, N)
    # create array of radii
    centers[1] = PeriodicVector(rand(Xs_distr, d), L)

    for i in 2:N
        new_center = PeriodicVector(rand(Xs_distr, d), L)
        overlaps = [are_overlapping(centers[j], new_center, Rs[j], Rs[i]) for j in 1:i-1]
        # overlaps = are_overlapping.(centers[1:i-1], [new_center], r)
        c=0
        
        while any(overlaps)
            new_center = PeriodicVector(rand(Xs_distr, d), L)
            c+=1
            # overlaps = are_overlapping.(centers[1:i-1], [new_center], r)
            overlaps = [are_overlapping(centers[j], new_center, Rs[j], Rs[i]) for j in 1:i-1]

            (c>max_tries) && error("Could not assign a new center to particle ", i, "    after ", max_tries, " tries")
        end
        centers[i] = new_center
    end
    return Rs, centers
end
generate_random_configuration(3, 10, 0.02, 5, 0.01)

function generate_random_configuration(d::Int64, N1::Int64, N2::Int64, ϕ::T, c::T=1.4, L::T=1.0; max_tries::Int64=5000, verbose::Bool=false ) where {T<:AbstractFloat}
    R1, R2 = radii_from_phi(d, ϕ, N1, N2, c, L)
    Rs, centers = generate_random_configuration(d, N1, R1, N2, R2, L; max_tries=max_tries, verbose=verbose)
    return R1, R2, centers
end
generate_random_configuration(3, 50, 50, 0.2, 1.0)


function generate_random_configuration(d::Int64, Rs::Vector{T}, L::T=1.0; max_tries::Int64=5000, verbose::Bool=false ) where {T<:AbstractFloat}
    N = length(Rs)
    φ = packing_fraction(d, Rs, L)
    verbose && println("Trying to generate a d-", d, " dimensional configuration of N = ", N, " spheres with packing fraction φ= ", φ)

    Xs_distr = Uniform(0, L)
    centers = Vector{SVector{d, PeriodicNumber{T}}}(undef, N)
    # create array of radii
    centers[1] = PeriodicVector(rand(Xs_distr, d), L)

    for i in 2:N
        new_center = PeriodicVector(rand(Xs_distr, d), L)
        overlaps = [are_overlapping(centers[j], new_center, Rs[j], Rs[i]) for j in 1:i-1]
        # overlaps = are_overlapping.(centers[1:i-1], [new_center], r)
        c=0
        
        while any(overlaps)
            new_center = PeriodicVector(rand(Xs_distr, d), L)
            c+=1
            # overlaps = are_overlapping.(centers[1:i-1], [new_center], r)
            overlaps = [are_overlapping(centers[j], new_center, Rs[j], Rs[i]) for j in 1:i-1]

            (c>max_tries) && error("Could not assign a new center to particle ", i, "    after ", max_tries, " tries")
        end
        centers[i] = new_center
    end
    return φ, centers
end
generate_random_configuration(2, 0.02*rand(10))