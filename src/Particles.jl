include("functions-for-periodic-vectors.jl")
using LinearAlgebra, Statistics, SpecialFunctions

const default_tol_force_equilibrium = 1e-12

"Define an abstract type for all the type of particles that will be used."
abstract type AbstractParticle{d, T} end # To define an abstract type for all type of particles to be used

#########################################################################################################
#########################################################################################################
### 'MonoParticle' type and related constructors and functions. Type meant to be used in monodisperse packings
#########################################################################################################
#########################################################################################################

#= 
NB:We are including the docstring of all the fields of 'MonoParticle', following the recommendations in https://docs.julialang.org/en/v1/manual/documentation/ . However, such usage produces an error in parametric struct's, just as 'MonoParticle' is (see, e.g. https://github.com/JuliaLang/julia/issues/39410). So we are leaving such docstring, but also providing the description of each field in the docstring of the struct itself.
=#

""" 
    MonoParticle{d, T}

Struct used to construct a particle in a d-dimensional space. 

It's meant to be used as components of *monodisperse* packings. Hence, while many of its 
relevant properties are stored, its radius is instead stored in a field in 'MonoPacking' type.

The quantities this struct has are:
1. X: The particle's center (as a StaticVector{d, PeriodicNumber{T}}})
2. contact_vecs: The set of contact vectors with the particle's neighbours (as a Vector containing SVector{d,T})
3. forces: The magnitude of contact forces (Vector{T})
4. neighbours: The indices of the particle's neighbours (as Vector{Int64})

Naturally, 'contact_vecs', 'forces', and 'neighbours' have the same length.

See also: [`MonoPacking`](@ref).
"""
mutable struct MonoParticle{d, T}  <: AbstractParticle{d, T} # a struct defining a particle in a d-dimensional MONODISPERSE packing
    "Position of particle's center (SVector{d, PeriodicNumber{T}}"
    X::SVector{d, PeriodicNumber{T}}   # particle's centre
    "Vector of contact vectors, pointing towards X. Array of {SVector{d,T}}"
    contact_vecs::Vector{SVector{d,T}} # array of contact vectors (should point towards X)
    "Vector with the magnitudes of the contact forces"
    forces::Vector{T}  # list of forces magnitudes
    "List of particles' indices with which  the particle is in contact"
    neighbours::Vector{Int64}  # list of the indices, corresponding to particles in contact
end

#### Format how a `MonoParticle` is shown in IO
function Base.show(io::IO, particle::MonoParticle{d,T}) where {d, T<:AbstractFloat}
    x = value.(particle.X)
    neighs = particle.neighbours
    forces = particle.forces

    println(io, d, "d MonoParticle\t centred at:  ", x)
    println(io, "neighbours: ", neighs)
    println(io, "contact forces: ", forces)
end

### Now we define some convenient constructors (and related functions) for 'MonoParticle'

#This function is suggested in the 'StaticArrays' documentation (https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/#Guide-1)
"Convert a matrix to a vector with SVector elements."
function svectors(M::Matrix{T}, ::Val{d}) where {T,d}
    size(M,1) == d || error("sizes mismatch")
    isbitstype(T) || error("use for bitstypes only")
    copy(reinterpret(SVector{d,T}, vec(M)))
end
svectors(reshape(collect(1:6), (2,3)), Val(2))


function MonoParticle(X::SVector{d, PeriodicNumber{T}}, contact_vecs::Matrix{T}, fs::Vector{T}, neighbours::Vector{Int64}) where {d, T<:AbstractFloat}

    dim, z = size(contact_vecs)
    (dim!=d) && error("There is a mismatch in the dimensionality!! Dimensions of particle's position: ", d, "\t dimensions of contact vectors: ", dim)

    if z==length(fs) && z==length(neighbours)
        return MonoParticle(X, svectors(contact_vecs, Val(d)), fs, neighbours)
    else
        error("Number of contacts mismatch! There are: ", z, " contact vectors \t", length(fs), " forces and \t ", length(neighbours), " contacts indices")
    end
end
MonoParticle(SVector{3}(PeriodicNumber.(rand(3), 0.5)), rand(3,10), rand(10), rand(1:100,10))

# The following two methods are useful if one's being forgetful (or a bit lazy) about defining 'X' as a SVector of 'PeriodicNumber' elements

function MonoParticle(X::Vector{PeriodicNumber{T}}, contact_vecs::Matrix{T}, fs::Vector{T}, neighbours::Vector{Int64}) where {T<:AbstractFloat}
    d = length(X)
    MonoParticle(SVector{d}(X), contact_vecs, fs, neighbours)
end
MonoParticle(PeriodicNumber.(rand(3), 0.5), rand(3,10), rand(10), rand(1:100,10))


MonoParticle(X::Vector{T}, L::T, contact_vecs::Matrix{T}, fs::Vector{T}, neighbours::Vector{Int64}) where {T<:AbstractFloat} = MonoParticle(PeriodicNumber.(X,L), contact_vecs, fs, neighbours)


#########################################################################################################
#########################################################################################################
### 'Particle' type and related constructors and functions. Type meant to be used in polydisperse packings
#########################################################################################################
#########################################################################################################

""" 
    Particle{d, T}
    This struct is used to construct a particle in a d-dimensional space. It's meant to be used as components of *bidisperse* or *polydisperse* packings. Hence, in this case also its radius is stored, along with other of its relevant properties are stored.
    The quantities this struct has are:
        1. X: The particle's center (as a StaticVector{d, PeriodicNumber{T}}})
        2. R: The particle's radius
        2. contact_vecs: The set of contact vectors with the particle's neighbours (as a Vector containing SVector{d,T})
        3. forces: The magnitude of contact forces (Vector{T})
        4. neighbours: The indices of the particle's neighbours (as Vector{Int64})
    Naturally, 'contact_vecs', 'forces', and 'neighbours' have the same length.
"""
mutable struct Particle{d, T}  <: AbstractParticle{d, T} # a struct defining a particle in a d-dimensional MONODISPERSE packing
    X::SVector{d, PeriodicNumber{T}}   # particle's centre
    R::T # particle's radius
    contact_vecs::Vector{SVector{d,T}} # array of contact vectors (should point towards X)
    forces::Vector{T}  # list of forces magnitudes
    neighbours::Vector{Int64}  # list of the indices, corresponding to particles in contact
end

#### Format how a `Particle` is shown in IO
function Base.show(io::IO, particle::Particle{d,T}) where {d, T<:AbstractFloat}
    x = value.(particle.X)
    neighs = particle.neighbours
    forces = particle.forces

    println(io, d, "d Particle\t centred at:  ", x, " \t and radius R = ", particle.R)
    println(io, "neighbours: ", neighs)
    println(io, "contact forces: ", forces)
end



function Particle(X::SVector{d, PeriodicNumber{T}}, R::T, contact_vecs::Matrix{T}, fs::Vector{T}, neighbours::Vector{Int64}) where {d, T<:AbstractFloat}

    z = size(contact_vecs, 2)

    if z==length(fs) && z==length(neighbours)
        Scvs = Vector{SVector{d,T}}(undef,z)
        for k in 1:z
            Scvs[k] = SVector{d}(contact_vecs[:, k])
        end

        return Particle(X, R, Scvs, fs, neighbours)
    else
        error("Number of contacts mismatch! There are: ", z, " contact vectors \t", length(fs), " forces and \t ", length(neighbours), " contacts indices")
    end
end


# The following two methods are useful if one's being forgetful (or a bit lazy) about defining 'X' as a SVector of 'PeriodicNumber' elements
function Particle(X::Vector{PeriodicNumber{T}}, R::T, contact_vecs::Matrix{T}, fs::Vector{T}, neighbours::Vector{Int64}) where {T<:AbstractFloat}

    d, z = size(contact_vecs)

    if z==length(fs) && z==length(neighbours)
        SX = SVector{d}(X)

        Scvs = Vector{SVector{d,T}}(undef,z)
        for k in 1:z
            Scvs[k] = SVector{d}(contact_vecs[:, k])
        end

        return Particle(SX, R, Scvs, fs, neighbours)
    else
        error("Number of contacts mismatch! There are: ", z, " contact vectors \t", length(fs), " forces and \t ", length(neighbours), " contacts indices")
    end
end
Particle(PeriodicNumber.(rand(3), 0.5), rand(), rand(3,10), rand(10), rand(1:100,10))


Particle(X::Vector{T}, R::T, L::T, contact_vecs::Matrix{T}, fs::Vector{T}, neighbours::Vector{Int64}) where {T<:AbstractFloat} = Particle(PeriodicNumber.(X,L), R, contact_vecs, fs, neighbours)


#########################################################################################################
#########################################################################################################
### Other functions; e.g. to compute volume
#########################################################################################################
#########################################################################################################

"""
    volume_d_ball(d::Int64, R::Real)

Compute the volume of a d-dimensional sphere of radius R (that defaults to 1)
"""
function volume_d_ball(d::Int64, R::Real=1)
    ((sqrt(π)*R)^d)/gamma(1.0 + 0.5d)
end

"""
    volume(P::Particle)

Compute the volume of 'P' (i.e. a hypersphere of d dimensions and radius R=P.R.
"""
volume(P::Particle) = volume_d_ball(length(P.X), P.R)

function check_for_overlaps(P1::Particle{d,T}, P2::Particle{d,T}, tolerance::T) where {d, T<:AbstractFloat}
    σij = P1.R + P2.R
    dist = norm(P1.X - P2.X)
    overlap=false; gap = 0.0
    if dist-σij < -tolerance
        gap = σji - dist
        overlap = true
    end
    return overlap, gap
end

# function check_for_overlaps(Xs::Vector{SVector{d, PeriodicNumber{T}}}, R::T, tolerance::T) where {d, T<:AbstractFloat}
#     σ = 2*R
#     N = length(Xs)
#     distances = distances_between_centers(Xs)
#     overlap=false; message = "No overlap is present"; particles = (0, 0)
#     for i in 1:N, j in i+1:N
#         if distances[i,j]-σ < -tolerance
#             gap = σ - distances[i,j]
#             Xi = value.(Xs[i])
#             Xj = value.(Xs[j])
#             message = "Overlap between particles $i and $j;\n Centers at: $Xi; \t $Xj. \n\t\tradius:  $R; \t Overlap = $gap ."

#             overlap = true
#             particles = (i, j)
#             break
#         end
#     end
#     return overlap, message, particles
# end

#########################################################################################################
#########################################################################################################
### Functions for generic (i.e. abstract) type of particle
#########################################################################################################
#########################################################################################################

"Return the coordination number (z) of 'P'."
get_coordination_number(P::AbstractParticle) = length(P.neighbours)
# get_coordination_number(P1); get_coordination_number(P3); 

"Return the norm of each of the contact vectors of 'P'."
norm_contact_vectors(P::AbstractParticle) = norm.(P.contact_vecs)
# norm_contact_vectors(P1); norm_contact_vectors(P3);

"""
Return whether 'P' is a rattler or not; that is, if it has less than d+1 contacts.

See also: [`is_stable_particle`](@ref).
"""
function is_a_rattler(P::AbstractParticle{d, T}) where {d, T<:AbstractFloat}
    z = get_coordination_number(P)
    if z<=d
        return true
    else
        return false
    end
end


"""
Return whether 'P' is a stable particle or not; that is, if it has more than d contacts.

See also: [`is_a_rattler`](@ref).
"""
is_stable_particle(P::AbstractParticle{d, T}) where {d, T<:AbstractFloat} = !is_a_rattler(P)
is_stable_particle(Particle(rand(3), rand(), 0.4, rand(3,3), rand(3), collect(1:3))); 


"Compute the sum of forces acting on a given particle. The output is a `StaticVector`"
function total_force(P::AbstractParticle{d, T}) where {d, T<:AbstractFloat}
    fs = P.forces
    ctc_vecs = P.contact_vecs
    # tot_f = @SVector zeros(size(P.X, 1))
    tot_f = @SVector zeros(d)
    for (k, f) in enumerate(fs)
        tot_f += f*ctc_vecs[k]
    end
    return tot_f
end



"""
Test whether the sum of forces in a set of particles is smaller than a given tolerance value.

By default such tolerance is `default_tol_force_equilibrium` (a `const`) and equals 1e-12.
"""
function force_equilibrium(particles::Vector{<:AbstractParticle} ; tol_mechanical_equilibrium::Float64=default_tol_force_equilibrium)
    all(x-> x <= tol_mechanical_equilibrium, norm.(total_force.(particles)))
end
