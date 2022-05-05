include("Particles.jl")
import Base: length, size

"Abstract struct for all the type of packings to be used."
abstract type AbstractPacking{d, T} end

#########################################################################################################
#########################################################################################################
### Functions for generic (i.e. any type of packings) packings 
#########################################################################################################
#########################################################################################################

get_positions(packing::AbstractPacking) = getfield.(packing.Particles, :X)
length(packing::AbstractPacking) = length(packing.Particles)
size(packing::AbstractPacking{d, T}) where {d, T<:AbstractFloat} = (d, length(packing))

"Obtain the set of stable particles in a packing."
function get_non_rattlers(packing::AbstractPacking)
    N = length(packing)
    return (1:N)[is_stable_particle.(packing.Particles)]
end

"Obtain the set of rattlers (i.e. unstable particles) in a packing."
function get_rattlers(packing::AbstractPacking)
    N = length(packing)
    return (1:N)[is_a_rattler.(packing.Particles)]
end

function get_coordination_number(packing::AbstractPacking; only_stable::Bool=false)
    if only_stable # so only the contacts of non-rattlers are considered
        non_rattlers = get_non_rattlers(packing)
        return get_coordination_number.(packing.Particles[non_rattlers])
    else 
        return get_coordination_number.(packing.Particles)
    end
end

"Test whether a given packing is isostatic or not. The output is a boolean"
function is_isostatic(packing::AbstractPacking{d, T}; only_stable::Bool=true)::Tuple{Bool, Int64, Int64} where {d, T<:AbstractFloat}
    Nnr = length(get_non_rattlers(packing)); # amount of non-rattlers (i.e. stable particles)
    N_dof = d*(Nnr-1)+1 # number of degrees of freedom
    Nc = 0.5*sum(get_coordination_number(packing; only_stable=only_stable)) # number of contacts; if 'only_stable=true' only non-rattlers are considered
    return Nc==N_dof, Nc, Nnr
end

"Output an array of `StaticVector`s corresponding to the sum of forces acting on each particle of the packing."
function total_force(packing::AbstractPacking{d, T})::Vector{SVector{d, T}} where {d, T<:AbstractFloat}
    total_force.(packing.Particles)
end


distances_between_centers(packing::AbstractPacking) = distances_between_centers(get_positions(packing))

function difference_in_packings(pack1::AbstractPacking, pack2::AbstractPacking)
    Xs1 = get_positions(pack1); Xs2 = get_positions(pack2)
    isostatic, Nc, Nnr = is_isostatic(pack1)
    if isostatic
        non_rattlers = get_non_rattlers(pack1)
        norms = norm.(Xs1-Xs2)[non_rattlers]
    else
        norms = norm.(Xs1-Xs2)
    end
    return mean(norms), std(norms)
end

#########################################################################################################
#########################################################################################################
### Monodisperse packing type and related constructors
#########################################################################################################
#########################################################################################################


#When  dealing with monodisperse packings it's more efficient to store the particles' radius as a field of
#the packing itself rather than of each particle. Of course, in case of polydisperse configurations, that would change.
"""
    MonoPacking{d,T}

Structure used to store all the relevant info of a *monodisperse* hard-spheres packing in d-dimensions.

Its fields are:
1. 'Particles': An vector of `MonoParticle` that contains all the particles that make the packing.
2. 'R': The radius of *all* particles.
3. 'mechanical_equilibrium': A boolean that specifies whether force balances is satisfied
all across the packing (i.e. for each particle).
4. 'jammed': A boolean that specifies whether a packing is jammed or not.

Note that a packing might not be necessarily jammed nor in mechanical_equilibrium. That is 
why specifying such fields is important.

See also: [`MonoParticle`](@ref), [`PolyPacking`](@ref).
"""
mutable struct MonoPacking{d,T} <: AbstractPacking{d, T}
    Particles::Vector{MonoParticle{d, T}} # Vector of particles, each with info of position, contacts, forces, and indices of contacts
    R::T  # The radius of all the particles in a packing
    mechanical_equilibrium::Bool  # boolean to indicate if the packing is in mechanical equilibrium (it should be, in general)
    jammed::Bool # Boolean to indicate if the packing corresponds to a jammed state or not
end
MonoPacking([MonoParticle(rand(3), rand(), rand(3,3), rand(3), collect(1:3)), MonoParticle(rand(3), rand(), rand(3,3), rand(3), collect(1:3))], rand(), false, false)

#### Format how a `MonoPacking` is shown in IO
function Base.show(io::IO, packing::MonoPacking{d,T}) where {d, T<:AbstractFloat}
    N = length(packing)
    iso, Nc, Nnr = is_isostatic(packing)
    fr = (N-Nnr)/Nnr
    println(io, d, "d Monodisperse packing\t of N= ", N, " particles \t of radius R= ", packing.R)
    println(io, Nnr, " stable particles; \t fraction of rattlers = ", round(fr, digits=3))
    println(io, "jammed: ", packing.jammed, "\t isostatic: ", iso ,"\t and in mechanical equilibrium: ", packing.mechanical_equilibrium)
end


#####################
# Now we define some useful constructors
#####################

# Empty packing
"Create an empty d-dimensional packing of N particles"
MonoPacking(d::Int64, N::Int64=1) = MonoPacking(Vector{MonoParticle{d, Float64}}(undef, N), 0.0, false, false)
MonoPacking(5)

@doc raw"""
    MonoPacking(Xs::Vector{SVector{d, PeriodicNumber{T}}}, contact_vecs::Vector{Vector{SVector{d,T}}}, fs::Vector{Vector{T}}, neighbours::Vector{Vector{Int64}}, R::T, jammed::Bool=false; <keyword arguments>)
    MonoPacking(Xs::Vector{SVector{d, PeriodicNumber{T}}}, contact_vecs::Vector{Matrix{T}}, fs::Vector{Vector{T}}, neighbours::Vector{Vector{Int64}}, R::T, jammed::Bool=false; <keyword arguments>)

Create a d-dimensional packing, where the attributes of each particles are inferred from the
 elements of the input arrays.

The number of particles (N) is inferred from the length of the input arrays. 
Then a Vector{MonoParticle{d,T}} of size N but undefined elements is constructed. Each of its
 elements is then defined by calling 'MonoParticle(Xs[i], contact_vecs[i], fs[i], neighbours[i])'.

The constructors also asses whether force balance for each particle is satisfied, within a 
given precision 'tol_mechanical_equilibrium' (that defaults to `default_tol_force_equilibrium=1e-12`).
When this condition is not met, it throws a warning, but the packing is created.
"""
MonoPacking()

function MonoPacking(Xs::Vector{SVector{d, PeriodicNumber{T}}}, contact_vecs::Vector{Vector{SVector{d,T}}}, fs::Vector{Vector{T}}, neighbours::Vector{Vector{Int64}}, R::T, jammed::Bool=false; tol_mechanical_equilibrium::T=default_tol_force_equilibrium, verbose::Bool=true) where {d, T<:AbstractFloat}
    N = length(Xs) # number of particles
    # th construct the packing all the input arrays should have the same length
    if N==length(fs) && N==length(contact_vecs) && N==length(neighbours) 
        particles = Vector{MonoParticle{d,T}}(undef, N) # empty array of particles
        mech_eq = true # default value of mechanical equilibrium condition

        # create particles with the elements of the input arrays
        for i in 1:N
            particles[i] = MonoParticle(Xs[i], contact_vecs[i], fs[i], neighbours[i])
        end

    # test whether the mechanical equilibrium condition is satisfied for all particles, at least within a 'tol_mechanical_equilibrium' precision
        any(norm.(total_force.(particles)).>tol_mechanical_equilibrium) && (mech_eq=false)

        # if this is not the case throw a warning, but create the packing anyway
        if !mech_eq && verbose
            fmismatch = maximum(norm.(total_force.(particles)))
            @warn "Force balance condition is NOT satisfied! Max force mismatch = $fmismatch ;\tCreating the packing anyway."
        end
        
        return MonoPacking(particles, R, mech_eq, jammed) # return the packing created

    else # if there is a mismatch in the dimensions of the input array throw an error (and some info)
        error("Number of dimensions mismatch! There are: ", N, " particles \t ", length(contact_vecs) ," contact vectors \t", length(fs), " forces and \t ", length(neighbours), " list of contacts indices")
    end
end

function MonoPacking(Xs::Vector{SVector{d, PeriodicNumber{T}}}, contact_vecs::Vector{Matrix{T}}, fs::Vector{Vector{T}}, neighbours::Vector{Vector{Int64}}, R::T, jammed::Bool=false; tol_mechanical_equilibrium::T=default_tol_force_equilibrium, verbose::Bool=true) where {d, T<:AbstractFloat}
    MonoPacking(Xs, svectors.(contact_vecs, Val(d)), fs, neighbours, R, jammed; tol_mechanical_equilibrium=tol_mechanical_equilibrium, verbose=verbose)
end


function MonoPacking(particles::Vector{MonoParticle{d, T}}, R::T, jammed::Bool=false;  tol_mechanical_equilibrium::T=default_tol_force_equilibrium, verbose::Bool=true) where {d, T<:AbstractFloat}
    force_balance = force_equilibrium(particles; tol_mechanical_equilibrium=tol_mechanical_equilibrium)
    if !force_balance && verbose
        fmismatch = maximum(norm.(total_force.(particles)))
        @warn "Force balance condition is NOT satisfied! Max force mismatch = $fmismatch ; \tCreating the packing anyway."
    end

    MonoPacking(particles, R, force_balance, jammed)
end



"""
    distances_between_particles(packing::MonoPacking)

Compute the distance between the *borders* of all pairs of particles in 'packing'.

For a given pair, the distance between the centers of the two particles is computed, by calling 
`distances_between_centers`, and then the diameter is subtracted. In this way, the result is 
AbstractFloatly the distance between the border of each pair of particles.

See also [`distances_between_centers`](@ref)
"""
function distances_between_particles(packing::MonoPacking)
    distances_between_centers(packing) .- (2*packing.R)
end

"""
    check_for_overlaps(packing::MonoPacking, tolerance::AbstractFloat)

Apply `check_for_overlaps` to all the particles in 'packing'.
"""
function check_for_overlaps(packing::MonoPacking{d, T}, tolerance::T) where {d, T<:AbstractFloat}
    R = packing.R
    Xs = get_positions(packing)
    check_for_overlaps(Xs, R, tolerance)
end

#################################################
# Other useful functions
################################################

"""
    packing_fraction(d::Int64, R::AbstractFloat, N::Int64, L=1.0)

Compute the packing fraction of N d-dimensional hyperspheres of the same radius, 'R', inside
 a box of size 'L'.

See also [`volume`](@ref), [`volume_d_ball`](@ref)
 """
function packing_fraction(d::Int64, R::T, N::Int64, L::T=1.0) where T<:AbstractFloat
    return N*volume_d_ball(d, R)/L^d
end
packing_fraction(3, 0.5, 2, 1.0)

"""
    packing_fraction(packing::MonoPacking{d, T})

Compute the packing fraction of a monodisperse packing.

See also [`volume`](@ref), [`volume_d_ball`](@ref)
"""
function packing_fraction(packing::MonoPacking{d, T}) where {d, T<:AbstractFloat}
    R = packing.R
    N = length(packing.Particles)
    L = packing.Particles[1].X[1].L
    packing_fraction(d, R, N, L)
end


#########################################################################################################
#########################################################################################################
### Polydisperse packing type and related constructors
#########################################################################################################
#########################################################################################################
"""
    PolyPacking{d,T}

Structure used to store all the relevant info of a *monodisperse* hard-spheres packing in d-dimensions.

Its fields are:
1. 'Particles': An vector of `Particle`s that contains all the particles that make the packing.
2. 'mechanical_equilibrium': A boolean that specifies whether force balances is satisfied
all across the packing (i.e. for each particle).
3. 'jammed': A boolean that specifies whether a packing is jammed or not.

Note that a packing might not be necessarily jammed nor in mechanical_equilibrium. That is 
why specifying such fields is important.

See also: [`Particle`](@ref), [`MonoPacking`(@ref)].
"""
mutable struct PolyPacking{d,T} <: AbstractPacking{d, T}
    Particles::Vector{Particle{d, T}} # Vector of particles, each with info of position, contacts, forces, and indices of contacts
    mechanical_equilibrium::Bool  # boolean to indicate if the packing is in mechanical equilibrium (it should be, in general)
    jammed::Bool # Boolean to indicate if the packing corresponds to a jammed state or not
end
PolyPacking([Particle(rand(3), 0.1, rand(), rand(3,3), rand(3), collect(1:3)), Particle(rand(3), 0.1, rand(), rand(3,3), rand(3), collect(1:3))],  false, false)


function get_radii(packing::PolyPacking{d,T}) where {d, T<:AbstractFloat}
    getfield.(packing.Particles, :R)
end
get_radii(PolyPacking([Particle(rand(3), 0.1, rand(), rand(3,3), rand(3), collect(1:3)), Particle(rand(3), 0.5, rand(), rand(3,3), rand(3), collect(1:3))],  false, false))

#### Format how a `PolyPacking` is shown in IO
function Base.show(io::IO, packing::PolyPacking{d,T}) where {d, T<:AbstractFloat}
    N = length(packing)
    iso, Nc, Nnr = is_isostatic(packing)
    fr = (N-Nnr)/Nnr
    println(io, d, "d Polydisperse packing\t of N= ", N, " particles \t of radii = ", get_radii(packing))
    println(io, Nnr, " stable particles; \t fraction of rattlers = ", round(fr, digits=3))
    println(io, "jammed: ", packing.jammed, "\t isostatic: ", iso ,"\t and in mechanical equilibrium: ", packing.mechanical_equilibrium)
end


#####################
# Now we define some useful constructors
#####################

# Empty packing
"Create an empty d-dimensional packing of N particles"
PolyPacking(d::Int64, N::Int64=1) = PolyPacking(Vector{Particle{d, Float64}}(undef, N), false, false)
PolyPacking(5)

@doc raw"""
    PolyPacking(Xs::Vector{SVector{d, PeriodicNumber{T}}}, contact_vecs::Vector{Vector{SVector{d,T}}}, fs::Vector{Vector{T}}, neighbours::Vector{Vector{Int64}}, jammed::Bool=false; <keyword arguments>)
    PolyPacking(Xs::Vector{SVector{d, PeriodicNumber{T}}}, contact_vecs::Vector{Matrix{T}}, fs::Vector{Vector{T}}, neighbours::Vector{Vector{Int64}}, jammed::Bool=false; <keyword arguments>)

Create a d-dimensional polydisperse packing, where the attributes of each particles are inferred from the
 elements of the input arrays.

The number of particles (N) is inferred from the length of the input arrays. 
Then a Vector{Particle{d,T}} of size N but undefined elements is constructed. Each of its
 elements is then defined by calling 'Particle(Xs[i], contact_vecs[i], fs[i], neighbours[i])'.

The constructors also asses whether force balance for each particle is satisfied, within a 
given precision 'tol_mechanical_equilibrium' (that defaults to `default_tol_force_equilibrium=1e-12`).
When this condition is not met, it throws a warning, but the packing is created.
"""
PolyPacking()

function PolyPacking(Xs::Vector{SVector{d, PeriodicNumber{T}}}, contact_vecs::Vector{Vector{SVector{d,T}}}, fs::Vector{Vector{T}}, neighbours::Vector{Vector{Int64}}, jammed::Bool=false; tol_mechanical_equilibrium::T=default_tol_force_equilibrium, verbose::Bool=true) where {d, T<:AbstractFloat}
    N = length(Xs) # number of particles
    # th construct the packing all the input arrays should have the same length
    if N==length(fs) && N==length(contact_vecs) && N==length(neighbours) 
        particles = Vector{Particle{d,T}}(undef, N) # empty array of particles
        mech_eq = true # default value of mechanical equilibrium condition

        # create particles with the elements of the input arrays
        for i in 1:N
            particles[i] = PolyParticle(Xs[i], contact_vecs[i], fs[i], neighbours[i])
        end

    # test whether the mechanical equilibrium condition is satisfied for all particles, at least within a 'tol_mechanical_equilibrium' precision
        any(norm.(total_force.(particles)).>tol_mechanical_equilibrium) && (mech_eq=false)

        # if this is not the case throw a warning, but create the packing anyway
        if !mech_eq && verbose
            fmismatch = maximum(norm.(total_force.(particles)))
            @warn "Force balance condition is NOT satisfied! Max force mismatch = $fmismatch ;\tCreating the packing anyway."
        end
        
        return PolyPacking(particles, mech_eq, jammed) # return the packing created

    else # if there is a mismatch in the dimensions of the input array throw an error (and some info)
        error("Number of dimensions mismatch! There are: ", N, " particles \t ", length(contact_vecs) ," contact vectors \t", length(fs), " forces and \t ", length(neighbours), " list of contacts indices")
    end
end

function PolyPacking(Xs::Vector{SVector{d, PeriodicNumber{T}}}, contact_vecs::Vector{Matrix{T}}, fs::Vector{Vector{T}}, neighbours::Vector{Vector{Int64}}, jammed::Bool=false; tol_mechanical_equilibrium::T=default_tol_force_equilibrium, verbose::Bool=true) where {d, T<:AbstractFloat}
    PolyPacking(Xs, svectors.(contact_vecs, Val(d)), fs, neighbours, jammed; tol_mechanical_equilibrium=tol_mechanical_equilibrium, verbose=verbose)
end


function PolyPacking(particles::Vector{Particle{d, T}}, jammed::Bool=false;  tol_mechanical_equilibrium::T=default_tol_force_equilibrium, verbose::Bool=true) where {d, T<:AbstractFloat}
    force_balance = force_equilibrium(particles; tol_mechanical_equilibrium=tol_mechanical_equilibrium)
    if !force_balance && verbose
        fmismatch = maximum(norm.(total_force.(particles)))
        @warn "Force balance condition is NOT satisfied! Max force mismatch = $fmismatch ; \tCreating the packing anyway."
    end

    PolyPacking(particles, force_balance, jammed)
end



"""
    distances_between_particles(packing::PolyPacking)

Compute the distance between the *borders* of all pairs of particles in 'packing'.

For a given pair, the distance between the centers of the two particles is computed, by calling 
`distances_between_centers`, and then the sum of their radii is subtracted. In this way, the result is 
really the distance between the border of each pair of particles.

See also [`distances_between_centers`](@ref)
"""
function distances_between_particles(packing::PolyPacking)
    N = length(packing)
    Rs = get_radii(packing)
    dists = distances_between_centers(packing).data
    for i in 1:N, j in i+1:N
        dists[i, j] -= Rs[i]+Rs[j]
    end
    Symmetric(dists)
end

"""
    check_for_overlaps(packing::PolyPacking{d,T}, tolerance::T)

Apply `check_for_overlaps` to all the particles in 'packing'.
"""
function check_for_overlaps(packing::PolyPacking{d, T}, tolerance::T) where {d, T<:AbstractFloat}
    Rs = get_radii(packing)
    Xs = get_positions(packing)
    check_for_overlaps(Xs, Rs, tolerance)
end
# function check_for_overlaps(packing::PolyPacking{d, T}, tolerance::T) where {d, T<:AbstractFloat}
#     parts = packing.Particles
#     N = length(packing)
#     overlap = false; message = "No overlap is present"; ovlp_indices = (0,0)
#     for i in 1:N, j in i+1:N
#         overlap, gap = check_for_overlaps(parts[i], parts[j], tolerance)
#         if overlap
#             Xi = value.(parts[i].X); Ri = parts[i].R
#             Xj = value.(parts[j].X); Rj = parts[j].R

#             message = "Overlap between particles $i and $j;\nCenters at: $Xi; \t $Xj. \n\t\tradii:  ($R1, $R2); \t Overlap = $gap ."
#             ovlp_indices = (i,j)
#             break
#         end
#     end
#     return overlap, message, ovlp_indices
# end


#################################################
# Other useful functions
################################################
"""
    packing_fraction(d::Int64, R1::T, N1::Int64, R2::T, N2::Int64, L::T=1.0) where T<:AbstractFloat

Compute the packing fraction of a configuration of d-dimensional hyperspheres, of which 'N1' have radius 'R1', and 'N2' have radius 'R2'. The configuration is assumed to be inside  a box of size 'L'.

See also [`volume`](@ref), [`volume_d_ball`](@ref)
 """
function packing_fraction(d::Int64, R1::T, N1::Int64, R2::T, N2::Int64, L::T=1.0) where T<:AbstractFloat
    return (N1*volume_d_ball(d, R1)+ N2*volume_d_ball(d, R2))/L^d
end
packing_fraction(3, 0.1, 5, 0.3, 8, 1.0)


"""
    packing_fraction(d::Int64, Rs::Vector{T}, L::T=1.0) where T<:AbstractFloat

Compute the packing fraction of N d-dimensional hyperspheres with radii 'Rs', inside
 a box of size 'L'.

See also [`volume`](@ref), [`volume_d_ball`](@ref)
 """
function packing_fraction(d::Int64, Rs::Vector{T}, L::T=1.0) where T<:AbstractFloat
    return sum(volume_d_ball.(d, Rs))/L^d
end
packing_fraction(3, rand(5), 1.0)

"""
    packing_fraction(packing::PolyPacking{d, T})

Compute the packing fraction of a monodisperse packing.

See also [`volume`](@ref), [`volume_d_ball`](@ref)
"""
function packing_fraction(packing::PolyPacking{d, T}) where {d, T<:AbstractFloat}
    Rs = get_radii(packing)
    L = packing.Particles[1].X[1].L
    packing_fraction(d, Rs, L)
end