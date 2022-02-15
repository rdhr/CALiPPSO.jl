include("Particles.jl")
import Base: length, size

"Abstract struct for all the type of packings to be used."
abstract type AbstractPacking{d, T} end

#########################################################################################################
#########################################################################################################
### Functions for generic (i.e. any type of packings) packings 
#########################################################################################################
#########################################################################################################
"Obtain the set of stable particles in a packing."
function get_non_rattlers(packing::AbstractPacking)
    N = length(packing.Particles)
    return (1:N)[is_stable_particle.(packing.Particles)]
end

"Obtain the set of rattlers (i.e. unstable particles) in a packing."
function get_rattlers(packing::AbstractPacking)
    N = length(packing.Particles)
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
function is_isostatic(packing::AbstractPacking{d, T}; only_stable::Bool=true)::Tuple{Bool, Int64, Int64} where {d, T<:Real}
    Nnr = length(get_non_rattlers(packing)); # amount of non-rattlers 
    N_dof = d*(Nnr-1)+1 # number of degrees of freedom
    Nc = 0.5*sum(get_coordination_number(packing; only_stable=only_stable)) # number of contacts; if 'only_stable=true' only non-rattlers are considered
    return Nc==N_dof, Nc, Nnr
end

"Output an array of `StaticVector`s corresponding to the sum of forces acting on each particle of the packing."
function total_force(packing::AbstractPacking{d, T})::Vector{SVector{d, T}} where {d, T<:Real}
    total_force.(packing.Particles)
end


get_positions(packing::AbstractPacking) = getfield.(packing.Particles, :X)
length(packing::AbstractPacking) = length(get_positions(packing))
size(packing::AbstractPacking{d, T}) where {d, T<:Real} = (d, length(packing))

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

See also: [`MonoParticle`](@ref).
"""
mutable struct MonoPacking{d,T} <: AbstractPacking{d, T}
    Particles::Vector{MonoParticle{d, T}} # Vector of particles, each with info of position, contacts, forces, and indices of contacts
    R::T  # The radius of all the particles in a packing
    mechanical_equilibrium::Bool  # boolean to indicate if the packing is in mechanical equilibrium (it should be, in general)
    jammed::Bool # Boolean to indicate if the packing corresponds to a jammed state or not
end
MonoPacking([MonoParticle(rand(3), rand(), rand(3,3), rand(3), collect(1:3)), MonoParticle(rand(3), rand(), rand(3,3), rand(3), collect(1:3))], rand(), false, false)

#### Format how a `MonoPacking` is shown in IO
function Base.show(io::IO, packing::MonoPacking{d,T}) where {d, T<:Real}
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

function MonoPacking(Xs::Vector{SVector{d, PeriodicNumber{T}}}, contact_vecs::Vector{Vector{SVector{d,T}}}, fs::Vector{Vector{T}}, neighbours::Vector{Vector{Int64}}, R::T, jammed::Bool=false; tol_mechanical_equilibrium::Float64=default_tol_force_equilibrium, verbose::Bool=true) where {d, T<:Real}
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

function MonoPacking(Xs::Vector{SVector{d, PeriodicNumber{T}}}, contact_vecs::Vector{Matrix{T}}, fs::Vector{Vector{T}}, neighbours::Vector{Vector{Int64}}, R::T, jammed::Bool=false; tol_mechanical_equilibrium::Float64=default_tol_force_equilibrium, verbose::Bool=true) where {d, T<:Real}
    MonoPacking(Xs, svectors.(contact_vecs, Val(d)), fs, neighbours, R, jammed; tol_mechanical_equilibrium=tol_mechanical_equilibrium, verbose=verbose)
end


# Some quantities for calling previously defined functions for the first time and compile them.
# N_comp=10; d_comp=2;
# Xs_comp = PeriodicVectors(rand(d_comp, N_comp))
# zs_comp = 2*rand(1:4, N_comp)

# cvec_comp = [rand(d_comp, zs_comp[i]) for i in 1:N_comp]
# st_cvec_comp = svectors.(cvec_comp, Val(d_comp))

# fs_comp = [rand(zs_comp[i]) for i in 1:N_comp]
# nghs_lists_comp = [rand(1:N_comp, zs_comp[i]) for i in 1:N_comp]

# packing_comp = MonoPacking(Xs_comp, st_cvec_comp, fs_comp, nghs_lists_comp, rand(), false; verbose=false)
# MonoPacking(Xs_comp, cvec_comp, fs_comp, nghs_lists_comp, rand(), false; verbose=false)
# force_equilibrium(packing_comp.Particles)
# println(packing_comp)


# get_coordination_number(packing_comp)
# get_non_rattlers(packing_comp)
# get_rattlers(packing_comp)
# is_isostatic(packing_comp)
# total_force(packing_comp)
# get_positions(packing_comp) == Xs_comp
# length(packing_comp); size(packing_comp)
# difference_in_packings(packing_comp, packing_comp)

# Finish calling previously defined functions


function MonoPacking(particles::Vector{MonoParticle{d, T}}, R::T, jammed::Bool=false;  tol_mechanical_equilibrium::Float64=default_tol_force_equilibrium, verbose::Bool=true) where {d, T<:Real}
    force_balance = force_equilibrium(particles; tol_mechanical_equilibrium=tol_mechanical_equilibrium)
    if !force_balance && verbose
        fmismatch = maximum(norm.(total_force.(particles)))
        @warn "Force balance condition is NOT satisfied! Max force mismatch = $fmismatch ; \tCreating the packing anyway."
    end

    MonoPacking(particles, R, force_balance, jammed)
end
# MonoPacking(packing_comp.Particles, rand(); verbose=false)


distances_between_centers(packing::AbstractPacking) = distances_between_centers(get_positions(packing))
# distances_between_centers(packing_comp)


"""
    distances_between_particles(packing::MonoPacking)

Compute the distance between the *borders* of all pairs of particles in 'packing'.

For a given pair, the distance between the centers of the two particles is computed, by calling 
`distances_between_centers`, and then the diameter is subtracted. In this way, the result is 
really the distance between the border of each pair of particles.

See also [`distances_between_centers`](@ref)
"""
function distances_between_particles(packing::MonoPacking)
    distances_between_centers(packing) .- (2*packing.R)
end
# distances_between_particles(packing_comp)

"""
    check_for_overlaps(packing::MonoPacking, tolerance::Float64)

Apply `check_for_overlaps` to all the particles in 'packing'.
"""
function check_for_overlaps(packing::MonoPacking, tolerance::Float64)
    R = packing.R
    Xs = get_positions(packing)
    check_for_overlaps(Xs, R, tolerance)
end
# check_comp = check_for_overlaps(packing_comp, 1e-8)
#################################################
# Other useful functions
################################################

"""
    packing_fraction(d::Int64, R::Real, N::Int64, L::Real=1.0)

Compute the packing fraction of N d-dimensional hyperspheres of the same radius, 'R', inside
 a box of size 'L'.

See also [`volume`](@ref), [`volume_d_ball`](@ref)
 """
function packing_fraction(d::Int64, R::Real, N::Int64, L::Real=1.0)::Float64
    return N*volume_d_ball(d, R)/L^d
end
packing_fraction(3, 0.5, 2, 1.0)

"""
    packing_fraction(packing::MonoPacking{d, T})

Compute the packing fraction of a monodisperse packing.

See also [`volume`](@ref), [`volume_d_ball`](@ref)
"""
function packing_fraction(packing::MonoPacking{d, T})::Float64 where {d, T}
    R = packing.R
    N = length(packing.Particles)
    L = packing.Particles[1].X[1].L
    packing_fraction(d, R, N, L)
end
# packing_fraction(packing_comp)

