module CALiPPSO

export produce_jammed_configuration!, precompile_main_function # main function!! and a function used to precompile it on a small system
export default_parameters # Dict of default parameters defined in this module
export generate_random_configuration, radius_from_phi # function to generate a low density, random initial configuration; and compute r from φ
export convergence_info, PeriodicNumber, MonoParticle, Particle, MonoPacking, PolyPacking, PeriodicVector, Packing # `struct`s defined in the package
export network_of_contacts, check_for_overlaps, PeriodicVectors, packing_fraction, get_non_rattlers, get_rattlers, is_isostatic, get_coordination_number, total_force, difference_in_packings # other useful functions; mostly for packings
export volume_d_ball, norm # these functions are needed for generating random initial packings, in the 'random_initial_conditions.jl' script. Nevertheless, it could be the case that they're also useful when analysing results, specially 'norm'

include("Packing.jl")
include("random_initial_conditions.jl")
include("functions-printing-output.jl")

using JuMP # Library for using high level interface for optimization and model creation/manipulation 
import JuMP.MOI as MOI
const max_threads = Int(round(Sys.CPU_THREADS/2)) # max number of processes to be used by any optimizer. By default, half the number of threads available in a computer

##############################################################################
##############################################################################
## Some possibly useful parameters with the different solvers
##############################################################################
##############################################################################

################################################
## We suggest to use Gurobi if you can get a license
################################################

# ### Parameters for Gurobi
# using Gurobi
# const default_tol_overlap=1e-8 # tolerance for identifying an Overlap. Given that Gurobi's precision is 10^-9, a larger value is needed.
# const default_tol_optimality = 1e-9; # optimality tolerances. This is the most precise value allowed by Gurobi
# const grb_env = Gurobi.Env() # to avoid printing license info every time a model is created or optimized 
# const default_optimizer = Gurobi.Optimizer(grb_env)
# #= 
# We are omitting all the standard output of Gurobi by setting the 'OutputFlag' to 0.
# The values of feasibility (FeasibilityTol) and optimality (OptimalityTol) chosen are the most stringent ones allowed by Gurobi.
# Instead, the choice of using the 3rd method (i.e. the barrier method) might lead, in some cases, to a slower performance in comparison with the default method '0', which is a concurrent implementation of simplex (primal and dual) as well as barrier solvers. However, for reproducibility and better characterization, restricting to only using barrier method is better. Of course, this can be changed without any problems.
# More info can be found here: https://www.gurobi.com/resource/parallelism-linear-mixed-integer-programming/ and in the corresponding section of Gurobi's manual: https://www.gurobi.com/documentation/9.1/refman/method.html
# =#
# const default_solver_attributes = Dict("OutputFlag" => 0, "FeasibilityTol" => default_tol_optimality, "OptimalityTol" => default_tol_optimality, "Method" => 3, "Threads" => max_threads)



################################################
## Suggested parameters if you want to use HiGHS solver
################################################

# # #### Parameters for HiGHS
# using HiGHS
# const default_tol_overlap=1e-8 # tolerance for identifying an Overlap.
# const default_tol_optimality = 1e-9
# const default_optimizer = HiGHS.Optimizer()
# const default_solver_attributes = Dict("small_matrix_value"=>0.1*default_tol_optimality, "primal_feasibility_tolerance" => default_tol_optimality, "dual_feasibility_tolerance" => default_tol_optimality, "solver" => "ipm", "ipm_optimality_tolerance"=>default_tol_optimality,  "highs_max_threads" => max_threads, "parallel"=>"on", "output_flag"=>false)



###############################################
# Suggested parameters if you want to use Clp solver
###############################################

# #### Parameters for Clp
# using Clp
# const default_tol_overlap=1e-8 # tolerance for identifying an Overlap.
# const default_tol_optimality = 1e-9
# const default_optimizer = Clp.Optimizer()
# const default_solver_attributes = Dict("PrimalTolerance"=>default_tol_optimality, "DualTolerance" => default_tol_optimality, "LogLevel" => 0, "SolveType" => 5)
# Model(()-> default_optimizer)


################################################
## By default CALiPPSO uses the GLPK solver with the following parameters
################################################

### Parameters for GLPK
using GLPK
const default_tol_overlap=1e-8 # tolerance for identifying an Overlap.
const default_tol_optimality = 0.1*default_tol_overlap
const default_optimizer = GLPK.Optimizer(;want_infeasibility_certificates=false, method=GLPK.MethodEnum(0))
const default_solver_attributes = Dict("msg_lev"=>GLPK.GLP_MSG_OFF, "tol_bnd"=>default_tol_optimality, "tol_dj"=>default_tol_optimality)

### First calls for compilation of model creation
empty!(Model(()-> default_optimizer))



################################################
## Suggested parameters if you want to use COSMO solver
# (**NOT** recommended to use these due to poor precision in identifying active dual variables)
################################################

# ### Parameters for COSMO
# using COSMO
# const default_tol_overlap=1e-5 # tolerance for identifying an Overlap.
# const default_tol_optimality = 0.1*default_tol_overlap
# const default_solver_attributes = Dict("verbose"=>false, "check_termination"=>500, "max_iter"=>10000, "eps_abs"=>default_tol_optimality, "eps_rel"=>default_tol_optimality, "adaptive_rho"=>false)


################################################
## Suggested parameters if you want to use Hypatia solver 
## (**NOT** recommended to use these due to poor precision in identifying active dual variables; although it passes the d=3 test)
################################################


# ### Parameters for Hypatia
# using Hypatia
# const default_tol_overlap=1e-8 # tolerance for identifying an Overlap.
# const default_tol_optimality = 0.1*default_tol_overlap
# const default_optimizer = Hypatia.Optimizer(;verbose = false, tol_abs_opt = default_tol_optimality, tol_feas = 0.1*default_tol_overlap, tol_infeas = 0.1*default_tol_overlap)
# const default_solver_attributes = Dict()



const dict_solver = Dict("default_optimizer"=> default_optimizer,  "default_solver_attributes"=> default_solver_attributes)
##############################################################################
##############################################################################
## Here we begin defining all the required functions and remaining constants
##############################################################################
##############################################################################
const default_max_iterations = 1000
const default_tol_Γ_convergence = 1e-12 # default value for determining the convergence of √Γ-1. NB: this value is only recommended if you're using Gurobi, HiGHS, or GLPK; for different solvers this might be too low so you'll have to change it accordingly
const default_tol_displacements_convergence = 1e-9 # default value for determining the convergence of |sᵢ|
const default_tol_zero_forces = 0.1*default_tol_optimality # any dual variable smaller than this value is considered a 0. Note that this is actually a loose threshold
# const default_tol_zero_forces = 1e-5 # any dual variable smaller than this value is considered a 0. Value recommended for COSMO, Tulip, and Hypatia
# const default_tol_zero_forces = eps() # NB: when using Gurobi, in *many* cases this can actually be set as 0.0 without affecting the accuracy with which the *real* active dual variables are found

const dict_precision = Dict("default_tol_overlap"=>default_tol_overlap, "default_tol_optimality"=>default_tol_optimality, "default_tol_zero_forces"=>default_tol_zero_forces, "default_tol_force_equilibrium"=>default_tol_force_equilibrium)
const dict_convergence = Dict("default_tol_Γ_convergence"=>default_tol_Γ_convergence, "default_tol_displacements_convergence"=>default_tol_displacements_convergence, "default_max_iterations"=>default_max_iterations)

const dict_other = Dict("default_threads"=>max_threads)

const default_parameters = Dict("Convergence parameters"=>dict_convergence, "Precision parameters"=>dict_precision, "Solver parameters"=>dict_solver, "Other parameters"=>dict_other)


#########################################################################################################
#########################################################################################################
### Some other needed (but secondary) struct's and functions
#########################################################################################################
#########################################################################################################

# struct for saving information about the convergence time, resources, etc. of the CALiPPSO algorithm
"
Struct to save convergence information of an CALiPPSO solution. For more info, see docs of its fields: `converged`, `iterations`, `time`, `memory`, `times_LP_optim`."
struct convergence_info
    "Whether or not the convergence criteria were met."
    converged::Bool # whether or not the CALiPPSO algorithm converged (possibly to a jammed packing)
    "Number of iterations to reach convergence (or stopping if exceeded `max_iters`"
    iterations::Int64  # number of iterations needed for convergence
    "Total execution time of the CALiPPSO algorithm."
    time::Float64 # total time of the CALiPPSO algorithm (since the while loop begins), as measured by 'time' (in seconds)
    "Memory allocated while running CALiPPSO."
    memory::Float64 # allocated memory, as measured by 'time' (in GB)
    "Times needed to optimize each LP instance."
    times_LP_optim::Vector{Float64} # list of solving times of each LP optimization step
end
convergence_info(false, 0, 0.0, 0.0, zeros(2))

"Test whether there is an overlap between particles. Also info about displacements 
(in the possible overlapping pair) is given."
function check_for_overlaps(Xs::Vector{SVector{d, PeriodicNumber{T}}}, R::T, S::Matrix{T}, t::Int64, possible_neighbours::Vector{Vector{Int64}}; tol_overlap::T=default_tol_overlap) where {d, T<:AbstractFloat}
    overlap, message, particles = check_for_overlaps(Xs, R, tol_overlap)
    if overlap
        i, j = particles
        nghs_i = possible_neighbours[i]; nghs_j = possible_neighbours[j]
        Si = S[:, i]; Sj = S[:, j];

        message *= "\n Info particle $i: S = $Si;\t particles inducing constraints on it $(nghs_i)"
        message *= "\n Info particle $j: S = $Sj;\t particles inducing constraints on it $(nghs_j)"
        error("At iteration $t: "*message)
    end
end

function check_for_overlaps(Xs::Vector{SVector{d, PeriodicNumber{T}}}, Rs::Vector{T}, S::Matrix{T}, t::Int64, possible_neighbours::Vector{Vector{Int64}}; tol_overlap::T=default_tol_overlap) where {d, T<:AbstractFloat}
    overlap, message, particles = check_for_overlaps(Xs, Rs, tol_overlap)
    if overlap
        i, j = particles
        nghs_i = possible_neighbours[i]; nghs_j = possible_neighbours[j]
        Si = S[:, i]; Sj = S[:, j];

        message *= "\n Info particle $i: S = $Si;\t particles inducing constraints on it $(nghs_i)"
        message *= "\n Info particle $j: S = $Sj;\t particles inducing constraints on it $(nghs_j)"
        error("At iteration $t: "*message)
    end
end

@doc raw"""
    check_for_overlaps(packing::AbstractPacking, t::Int64, possible_neighbours::Vector{Vector{Int64}}, jammed::Bool; tolerance=default_tol_overlap)
Check for overlaps in all the particles of a given packing, obtained after CALiPPSO converged.
"""
function check_for_overlaps(packing::AbstractPacking, t::Int64, possible_neighbours::Vector{Vector{Int64}}, jammed::Bool; tolerance=default_tol_overlap)
    printstyled("Checking for overlaps after convergence...\n", color=:yellow, bold=false)

    overlap, message, particles_indices =check_for_overlaps(packing, tolerance)
    if overlap
        i, j = particles_indices
        nghs_i = possible_neighbours[i]; nghs_j = possible_neighbours[j]
        message *= "\n Particles inducing constraints on $i : $(nghs_i)\n Particles inducing constraints on $j: $(nghs_j)"
        
        jammed ? error("At iteration $t: "*message) : error("At iteration $(t+1): "*message)
    else
        printstyled("No overlaps found! :D\n", color=:green, bold=true)
    end
end


#This function is no longer used, but might be useful for analysing the CALiPPSO behaviour
"Compute the norm of each of the displacement vectors."
function norm_displacements(S::Matrix{Float64})::Vector{Float64} 
    d, N = size(S)
    (d>=N) && (@warn "There is likely an error in the dimensions of displacements matrix! Dimensionality: $d \t N = $N")

    return norm.(svectors(S, Val(d)))
end
norm_displacements(rand(4,50))


"Empty an *optimizer* in case it still is associated with a model"
CALiPPSO.empty!(optimizer::MOI.AbstractOptimizer) = MOI.empty!(optimizer)
CALiPPSO.empty!(default_optimizer)

#########################################################################################################
#########################################################################################################
### Functions to generate and solve the jamming LOP instances
#########################################################################################################
#########################################################################################################

"""
    bounds_and_cutoff(sqrΓ::T, R::T, ℓ0::T, d::Int64; thresholds::Vector{T}=[5e-4, 1e-5], sbound::T=0.01)

Obtain displacements' bound and radius of influence (or cutoff) for particles of radius 'R' 
and after an inflation factor with square root 'sqrΓ'. In case of a polydisperse packing
the largest radius is used.

# Output
1. The cutoff for assigning constraints, i.e. the radius of influence, ℓ.
2. The bounds to be imposed on the each component of the displacements.

The value of the displacements' bound and ℓ are chosen depending on how close the input 
value of the growth factor, Γ (or better √Γ), is to 1. Thus, as explained in the paper, 
even when the final φ_J is not known, the displacements' bounds and ℓ are estimated with 
the value of Γ from the previous LP optimization (or an initial suitable guess if needed).
Three different criteria are selected for large, small, and very small values of sqrΓ-1, 
respectively

The other arguments are the current value of 'R', an upper bound for the radius of influence 
'ℓ0', and the dimensionality of the system, 'd'. 
Keywords arguments are used to control when different criteria are triggered 'thresholds', 
and the fraction of 'R' used as displacement bound when Γ is already very close to 1, 
'sbound' (default 0.01).
"""
function bounds_and_cutoff(sqrΓ::T, R::T, ℓ0::T, d::Int64; thresholds::Tuple{T, T}=(5e-4, 1e-5), sbound::T=0.01) where {T<:Real}
    th1, th2 = thresholds

    if sqrΓ > 1+th1 # For large or intermediate values of inflation factor, choose cutoff and bounds loosely
        s_bound = (0.5ℓ0 - max(1.5, sqrΓ)*R)/sqrt(d)
        return ℓ0, s_bound
    elseif th1 >= sqrΓ-1 > th2 # For small values, the cutoff should be only neighbours within two shells
        return min(4R, ℓ0), 0.1R
    else  # For even smaller values, just a fraction outside the first shell of neighbours should be enough
        # also, displacements are bound to a very small fraction of R
        return min(2.7*R, ℓ0), sbound*R
    end
end
bounds_and_cutoff(5.0, 0.52, 3.0, 3)


"""
    add_non_overlapping_constraints!(model::JuMP.Model, Xs::Vector{SVector{d, PeriodicNumber{T}}}, R::T, ℓ::T, images::Vector{SVector{d, T}})
    add_non_overlapping_constraints!(model::JuMP.Model, Xs::Vector{SVector{d, PeriodicNumber{T}}}, Rs::Vector{T}, ℓ::T, images::Vector{SVector{d, T}})

Assign the linearized non-overlapping constraints to 'model', according to particle 
positions, 'Xs', and radius, 'R' (or radii 'Rs' in polydisperse systems).

A pair of particles is included in the set of constraints only if their distance is smaller 
than the cutoff 'ℓ'. Because the MIC distance is considered, also the set of virtual 
images 'images' should also be provided as input.

# Output
1. 'constraints': a Vector whose elements are themselves arrays of 'ConstraintRef' type.
2. `nearby_particles_list`: a vector whose i-th entry is the list of particles' indices that induce a constraint on particle i.

Note that the constraint associated to the pair (i,j) is only counted once (is associated 
to 'i' if j>i; or to 'j' otherwise). Thus the i-th entry of the output arrays only contains
constraints and indices associated to particles of index greater than 'i'.

See also [`MIC_distance`](@ref), [`MIC_vector`](@ref), [`@constraint`](https://jump.dev/JuMP.jl/stable/reference/constraints/#JuMP.@constraint), 
[`solve_LP_instance!`](@ref), [`fine_tune_forces!`](@ref).
"""
function add_non_overlapping_constraints!(model::JuMP.Model, Xs::Vector{SVector{d, PeriodicNumber{T}}}, R::T, ℓ::T, images::Vector{SVector{d, T}}, distances::Symmetric{T, Matrix{T}}) where {d, T<:AbstractFloat}
    #TODO: check if the method of distances_between_centers(centers, images) yields a more efficient function. In such case, the indices of images are also computed from the beginning
    
    N = length(Xs); # number of particles
    # access the design variables of 'model' (S corresponds to the d*N array of the displacements; the inflation factor determines how much the radius of each particle is increased)
    S⃗ = model[:S];  Γ = model[:inflation_factor]

    # Because we're using ordered pairs [i.e. (i,j) with (i<j)] to enumerate the constraints, only N-1 entries are needed in the following arrays
    constraints = Vector{Vector{ConstraintRef}}(undef, N-1) # initialization of vector for storing the lists of constraints
    nearby_particles_list = Vector{Vector{Int64}}(undef, N-1) # initialization of vector for storing the lists of particles inducing constraints

    # We will only consider nearby particles --i.e. particles inside a neighbourhood of size 'ℓ' around-- as possible neighbours
    # And hence, the non-overlapping constraint will only be applied to this smaller set of particles
    for i in 1:N-1::Int
        #NB. We reduce the number of constraints, by considering neighbours with index >i.
        neigh_parts = collect(i+1:N)[distances[i+1:N, i] .<= ℓ]  # indices of particles whose distance to i is smaller than ℓ
        Nn = length(neigh_parts) # number of particles nearby i
        nearby_particles_list[i] = neigh_parts; # store the list of nearby particles
        
        constrs_part_i = Vector{ConstraintRef}(undef, Nn) # temporary array to store the constraints associated to particle i
        # this loop adds the relevant constraints of particle i to 'model' AND stores them in 'constrs_part_i'
        for (num, j) in enumerate(neigh_parts)
            # to determine the vector joining Xi with a given neighbour we used the MIC distance
            Xij, mic_image = MIC_vector(Xs[i], Xs[j], images)
            constrs_part_i[num] = @constraint(model, 2*dot(Xij, S⃗[:, j]-S⃗[:, i]) + Γ*(2R)^2 <= norm(Xij)^2 )
        end
        constraints[i] = constrs_part_i #store the list of constraints of particle i as the i-th entry of 'constraints'

    end
    return constraints, nearby_particles_list
end

function add_non_overlapping_constraints!(model::JuMP.Model, Xs::Vector{SVector{d, PeriodicNumber{T}}}, Rs::Vector{T}, ℓ::T, images::Vector{SVector{d, T}}, distances::Symmetric{T, Matrix{T}}) where {d, T<:AbstractFloat}
    #TODO: check if the method of distances_between_centers(centers, images) yields a more efficient function. In such case, the indices of images are also computed from the beginning
    
    N = length(Xs); # number of particles
    # access the design variables of 'model' (S corresponds to the d*N array of the displacements; the inflation factor determines how much the radius of each particle is increased)
    S⃗ = model[:S];  Γ = model[:inflation_factor]

    # Because we're using ordered pairs [i.e. (i,j) with (i<j)] to enumerate the constraints, only N-1 entries are needed in the following arrays
    constraints = Vector{Vector{ConstraintRef}}(undef, N-1) # initialization of vector for storing the lists of constraints
    nearby_particles_list = Vector{Vector{Int64}}(undef, N-1) # initialization of vector for storing the lists of particles inducing constraints

    # We will only consider nearby particles --i.e. particles inside a neighbourhood of size 'ℓ' around-- as possible neighbours
    # And hence, the non-overlapping constraint will only be applied to this smaller set of particles
    for i in 1:N-1::Int
        #NB. We reduce the number of constraints, by considering neighbours with index >i.
        neigh_parts = collect(i+1:N)[distances[i+1:N, i] .<= ℓ]  # indices of particles whose distance to i is smaller than ℓ
        Nn = length(neigh_parts) # number of particles nearby i
        nearby_particles_list[i] = neigh_parts; # store the list of nearby particles
        
        constrs_part_i = Vector{ConstraintRef}(undef, Nn) # temporary array to store the constraints associated to particle i
        # this loop adds the relevant constraints of particle i to 'model' AND stores them in 'constrs_part_i'
        for (num, j) in enumerate(neigh_parts)
            # to determine the vector joining Xi with a given neighbour we used the MIC distance
            Xij, mic_image = MIC_vector(Xs[i], Xs[j], images)
            constrs_part_i[num] = @constraint(model, 2*dot(Xij, S⃗[:, j]-S⃗[:, i]) + Γ*(Rs[i]+Rs[j])^2 <= norm(Xij)^2 )
        end
        constraints[i] = constrs_part_i #store the list of constraints of particle i as the i-th entry of 'constraints'

    end
    return constraints, nearby_particles_list
end


"""
    solve_LP_instance!(LP_model::Model, Xs::Vector{SVector{d, PeriodicNumber{T}}}, R::T, ℓ::T, s_bound::T, images::Vector{SVector{d, T}}, distances::Symmetric{T, Matrix{T}}; verbose_LP_info::Bool=false, dig_S::Int=13) where {d, T<:Float64}
    solve_LP_instance!(LP_model::Model, Xs::Vector{SVector{d, PeriodicNumber{T}}}, Rs::Vector{T}, ℓ::T, s_bound::T, images::Vector{SVector{d, T}}, distances::Symmetric{T, Matrix{T}}; verbose_LP_info::Bool=false, dig_S::Int=4) where {d, T<:Float64}
    
Optimize a LP instance defined by the particles positions `Xs` and radius `R` (or radii `Rs` in polydisperse systems).

The function needs as input a JuMP model (`LP_model`), then it defines the required variables (`S⃗` and `Γ`), 
assigns all the non-overlapping constraints, the bounds on the displacements, and performs the optimization.
The constraints are assigned through `add_non_overlapping_constraints!`, so only ordered 
pairs are considered when assigning constraints. This implies that each constraint (and the 
associated indices of particles) is contained only once in the output arrays. It also 
implies that the cutoff distance, `ℓ`, along with the set of virtual `images` should be given as input.
For performance reasons, this function works assuming `distances` is a matrix whose entries are the *updated*
distances between all particles pairs.

# Output
1. The value of the optimal displacements.
2. The optimal value of the inflation factor
3. A vector containing the list of constraints of each particle.
4. A vector containing the list of particles' indices that induce constraints on each particle.
5. The time required to perform the LP optimization.

# Other Arguments
1. 

# Keyword arguments
- `verbose_LP_info::Bool=false`: a boolean to control whether or not to print info of ℓ and `sbound`
- `dig_S::Int=5`: Number of significant digits to use when printing info about the LP optimization


See also [`add_non_overlapping_constraints!`](@ref), [`@constraint`](https://jump.dev/JuMP.jl/stable/reference/constraints/#JuMP.@constraint), [`optimize!`](https://jump.dev/JuMP.jl/stable/reference/solutions/#JuMP.optimize!),
[`produce_jammed_configuration!`](@ref), [`fine_tune_forces!`](@ref).
"""
function solve_LP_instance!(LP_model::Model, Xs::Vector{SVector{d, PeriodicNumber{T}}}, R::T, ℓ::T, s_bound::T, images::Vector{SVector{d, T}}, distances::Symmetric{T, Matrix{T}};
    verbose_LP_info::Bool=false, dig_S::Int=4) where {d, T<:Float64}
        
    N = length(Xs); # number of particles
   
    @variable(LP_model, inflation_factor>=1) # Declare inflation factor design variable
    @variable(LP_model, -s_bound <= S[μ=1:d, i=1:N] <= s_bound); # Declare displacements design variables (BOUNDED)
    fix.(S[:, N], 0.0; force=true) # Fix one displacement vector so the center of mass remains constant
    
    # add the full list of non-overlapping to the model
    constraints, possible_neighs  = add_non_overlapping_constraints!(LP_model, Xs, R, ℓ, images, distances) 
    
    @objective(LP_model, Max, 1*inflation_factor) # define the objective of the model (i.e. maximize the inflation factor)
    optimize!(LP_model) # solve the LP instance, i.e. find the optimal solution (consistent with the constraints imposed above)
    t_solve = solve_time(LP_model) #store the  time required for the optimization
    
    Γ = JuMP.value(inflation_factor)
    status = termination_status(LP_model)
    S⃗ = JuMP.value.(S)
    max_Si_print = round(maximum(abs.(S⃗)), sigdigits=dig_S)
    
    # if verbose option is true, show info of the values of ℓ and bound on displacements used when creating the LP instance
    if verbose_LP_info
        println("LP instance generated with (ℓ, ℓ/R) = ", round.([ℓ, ℓ/R], sigdigits=dig_S), "\t and bound on (|s_i|, |s_i|/R) = ", round.([s_bound, s_bound/R], sigdigits=dig_S))
        println("Status: ", status, "; \t time solve LP instance (s) = ", round(t_solve, digits=4) )
        printstyled("\t Optimal values: \t √Γ-1 = ", sqrt(Γ)-1, ",\t max |sᵢ| (all particles) = ", max_Si_print, "\n", bold=true)
    end
    
    return S⃗, Γ, constraints, possible_neighs, t_solve
end

function solve_LP_instance!(LP_model::Model, Xs::Vector{SVector{d, PeriodicNumber{T}}}, Rs::Vector{T}, ℓ::T, s_bound::T, images::Vector{SVector{d, T}}, distances::Symmetric{T, Matrix{T}};
    verbose_LP_info::Bool=false, dig_S::Int=4) where {d, T<:Float64}

    N = length(Xs); # number of particles
        
    @variable(LP_model, inflation_factor>=1) # Declare inflation factor design variable
    @variable(LP_model, -s_bound <= S[μ=1:d, i=1:N] <= s_bound); # Declare displacements design variables (BOUNDED)
    fix.(S[:, N], 0.0; force=true) # Fix one displacement vector so the center of mass remains constant
    
    # add the full list of non-overlapping to the model
    constraints, possible_neighs  = add_non_overlapping_constraints!(LP_model, Xs, Rs, ℓ, images, distances) 
    
    @objective(LP_model, Max, 1*inflation_factor) # define the objective of the model (i.e. maximize the inflation factor)
    optimize!(LP_model) # solve the LP instance, i.e. find the optimal solution (consistent with the constraints imposed above)
    t_solve = solve_time(LP_model) #store the  time required for the optimization
    
    Γ = JuMP.value(inflation_factor)
    status = termination_status(LP_model)
    S⃗ = JuMP.value.(S)
    max_Si_print = round(maximum(abs.(S⃗)), sigdigits=dig_S)

    # if verbose option is true, show info of the values of ℓ and bound on displacements used when creating the LP instance
    if verbose_LP_info
        R_max = maximum(Rs)
        println("LP instance generated with (ℓ, ℓ/Rₘₐₓ) = ", round.([ℓ, ℓ/R_max], sigdigits=dig_S), "\t and bound on (|s_i|, |s_i|/Rₘₐₓ) = ", round.([s_bound, s_bound/R_max], sigdigits=dig_S))
        println("Status: ", status, "; \t time solve LP instance (s) = ", round(t_solve, digits=4) )
        printstyled("\t Optimal values: \t √Γ-1 = ", sqrt(Γ)-1, ",\t max |sᵢ| (all particles) = ", max_Si_print, "\n", bold=true)
    end
    
    return S⃗, Γ, constraints, possible_neighs, t_solve
end


"""
    network_of_contacts(Xs::Vector{SVector{d, PeriodicNumber{T}}}, constraints::Vector{Vector{ConstraintRef}}, neighbours_list::Vector{Vector{Int64}}, images::Vector{SVector{d, T}}) where {d, T<:Float64}

Construct the set of contact vectors, forces magnitudes and list of interacting particles of the full configuration.

In a first loop, the pairs of "interacting" particles are identified from the *non-zero* dual variables 
associated to the 'constraints' that are saturated, *once the LP instance has been optimized*.
So every time `shadow_price` returns a non-zero value, an ordered pair of interacting particles is 
stored, with the corresponding value of the dual variable identified as the contact force.
Finally, the corresponding contact vector for the pair is obtained by calling `MIC_vector` on 
the positions ('Xs') of the particles involved. This is why 'images' should also be provided 
so the MIC contact vector is computed more rapidly.

In the second loop the contact vectors, forces, etc. of the complementary pairs are stored. 
That is, if in the first loop only pairs (i,j) with i<j are considered, in the second one 
we assume i>j.

# Output
- `all_contact_vectors`: A vector whose elements are vectors containing SVector{d,T} entries. So the i-th element is the set of contact vectors of the i-th particle.
- `forces_dual`: A Vector{Vector{Float64}} containing the forces magnitudes acting on each particle. So its i-th element is the list of the magnitude of the forces acting on particle i.
- `particles_dual_contact`: A Vector{Vector{Int64}} containing the indices of particles in contact with each particle. So its i-th element is the list of indices of particles in contact with the i-th particle.

See also [`add_non_overlapping_constraints!`](@ref), [`@constraint`](https://jump.dev/JuMP.jl/stable/reference/constraints/#JuMP.@constraint), [`optimize!`](https://jump.dev/JuMP.jl/stable/reference/solutions/#JuMP.optimize!), [`solve_LP_instance!`](@ref).
"""
function network_of_contacts(Xs::Vector{SVector{d, PeriodicNumber{T}}}, constraints::Vector{Vector{ConstraintRef}}, neighbours_list::Vector{Vector{Int64}}, images::Vector{SVector{d, T}}; zero_force::T=default_tol_zero_forces) where {d, T<:AbstractFloat}
    
    N = length(Xs) ; 

    particles_dual_contact = Vector{Vector{Int64}}(undef, N)  # particles linked by dual variables
    forces_dual = Vector{Vector{Float64}}(undef, N) # values of forces magnitudes, from dual variables

    all_contacts = Vector{Vector{SVector{d, T}}}(undef, N)

    # In this first loop we construct the set of contact forces, neighbours, etc. as ordered pairs ([i,j] with i<j)
    for i in 1:N-1::Int
        neighs=neighbours_list[i]; #list of particles in neighbourhood of 'i'; the particles in this array are only *possible* neighbours
        real_contacts = findall(x-> x>zero_force, shadow_price.(constraints[i])) # indices of non-zero dual variables of the respective list of constraints of particle i
        m=length(real_contacts);  # number of contacts with finite dual variable 
    
        particles_dual_contact[i] = neighs[real_contacts] # storing the indices of real contacts
        forces_dual[i] = shadow_price.(constraints[i])[real_contacts]

        contact_vecs = Vector{SVector{d, T}}(undef, m)

        for (ind, j) in enumerate(particles_dual_contact[i])
            contact_vecs[ind], = MIC_vector(Xs[i], Xs[j], images)
        end
        all_contacts[i] = contact_vecs
    end
    # This ends the process of obtaining  the list of contacts (as ordered pairs [i,j], for i<j). This means 'reaction' forces are NOT counted
    # Create the (empty) lists of contact vectors, forces,etc. of particle N,  so they can be referenced to, in the next loop
    all_contacts[N] = SVector{d, T}[] ; particles_dual_contact[N] = Vector{Int64}[] ;  forces_dual[N] = Vector{Float64}[]

    # Here we also obtain the contacts of particles for which i>j; i.e. counting reaction forces
    for i in 1:N
        for j in (1:i-1)
            f_ind = findfirst(isequal(i), particles_dual_contact[j])
            if f_ind !== nothing
                push!(particles_dual_contact[i], j) 
                push!(forces_dual[i], forces_dual[j][f_ind]) 

                contact_vec, img_ind = MIC_vector(Xs[i], Xs[j], images)
                push!(all_contacts[i], contact_vec)
            end
        end
    end
    return all_contacts, forces_dual, particles_dual_contact
end

function network_of_contacts(Xs::Vector{SVector{d, PeriodicNumber{T}}}, Rs::Vector{T}, constraints::Vector{Vector{ConstraintRef}}, neighbours_list::Vector{Vector{Int64}}, images::Vector{SVector{d, T}}; zero_force::T=default_tol_zero_forces) where {d, T<:AbstractFloat}
    
    N = length(Xs) ; 

    particles_dual_contact = Vector{Vector{Int64}}(undef, N)  # particles linked by dual variables
    forces_dual = Vector{Vector{Float64}}(undef, N) # values of forces magnitudes, from dual variables

    all_contacts = Vector{Vector{SVector{d, T}}}(undef, N)

    # In this first loop we construct the set of contact forces, neighbours, etc. as ordered pairs ([i,j] with i<j)
    for i in 1:N-1::Int
        neighs=neighbours_list[i]; #list of particles in neighbourhood of 'i'; the particles in this array are only *possible* neighbours
        real_contacts = findall(x-> x>zero_force, shadow_price.(constraints[i])) # indices of non-zero dual variables of the respective list of constraints of particle i
        m=length(real_contacts); #N_contacts[i] = m # number of dual contacts (storing it)
    
        particles_dual_contact[i] = neighs[real_contacts] # storing the indices of real contacts
        λs_i = shadow_price.(constraints[i])[real_contacts]

        contact_vecs = Vector{SVector{d, T}}(undef, m)

        for (ind, j) in enumerate(particles_dual_contact[i])
            contact_vec, img_ind = MIC_vector(Xs[i], Xs[j], images)
            contact_vecs[ind] = normalize(contact_vec)
            λs_i[ind] *= (Rs[i]+Rs[j])
        end
        all_contacts[i] = contact_vecs
        forces_dual[i] = λs_i
    end
    # This ends the process of obtaining  the list of contacts (as ordered pairs [i,j], for i<j). This means 'reaction' forces are NOT counted
    # Create the (empty) lists of contact vectors, forces,etc. of particle N,  so they can be referenced to, in the next loop
    all_contacts[N] = SVector{d, T}[] ; particles_dual_contact[N] = Vector{Int64}[] ;  forces_dual[N] = Vector{Float64}[]

    # Here we also obtain the contacts of particles for which i>j; i.e. counting reaction forces
    for i in 1:N
        for j in (1:i-1)
            f_ind = findfirst(isequal(i), particles_dual_contact[j])
            if f_ind !== nothing
                push!(particles_dual_contact[i], j) 
                push!(forces_dual[i], forces_dual[j][f_ind]) 

                contact_vec, img_ind = MIC_vector(Xs[i], Xs[j], images)
                push!(all_contacts[i], normalize(contact_vec))
            end
        end
    end
    return all_contacts, forces_dual, particles_dual_contact
end


@doc raw"""
    MonoPacking(Xs::Vector{SVector{d, PeriodicNumber{T}}}, constraints::Vector{Vector{ConstraintRef}}, neighbours_list::Vector{Vector{Int64}}, R::T, images::Vector{SVector{d, T}}, jammed::Bool=false; <keyword arguments>) where {d, T<:AbstractFloat}

Construct a `MonoPacking` from the set of particles' position ('Xs'), set of all constraints defined in the LP model ('constraints'), list of *possible* neighbours ('neighbours_list'), and virtual images ('images') needed for the MIC contact vectors. "
"""
MonoPacking(Xs::Vector, constraints::Vector, neighbours::Vector, R::Real)
function MonoPacking(Xs::Vector{SVector{d, PeriodicNumber{T}}}, constraints::Vector{Vector{ConstraintRef}}, neighbours_list::Vector{Vector{Int64}}, R::T, images::Vector{SVector{d, T}}, jammed::Bool=false; tol_mechanical_equilibrium::T=default_tol_force_equilibrium, zero_force::T=default_tol_zero_forces, verbose::Bool=true) where {d, T<:AbstractFloat}
    
    all_contacts, forces_dual, particles_dual_contact = network_of_contacts(Xs, constraints, neighbours_list, images, zero_force=zero_force)
        
    particles = [MonoParticle(Xs[i], all_contacts[i], forces_dual[i], particles_dual_contact[i]) for i in eachindex(Xs)]

    return MonoPacking(particles, R, jammed; tol_mechanical_equilibrium=tol_mechanical_equilibrium, verbose=verbose)
end

@doc raw"""
    PolyPacking(Xs::Vector{SVector{d, PeriodicNumber{T}}}, Rs::Vector{T}, constraints::Vector{Vector{ConstraintRef}}, neighbours_list::Vector{Vector{Int64}}, R::T, images::Vector{SVector{d, T}}, jammed::Bool=false; <keyword arguments>) where {d, T<:AbstractFloat}

Construct a `PolyPacking` from the set of particles' position ('Xs'), radii ('Rs'), set of all constraints defined in the LP model ('constraints'), list of *possible* neighbours ('neighbours_list'), and virtual images ('images') needed for the MIC contact vectors. "
"""
function PolyPacking(Xs::Vector{SVector{d, PeriodicNumber{T}}}, constraints::Vector{Vector{ConstraintRef}}, neighbours_list::Vector{Vector{Int64}}, Rs::Vector{T}, images::Vector{SVector{d, T}}, jammed::Bool=false; tol_mechanical_equilibrium::T=default_tol_force_equilibrium, zero_force::T=default_tol_zero_forces, verbose::Bool=true) where {d, T<:AbstractFloat}
    
    all_contacts, forces_dual, particles_dual_contact = network_of_contacts(Xs, Rs, constraints, neighbours_list, images, zero_force=zero_force)

    particles = [Particle(Xs[i], Rs[i], all_contacts[i], forces_dual[i], particles_dual_contact[i]) for i in eachindex(Xs)]

    return PolyPacking(particles, jammed; tol_mechanical_equilibrium=tol_mechanical_equilibrium, verbose=verbose)
end

"""
    update_packing_forces!(Packing::MonoPacking{d,T}, constraints::Vector{Vector{ConstraintRef}}, neighbours_list::Vector{Vector{Int64}}, images::Vector{SVector{d, T}}; tol_mechanical_equilibrium::Float64=default_tol_force_equilibrium) where {d, T<:Float64}
    update_packing_forces!(Packing::PolyPacking{d,T}, constraints::Vector{Vector{ConstraintRef}}, neighbours_list::Vector{Vector{Int64}}, images::Vector{SVector{d, T}}; tol_mechanical_equilibrium::Float64=default_tol_force_equilibrium) where {d, T<:Float64}

Update the `forces` field of each `MonoParticle` in 'Packing',  (or the analogous field of 
each `Particle` in case of a polydispersity) from a new set of 'constraints'.

The arguments needed are the ones needed to call `network_of_contacts`, which is the main 
function used here. The kwarg 'tol_mechanical_equilibrium' defines the tolerance to assess 
whether the force balance condition is satisfied in each particle, and consequently update 
the `mechanical_equilibrium` field of 'Packing'.

See also [`network_of_contacts`](@ref), [`MonoPacking`](@ref), [`PolyPacking`](@ref), [`total_force`](@ref),
[`fine_tune_forces!`](@ref).
"""
function update_packing_forces!(Packing::MonoPacking{d,T}, constraints::Vector{Vector{ConstraintRef}}, neighbours_list::Vector{Vector{Int64}}, images::Vector{SVector{d, T}}; tol_mechanical_equilibrium::Float64=default_tol_force_equilibrium, zero_force::Float64= default_tol_zero_forces) where {d, T<:Float64}

    particles = Packing.Particles
    Xs = get_positions(Packing)

    all_contacts, all_forces, all_neighbours = network_of_contacts(Xs, constraints, neighbours_list, images; zero_force=zero_force)
    for (i, particle) in enumerate(particles)
        particle.contact_vecs = all_contacts[i]
        particle.forces = all_forces[i]
        particle.neighbours = all_neighbours[i]
    end

    # test whether the mechanical equilibrium condition is satisfied for all particles, at least within a 'tol_mechanical_equilibrium' precision
    f_mismatch = maximum(norm.(total_force(Packing)))
    (f_mismatch>tol_mechanical_equilibrium) ? (mech_eq=false) : (mech_eq=true)
    # if this is not the case throw a warning,
    mech_eq || ( @warn "Force balance condition is NOT satisfied! Max force mismatch = $f_mismatch.")

    Packing.mechanical_equilibrium = mech_eq # update the field assessing the force balance
end


function update_packing_forces!(Packing::PolyPacking{d,T}, constraints::Vector{Vector{ConstraintRef}}, neighbours_list::Vector{Vector{Int64}}, images::Vector{SVector{d, T}}; tol_mechanical_equilibrium::Float64=default_tol_force_equilibrium, zero_force::Float64= default_tol_zero_forces) where {d, T<:Float64}

    particles = Packing.Particles
    Xs = get_positions(Packing)
    Rs = get_radii(Packing)

    all_contacts, all_forces, all_neighbours = network_of_contacts(Xs, Rs, constraints, neighbours_list, images; zero_force=zero_force)
    for (i, particle) in enumerate(particles)
        particle.contact_vecs = all_contacts[i]
        particle.forces = all_forces[i]
        particle.neighbours = all_neighbours[i]
    end

    # test whether the mechanical equilibrium condition is satisfied for all particles, at least within a 'tol_mechanical_equilibrium' precision
    f_mismatch = maximum(norm.(total_force(Packing)))
    (f_mismatch>tol_mechanical_equilibrium) ? (mech_eq=false) : (mech_eq=true)
    # if this is not the case throw a warning,
    mech_eq || ( @warn "Force balance condition is NOT satisfied! Max force mismatch = $f_mismatch.")

    Packing.mechanical_equilibrium = mech_eq # update the field assessing the force balance
end


"""
    fine_tune_forces!(Packing::MonoPacking{d, T}, force_mismatch::T, sqrΓ::T, ℓ0::T, images::Vector{SVector{d, T}};      
    <keyword arguments> ) where {d, T<:Float64}

Update the forces of all particles in 'Packing' so the global force balance condition is more closely satisfied.

This function is needed because when CALiPPSO has converged (or terminated due to reaching 
maximum iterations) it might happen that mechanical equilibrium is rather inaccurate. This 
is caused by the displacements in the last LP optimization. So the total force on particles 
whose position was updated is about as large as the convergence tolerance on the displacements. 
Thus, to produce a packing that satisfies force balance much more precisely an extra LP 
optimization is done (so the dual variables are recalculated with the final positions), but 
the particles' position and radius are *NOT* updated.
Also importantly, to avoid the introduction of extra constraints that might ruin force balance 
(mainly when CALiPPSO terminated without producing a jammed packing), the optimization performed 
here is done *without* bounding the displacements.

# Output
Besides updating the `forces` of each particle in 'Packing', this function produces the following output
1. `t_solve`: The time elapsed during the extra LP optimization.
2. `isostatic`: A boolean that specifies whether the updated packing is isostatic or not.
3. `Nc`: The updated number of contacts.
4. `Nnr`: The updated number of non-rattlers.

# Keyword arguments
- `solver::Symbol=:Gurobi`: the solver used to call the optimizer when creating the JuMP model.
- `solver_attributes::Dict=default_solver_attributes`: The attributes (i.e. parameters) passed to the solver after creating model, using `set_optimizer_attributes`.
- `solver_args=default_args`: The arguments passed to `Optimizer` of the chosen solver. It should be either `nothing` or a `NamedTuple`. Choose the former if testing a solver other than Gurobi, GLPK, or Hypatia.
- `thresholds::Tuple{T, T}=(5e-4, 1e-5)`: thresholds that define the different criteria to determine the radius of influence, ℓ, and the displacements' bound when calling `bounds_and_cutoff`.
- `tol_mechanical_equilibrium=default_tol_force_equilibrium`: The tolerance to determine whether force balance is fulfilled in each particle.


See also [`add_non_overlapping_constraints!`](@ref), [`optimize!`](https://jump.dev/JuMP.jl/stable/reference/solutions/#JuMP.optimize!),
[`produce_jammed_configuration!`](@ref), [`solve_LP_instance!`](@ref).
"""
function fine_tune_forces!(Packing::MonoPacking{d, T}, LP_model::Model, force_mismatch::T, ℓ::T, images::Vector{SVector{d, T}}, distances::Symmetric{T, Matrix{T}};
    tol_mechanical_equilibrium::Float64=default_tol_force_equilibrium, zero_force::T = default_tol_zero_forces) where {d, T<:Float64}

    Xs = get_positions(Packing) ; N = length(Xs)        
    R = Packing.R
    rattlers = get_rattlers(Packing)

    printstyled("\nForce balance condition not met yet; max force mismatch = ", force_mismatch, "\n", color=:yellow, bold=true)
    println("Performing a last LP optimization to improve force balance.")
   
    @variable(LP_model, inflation_factor) # Declare growth factor design variable
    @variable(LP_model, S[μ=1:d, i=1:N]); # Declare displacements design variables (UNbounded)
    fix.(S[:, N], 0.0; force=true) # Fix one displacement vector so the center of mass remains constant
    
    constraints, possible_neighs  = add_non_overlapping_constraints!(LP_model, Xs, R, ℓ, images, distances)
    
    @objective(LP_model, Max, 1*inflation_factor)
    
    println("LP instance generated with (ℓ, ℓ/R) = ", [ℓ, ℓ/R], "\t and UNbounded displacements")
    
    optimize!(LP_model)
    t_solve = solve_time(LP_model) #store the  time required for the optimization
    status = termination_status(LP_model)
    #NB: We're NOT updating the positions nor the size. This is just to fine tune the value of the dual variables so a better estimation of the forces is achieved
    update_packing_forces!(Packing, constraints, possible_neighs, images; tol_mechanical_equilibrium=tol_mechanical_equilibrium, zero_force=zero_force)

    printstyled("Max optimal displacement = ", maximum(abs.(JuMP.value.(S[:, get_non_rattlers(Packing)]))), "\t √Γ-1 = ", sqrt(JuMP.value(inflation_factor))-1, "\n", bold=true)
    println("Status of last LP optimization: ", status)

    isostatic, Nc, Nnr = is_isostatic(Packing)
    if !isostatic
        printstyled("Obtained a NON-isostatic packing...\t", color=:yellow)
        rattlers = get_rattlers(Packing)
        zs = get_coordination_number(Packing)
        println("Total z in rattlers = ", sum(zs[rattlers]) )
    end
    
    printstyled("Force mismatch after final optimization = ", maximum(norm.(total_force(Packing))), "\n\n", bold=true)

    return t_solve, isostatic, Nc, Nnr
end


function fine_tune_forces!(Packing::PolyPacking{d, T}, LP_model::Model, force_mismatch::T, ℓ::T, images::Vector{SVector{d, T}}, distances::Symmetric{T, Matrix{T}};
    tol_mechanical_equilibrium::Float64=default_tol_force_equilibrium, zero_force::T = default_tol_zero_forces) where {d, T<:Float64}

    Xs = get_positions(Packing) ; N = length(Xs)        
    Rs = get_radii(Packing)
    rattlers = get_rattlers(Packing)

    printstyled("\nForce balance condition not met yet; max force mismatch = ", force_mismatch, "\n", color=:yellow, bold=true)
    println("Performing a last LP optimization to improve force balance.")
    
    @variable(LP_model, inflation_factor) # Declare growth factor design variable
    @variable(LP_model, S[μ=1:d, i=1:N]); # Declare displacements design variables (UNbounded)
    fix.(S[:, N], 0.0; force=true) # Fix one displacement vector so the center of mass remains constant
    
    constraints, possible_neighs  = add_non_overlapping_constraints!(LP_model, Xs, Rs, ℓ, images, distances)
    
    @objective(LP_model, Max, 1*inflation_factor)
    
    println("LP instance generated with (ℓ, ℓ/Rₘₐₓ) = ", [ℓ, ℓ/maximum(Rs)], "\t and UNbounded displacements")
    
    optimize!(LP_model)
    t_solve = solve_time(LP_model) #store the  time required for the optimization
    status = termination_status(LP_model)
    #NB: We're NOT updating the positions nor the size. This is just to fine tune the value of the dual variables so a better estimation of the forces is achieved
    update_packing_forces!(Packing, constraints, possible_neighs, images; tol_mechanical_equilibrium=tol_mechanical_equilibrium, zero_force=zero_force)

    printstyled("Max optimal displacement = ", maximum(abs.(JuMP.value.(S[:, get_non_rattlers(Packing)]))), "\t √Γ-1 = ", sqrt(JuMP.value(inflation_factor))-1, "\n", bold=true)
    println("Status of last LP optimization: ", status)

    isostatic, Nc, Nnr = is_isostatic(Packing)
    if !isostatic
        printstyled("Obtained a NON-isostatic packing...\t", color=:yellow)
        rattlers = get_rattlers(Packing)
        zs = get_coordination_number(Packing)
        println("Total z in rattlers = ", sum(zs[rattlers]) )
    end
    
    printstyled("Force mismatch after final optimization = ", maximum(norm.(total_force(Packing))), "\n\n", bold=true)

    return t_solve, isostatic, Nc, Nnr
end

## This function is required because the force balance condition is monitored throughout the CALiPPSO process.
"Compute the sum of forces on each particle from the full set of contact vectors and magnitudes of contact forces."
function total_force(all_contact_vectors::Vector{Vector{SVector{d,T}}}, all_forces::Vector{Vector{T}}) where {d, T<:Float64}
    N=length(all_forces)
    total_forces = Vector{SVector{d, T}}(undef, N)
    for i in 1:N
        total_force = @SVector zeros(d)
        contact_vecs = all_contact_vectors[i]
        forces = all_forces[i]
        for (k, f) in enumerate(forces)
            total_force += f*contact_vecs[k]
        end
        total_forces[i] = total_force
    end
    return total_forces
end


#= 
This function is called, *mainly*, to monitor the process of CALiPPSO.
But identifying the rattlers is very important because since they can move a lot, even *at* jamming, the convergence criterion of 
a very small displacements should only be applied to stable particles (i.e. to non-rattlers).
=#
"Return the set of stable particles, its amount, and the coordination number of *all* particles."
function obtain_non_rattlers(constraints::Vector{Vector{ConstraintRef}}, neighbours_list::Vector{Vector{Int64}}, d::Int64; zero_force::Float64=default_tol_zero_forces)
    N = length(neighbours_list) + 1  # we need to add an extra component since only the first (N-1) particles have reference to constraints
    particles_in_contact = Vector{Vector{Int64}}(undef, N)
    coord_nums=zeros(Int64, N);  #list of coordination numbers

    for (i, neighs) in enumerate(neighbours_list)
        # neighs=neighbours_list[i]; # list of particles in neighbourhood of 'i';
        real_contacts = findall(x-> x>zero_force, shadow_price.(constraints[i])) # indices of non-zero dual variables of the respective list of constraints of particle i
        m=length(real_contacts); coord_nums[i] = m
        particles_in_contact[i] = neighs[real_contacts] # indices of particles in (linear) contact with particle i
    end

    for i in 1:N::Int
        c = 0
        for j in (1:i-1)
            (i ∈ particles_in_contact[j]) && (c+=1)
        end
        coord_nums[i]+=c
    end

    non_rattlers=collect(1:N)[coord_nums .> d]; Nnr=length(non_rattlers)

    return non_rattlers, Nnr, coord_nums
end


"For a given configuration, return the largest value of sum of forces per particle (i.e. the largest force mismatch). It also returns a small sample of the smallest forces in the configuration."
function monitor_force_balance(Xs::Vector{SVector{d, PeriodicNumber{T}}}, constraints::Vector{Vector{ConstraintRef}}, neighbours_list::Vector{Vector{Int64}}, images::Vector{SVector{d, T}}; zero_force::T=default_tol_zero_forces, sample_size::Int64=10) where {d, T<:Float64}
    contacts, forces, neighs_lists= network_of_contacts(Xs, constraints, neighbours_list, images, zero_force=zero_force)
    all_f_mags = reduce(vcat, forces); 
    (sort!(all_f_mags); unique!(all_f_mags)); 
    small_fs = all_f_mags[1:min(sample_size, length(all_f_mags))]
    f_mismatch = maximum(norm.(total_force(contacts, forces)))
    return f_mismatch, small_fs
end

function monitor_force_balance(contact_vectors::Vector{Vector{SVector{d,T}}}, forces::Vector{Vector{T}}; sample_size::Int64=10) where {d, T<:Float64}
    
    all_f_mags = reduce(vcat, forces); 
    (sort!(all_f_mags); unique!(all_f_mags)); 
    small_fs = all_f_mags[1:min(sample_size, length(all_f_mags))]
    f_mismatch = maximum(norm.(total_force(contact_vectors, forces)))
    return f_mismatch, small_fs
end


"For a given configuration, return the mismatch in the constraint imposed by ∂ℒ/∂Γ = 0 or, equivalently, ∑₍ᵢⱼ₎λᵢⱼ (σᵢⱼ)² = 1"
function monitor_Γ_constraint(forces::Vector{Vector{T}}, R::T) where {T<:Float64}
    sum_λs = 0.0
    for fs in forces
        sum_λs += sum(fs)
    end
    Γ_mismatch = 2-(4R^2)*sum_λs
    return Γ_mismatch
end


function monitor_Γ_constraint(forces::Vector{Vector{T}}, Rs::Vector{T}, contacts_lists::Vector{Vector{Int64}}) where {T<:Float64}
    sum_λs = 0.0
    for (i, fs) in enumerate(forces)
        neighs = contacts_lists[i]
        for (ind, j) in enumerate(neighs)
            if j>i
                σij = Rs[i]+Rs[j]
                sum_λs += fs[ind]*σij^2
            end
        end
    end
    Γ_mismatch = 1-sum_λs
    return Γ_mismatch
end



"""
    function update_distances!(S⃗_cum::Vector{SVector{d, T}}, s_update::T, Xs::Vector{SVector{d, PeriodicNumber{T}}}, distances::Symmetric{T, Matrix{T}}; verbose=Bool=true) where {d, T<:AbstractFloat}
    
Update (in place) the list of distances between pairs of particles, depending on the cumulative displacements `S⃗_cum` and a threshold `s_update`.
If less than 10% of the particles have a displacement, whose norm is smaller than `s_update`, only the list of distances of such particles is recalculated.
Otherwise, the distances between of all pairs is recomputed.

For convenience, if this function is called with `s_update=0.0`, the distances are always recalculated, without the need to first find the particles which have moved a lot. This feature might be useful if you want to be sure that before each LP instance all the relevant constrictions are considered.
"""
function update_distances!(S⃗_cum::Vector{SVector{d, T}}, s_update::T, Xs::Vector{SVector{d, PeriodicNumber{T}}}, distances::Symmetric{T, Matrix{T}}; verbose=Bool=true) where {d, T<:AbstractFloat}

    if s_update ==0.0 # this option is useful in case you want to update the list of pairs' distances at every iteration
        copy!(distances, distances_between_centers(Xs))
        # copy!(S⃗_cum, svectors(zeros(d, N-1), Val(d)))

    else # this is the procedure if the list of distances should be updated only when the cumulative displacements are larger than a finite 's_update'
        norms_S = norm.(S⃗_cum)
        N = length(Xs)
        indices_to_change = findall(S -> S ≥ s_update, norms_S)
        n_to_change = length(indices_to_change)
        if n_to_change ≥ N/10  # so, if at least 10% of the configuration has moved a lot, update ALL the distances
            verbose && printstyled("Updating ALL distances because ", n_to_change, " particles moved more than \ts_update = ", s_update, "\n", color=:yellow)
            verbose && printstyled("__________________________________________________________________________\n\n", color=:yellow)
            copy!(distances, distances_between_centers(Xs))
            copy!(S⃗_cum, svectors(zeros(d, N-1), Val(d)))

        elseif 0 < n_to_change #instead, if less than 10% has changed, only update the relevant particles
            inds_sort = sortperm(norms_S[indices_to_change], rev=true)[1:min(10, n_to_change)]
            
            if verbose
                printstyled("Updating distances of ", n_to_change, " particles. Threshold for updating s_update = ", s_update, "\n", color=:yellow)
                println("Some of the cumulative displacements are: \t ", norms_S[indices_to_change][inds_sort])
                println("Corresponding to particles: \t ", indices_to_change[inds_sort])
                printstyled("__________________________________________________________________________\n\n", color=:yellow)
            end
            
            distances_between_centers!(distances, Xs, indices_to_change)
            for i in indices_to_change
                S⃗_cum[i] = @SVector zeros(d)
            end

        elseif 0 == n_to_change && verbose
            inds_sort = sortperm(norms_S,rev=true)[1:min(10, N-1)]
            large_disps = norms_S[inds_sort]
            printstyled("No need to update distances.  Threshold for updating s_update = ", s_update, "\n", color=:yellow)
            println("Some of the largest cumulative displacements are: \t", large_disps)
            println("Corresponding to particles: \t", inds_sort)
            verbose && printstyled("__________________________________________________________________________\n\n", color=:yellow)
        end
    end
end


##########################################################################
##########################################################################
## Finally, we define the main function to use CALiPPSO to produce a jammed configuration of hard-spheres
##########################################################################
##########################################################################
@doc raw"""
    produce_jammed_configuration!(Xs::Vector{SVector{d, PeriodicNumber{T}}}, R::T; <keyword arguments>) where {d, T<:Float64, I<:Int64}
    produce_jammed_configuration!(Xs::Matrix{T}, R::T, L::T=1.0; <keyword arguments>) where {T<:Float64, I<:Int64}

Use CALiPPSO to generate a jammed packing from a configuration of hard-spheres with positions 'Xs' and radius 'R'.

Because this function is meant to work only with hard-spheres, clearly the initial configuration must be such that 
no overlaps are present.
The dimensionality of the system, 'd', is automatically inferred from the size of the 
`SVector`s forming the positions vector. Besides, the periodic boundary conditions are 
taken into account given that each `SVector` contains elements of type [`PeriodicNumber`](@ref). (If 
the input is of type `Matrix{T}`, it is converted to `Vector{SVector{d,T}}` before CALiPPSO begins.)

The core of this function is the so called [main loop](@ref mainloop), whose essential step consists in using [`solve_LP_instance!`](@ref)
to obtain the maximal inflation factor (``\Gamma``) and set of optimal particles' displacements 
(``\vec{\mathbf{s}}^\star = \{\mathbf{s}_i^\star\}_{i=1}^N`` --denoted as S⃗ in our scripts.) 
from a given configuration. The particles' size and position are updated, and a new LP instance is 
created and solved. (Each of these LP optimizations is considered a single iteration of CALiPPSO.) 
As explained in our paper, convergence is said to be achieved when *both*, the packing 
fraction cannot be further increased (i.e. ``\Gamma^\star=1``), and (non-rattler) particles cannot be 
further displaced (``\vec{\mathbf{s}}^\star = 0``). 
In practice, a tolerance should be allowed in both quantities, 
and this is controlled by `tol_Γ_convergence` and `tol_S_conv` as described below.
Another reason why the main loop might be terminated is that the maximal number of iterations 
has been exceeded (see below).

In case the packing thus obtained does not satisfy force balance (within a given 
precision), [`fine_tune_forces!`](@ref) is called on such packing. In this way, the final packing is 
guaranteed to be in mechanical equilibrium, within the same precision.


# Output
1. `final_packing`: A [`MonoPacking{d,T}`](@ref) corresponding to a jammed state (unless `max_iters` exceeded).
2. `conv_info`: A `convergence_info`](@ref) `struct` storing the termination status of CALiPPSO and other useful information.
3. `Γs_vs_t`: An array containing the optimal values of `\sqrt{\Gamma}` obtained after each LP optimization.
4. `smax_vs_t`: An analogous array of the largest displacement (in absolute value) of *stable* particles.
5. `iso_vs_t`: A boolean vector whose elements indicate whether the preliminary configurations were isostatic or not.


# Keyword arguments

## Arguments needed for calling [`bounds_and_cutoff`](@ref)
- `ℓ0::T=4*R`: Upper bound for the radius of influence.
- `sqrΓ0::Real=1.01`: Initialization value of ``\sqrt{\Gamma}``.
- `thresholds_bounds::Tuple{T, T}=(5e-4, 1e-5)`: Thresholds that determine the different behaviour of `bounds_and_cutoff`.
- `sbound::T=0.01`: Fraction of 'R' used for very small inflation factor; see `bounds_and_cutoff`.

## Arguments that determine `produce_jammed_configuration!` termination criteria
The list of default values is specified in [this part](@ref list-defaults) of the documentation.
- `max_iters::I=1000`: Maximum number iterations of the main loop; that is, the maximum number of LP optimizations allowed.
- `tol_Γ_convergence::T=default_tol_Γ_convergence`: determines the convergence criterion of the packing fraction as ``\sqrt{\Gamma^\star}-1 \leq `` `tol_Γ_convergence`. 
- `tol_S_convergence::T=default_tol_displacements_convergence`: determines the convergence criterion for the displacements as ``\max |\mathbf{s}_{i,\mu}^\star|_{i=1,\dots,N}^{\mu=1,\dots, d} \leq `` `<=tol_S_conv`.
- `non_iso_break::I=max_iters`: Number of *consecutive* non-isostatic solutions allowed before `produce_jammed_configuration!` terminates. The reason is that it is very likely that the final configuration will also be non-isostatic (specially if beginning from a highly compressed state). Note however that every time an isostatic configuration is obtained, this counter resets to 0.

## Arguments for controlling the behaviour of the solver/optimizer

- `solver::Module=GLPK`: The solver (*i.e.* the package or library) employed used by JuMP to solve each LP instance. 
- `solver_attributes::Dict=default_solver_attributes`: The attributes used by the solver. It's used to control some of its features, such as precision, iteration limits, etc. But it depend on which solver is used.
- `solver_args=default_args`: Arguments passed to the `solver.Optimizer` function. It is also used to control the parameters of the optimizer.

More detailed information is available [in this part](@ref changing_the_solver) of the documentation.


## Arguments to control the tolerance of mechanical equilibrium, overlaps identification, etc.
The list of default values is specified in [this part](@ref list-defaults) of the documentation.

- `tol_mechanical_equilibrium::Float64=default_tol_force_equilibrium`: tolerance to test whether a packing satisfies the force balance condition.
- `zero_force::T=default_tol_zero_forces`: threshold for identifying a force as non-zero.
- `tol_overlap::T=default_tol_overlap:`: tolerance with which overlaps are identified. That is, if particles are overlapping, but are doing so by a quantity smaller than `default_tol_overlap` (default 1e-8), no error is thrown. This rather loose tolerance is related to the maximal precision available with Gurobi.

## Arguments controlling the screen printed output
- `verbose::Bool=true`: Control whether some info about the progress of CALiPPSO and the final packing is printed out or not.
- `non_iso_warn::Bool=false`: Print some information whenever a non-isostatic *preliminary* configuration is obtained, after solving a LOP.
- `monitor_step::I=10`: How often info about the progress in the main main loop should be printed out; will only take effect if `verbose=true`. See [`print_monitor_progress`](@ref) for more information.
- `initial_monitor::I=monitor_step`: print info about the main loop progress during this amount of initial LP optimizations; will only take effect if `verbose=true`.

## Arguments for performing overlaps checks and updating list of distances
See [`check_for_overlaps`](@ref) and [`update_distances!`](@ref) for more information
- `interval_overlaps_check::I=10`: interval of LP optimizations at which it is verified that no overlaps are present in the system.
- `initial_overlaps_check::I=initial_monitor`: number of initial LP optimizations at which it is verified that no overlaps are present; given that for low density initial configurations, CALiPPSO might produce rather large displacements, it is always convenient to keep track of the overlaps in such initial stage.
- `ratio_sℓ_update::T=0.1`: This quantity determines the threshold (`s_update`) after which the list of distances of a particles is updated. Such threshold is calculated as `s_update=ratio_sℓ_update*ℓ/sqrt(d)`, with `ℓ` determined by `bounds_and_cutoff`.

"""
function produce_jammed_configuration!(Xs::Vector{SVector{d, PeriodicNumber{T}}}, R::T;
    ℓ0::T=4*R, sqrΓ0::Real=1.01, thresholds_bounds::Tuple{T, T}=(5e-4, 1e-5), sbound::T=0.01,
    optimizer::MOI.AbstractOptimizer=default_optimizer, solver_attributes::Dict=default_solver_attributes, add_bridges::Bool=false,
    max_iters::I=default_max_iterations, tol_Γ_convergence::T=default_tol_Γ_convergence, tol_S_convergence::T=default_tol_displacements_convergence, tol_mechanical_equilibrium::Float64=default_tol_force_equilibrium, zero_force::T=default_tol_zero_forces,
    tol_overlap::T=default_tol_overlap, non_iso_break::I=2*max_iters, 
    ratio_sℓ_update::T=0.1,
    verbose::Bool=true, non_iso_warn::Bool=false, monitor_step::I=10, initial_monitor::I=monitor_step, interval_overlaps_check::I=10, initial_overlaps_check::I=initial_monitor) where {d, T<:Float64, I<:Int}

    N = length(Xs)
    L = Xs[1][1].L

    #Sanity check that no overlaps are present in the initial configuration
    overlap, message, particles = check_for_overlaps(Xs, R, tol_overlap)
    overlap && error("Overlap found in the initial configuration!!\n"*message) 
    
    verbose && printstyled("\n Producing a jammed configuration of $N particles inside a $d-dimensional box of size $L:\n\n", bold=true)
    config_images = generate_system_images(d, L)
    #Initialization values
    sqrΓ=sqrΓ0; iter=0;  max_Si=1.0;  non_iso_count = 0; converged = false
    Nnr_iter = 0; 
    S⃗_cum = svectors(zeros(d, N-1), Val(d))
    distances = distances_between_centers(Xs)
    #Arrays to store the output
    constraints = Vector{Vector{ConstraintRef}}(undef, N-1); possible_neighs = Vector{Vector{Int64}}(undef, N-1)
    Γs_vs_t = zeros(max_iters+1); smax_vs_t = L*ones(max_iters+1); iso_vs_t = falses(max_iters+1); solve_ts = similar(Γs_vs_t); # +1 because it could be that an additional LP step might be needed after main loop to guarantee force balance

    #### Initialize linear model, using the optimizer passed as kwarg, as well as the attributes
    if add_bridges
        LP_model = Model(() -> optimizer) # create a model with the solver defined by 'optimizer'
    else # this is the choice by default
        LP_model = Model(() -> optimizer; add_bridges = false) # create a model with the solver defined by 'optimizer' without bridges to improve performance; it doesn't work with COSMO
        # LP_model = direct_model(optimizer) # this way of creating models is more efficient but does *not* Clp
    end

    set_optimizer_attributes(LP_model, solver_attributes...) # pass the optimizer attributes to the model
    
    set_string_names_on_creation(LP_model, false) # disable string names for better performance
   
    #######################################################################################
    # **** Here CALiPPSO begins, in order to obtain the jammed configuration ****
    # We're timing the full mainloop process and also storing the allocated memory using '@timed'.
    #######################################################################################
    exec_stats = @timed while !converged

        iter+=1
        if iter>max_iters 
            #=if the maximal number of iterations has been reached, print the current status of inflation factor (sqrΓ) and maximum displacement (max_Si).
            Then terminate the main loop. No errors are thrown. =#
            print_failed_max_iters(max_iters, sqrΓ, tol_Γ_convergence, max_Si, tol_S_convergence)
            break
        end

        # when to show info about the current LP instance, solution, and other variables to monitor the progress of CALiPPSO
        if verbose && any([iter<=initial_monitor,  iter%monitor_step==0, iter%interval_overlaps_check==0, iter<=initial_overlaps_check])
            verbose_LP_info = true
            println("Iteration: ", iter)
        else
            verbose_LP_info = false
        end

        ℓ, s_bound = bounds_and_cutoff(sqrΓ, R, ℓ0, d; thresholds=thresholds_bounds, sbound=sbound)
        s_update = ratio_sℓ_update*ℓ/sqrt(d)

        # obtain the optimal value of displacements, inflation factor, relevant constraints, and other quantities from the LP optimization
        S⃗, Γ, constraints, possible_neighs, t_solve = solve_LP_instance!(LP_model, Xs, R, ℓ, s_bound, config_images, distances; verbose_LP_info=verbose_LP_info)
        
        # track the force balance in the current iteration or not 
        if verbose && (iter<=initial_monitor || iter%monitor_step==0)
            contacts_vecs_iter, forces_iter, neighs_lists= network_of_contacts(Xs, constraints, possible_neighs, config_images, zero_force=zero_force)
            f_mismatch, small_fs = monitor_force_balance(contacts_vecs_iter, forces_iter)
            Γ_mismatch = monitor_Γ_constraint(forces_iter, R)
            # f_mismatch, small_fs = monitor_force_balance(Xs, constraints, possible_neighs, config_images; zero_force=zero_force)
            # Γ_mismatch = monitor_Γ_constraint(Xs, R, constraints, possible_neighs, config_images, zero_force=zero_force)
        end

        ### **** UPDATE CONFIGURATION WITH THE VALUES FROM LP OPTIMIZATION ****
        sqrΓ = sqrt(Γ); R*=sqrΓ; # update particles' size

        # update particles' position with optimal displacements
        for i in 1:N-1 #Since S⃗[:, N] = 0
            # Xs[i] += SVector{d}(S⃗[:, i]) 
            # It seems  that converting each column of S⃗ to an SVector takes more allocations than simply summing the elements without conversion
            Xs[i] += S⃗[:, i]
            S⃗_cum[i] += S⃗[:, i] # To monitor the total displacement of each particle
        end

        Γs_vs_t[iter] = sqrΓ; solve_ts[iter] = t_solve # store the optimal √Γ and solving time
        
        # obtain non-rattlers (important so only the displacements of these particles are considered)
        non_rattlers_iter, Nnr_iter, zs_iter = obtain_non_rattlers(constraints, possible_neighs, d; zero_force=zero_force)
        # update and store value of maximum displacement (max_Si) 
        (Nnr_iter>0) ? (max_Si = maximum(abs.(S⃗[:, non_rattlers_iter]))) : (max_Si = maximum(abs.(S⃗)))
        smax_vs_t[iter] = max_Si 

        # What follows are just checks for overlaps, verifying isostaticity, printing some output at given steps, etc.
        iso_cond = 0.5*sum(zs_iter[non_rattlers_iter]) == (d*(Nnr_iter-1)+1) # test whether the configuration is isostatic 
        iso_vs_t[iter] = iso_cond; # store such info 
       
        if verbose && (iter<=initial_monitor || iter%monitor_step==0)
            print_monitor_progress(max_Si, R, N, L, d, zs_iter, non_rattlers_iter, length.(constraints), iso_cond, Γ_mismatch, f_mismatch, small_fs)
        end

        # determine whether distances between centers should be updated
        update_distances!(S⃗_cum, s_update, Xs, distances, verbose=verbose_LP_info)

        # count how many NON-isostatic *consecutive* solutions have been produced and print info if currently NON-isostatic
        if iso_cond
            non_iso_count = 0
        else
            non_iso_count +=1
            
            if non_iso_warn && (!verbose || !(iter<=initial_monitor || iter%monitor_step==0)) # in case the non-isostatic package was obtained in an iteration for which no printed output was expected
                f_mismatch, small_fs = monitor_force_balance(Xs, constraints, possible_neighs, config_images; zero_force=zero_force)
                bound_si = bounds_and_cutoff(sqrΓ, R, ℓ0, d; thresholds=thresholds_bounds, sbound=sbound)[2]

                print_non_isostatic(d, zs_iter, non_rattlers_iter, iter, max_Si, bound_si, f_mismatch, small_fs)

            elseif verbose_LP_info
                print_non_isostatic(d, zs_iter, non_rattlers_iter)
            end
        end
        

        if non_iso_count>=non_iso_break
            printstyled("Aborting the jamming process because too many iterations (=$non_iso_break) yield NON-isostatic configurations!!\n", color=:red, bold=true)
            # converged = false
            break
        end

        # Check for overlaps during initial iterations, which might be the ones causing some troubles, or after a given interval
        # Note that these are temporary configurations.
        if iter%interval_overlaps_check==0 || iter<=initial_overlaps_check
            check_for_overlaps(Xs, R, S⃗, iter, possible_neighs; tol_overlap=tol_overlap)
        end

        if (sqrΓ-1>tol_Γ_convergence || max_Si>tol_S_convergence)
            empty!(LP_model) # so the model can be created afresh in the next iteration
        else
            converged=true
        end
    end 
    ###################
    # THIS ENDS THE MAIN LOOP OF CALiPPSO
    ###################
    t_exec = exec_stats.time # running time of the main loop
    bytes = exec_stats.bytes # memory allocated in the main loop

    # The rest of the function analyses and stores the output of the CALiPPSO process, construct the final packing, and performs some last checks.
    if converged # Check whether the final configuration corresponds to a jammed stated (i.e. convergence was achieved)
        jammed = true
        status_LP = termination_status(LP_model)
        print_converged(iter, status_LP, sqrΓ, max_Si, R, N, L, d, length.(constraints), Nnr_iter)
    else # or CALiPPSO was terminated because 'max_iters' was reached, or too many consecutive non-isostatic solutions were obtained.
        iter-=1
        jammed = false
    end
    
    # CONSTRUCT THE FINAL PACKING
    final_packing = MonoPacking(Xs, constraints, possible_neighs, R, config_images, jammed, tol_mechanical_equilibrium = tol_mechanical_equilibrium, zero_force=zero_force,verbose=verbose) 

    isostatic, Nc, Nnr = is_isostatic(final_packing) # test whether such final packing is isostatic
    force_mismatch = maximum(norm.(total_force(final_packing))) # test whether such final packing is in mechanical equilibrium

    ## if there is a (relatively) large force mismatch, perform a last LP optimizations so the force balance condition is met more accurately
    if force_mismatch>tol_mechanical_equilibrium
        empty!(LP_model) 
        ℓ, s_bound = bounds_and_cutoff(sqrΓ, R, ℓ0, d; thresholds=thresholds_bounds, sbound=sbound)

        fine_tune_stats = @timed begin
            iter+=1
            t_solve, isostatic, Nc, Nnr,  = fine_tune_forces!(final_packing, LP_model, force_mismatch, ℓ, config_images, distances; tol_mechanical_equilibrium=tol_mechanical_equilibrium, zero_force=zero_force)

            solve_ts[iter] = t_solve; Γs_vs_t[iter] = 1.0; smax_vs_t[iter] = 0.0; iso_vs_t[iter] = isostatic # store time, √Γ, and if isostatic of this last iteration.
        end
        t_exec += fine_tune_stats.time
        bytes += fine_tune_stats.bytes
    end
    empty!(LP_model)

    memory = bytes/(1024^3); # allocated memory in GB
    conv_info = convergence_info(converged, iter, t_exec, memory, solve_ts[1:iter]) # store info about the CALiPPSO process and termination status

    # when verbose output is on, print some extra info about the final packing.
    verbose && print_info_convergence(final_packing, isostatic, t_exec, memory)

    # check for overlaps in the FINAL packing
    check_for_overlaps(final_packing, iter, possible_neighs, jammed; tolerance=tol_overlap)

    return final_packing, conv_info, Γs_vs_t[1:iter], smax_vs_t[1:iter], iso_vs_t[1:iter]
end


function produce_jammed_configuration!(Xs::Matrix{T}, R::T, L::T=1.0; kwargs...) where {T<:Float64}
    produce_jammed_configuration!(PeriodicVectors(Xs, L), R; kwargs...)
end

function produce_jammed_configuration!(Xs::Vector{SVector{d, PeriodicNumber{T}}}, Rs::Vector{T};
    ℓ0::T=4.0*maximum(Rs), sqrΓ0::Real=1.01, thresholds_bounds::Tuple{T, T}=(5e-4, 1e-5), sbound::T=0.01,
    optimizer::MOI.AbstractOptimizer=default_optimizer, solver_attributes::Dict=default_solver_attributes, add_bridges::Bool=false, max_iters::I=default_max_iterations, tol_Γ_convergence::T=default_tol_Γ_convergence, tol_S_convergence::T=default_tol_displacements_convergence, tol_mechanical_equilibrium::Float64=default_tol_force_equilibrium, zero_force::T=default_tol_zero_forces,
    tol_overlap::T=default_tol_overlap, non_iso_break::I=max_iters, 
    ratio_sℓ_update::T=0.1,
    verbose::Bool=true, non_iso_warn::Bool=false, monitor_step::I=10, initial_monitor::I=monitor_step, interval_overlaps_check::I=10, initial_overlaps_check::I=initial_monitor) where {d, T<:Float64, I<:Int}

    N = length(Xs)
    L = Xs[1][1].L
    Rmax, i_max = findmax(Rs) # initial value of the largest radius, and its respective index. Note that 'i_max' should remain unchanged

    #Sanity check that no overlaps are present in the initial configuration
    overlap, message, particles = check_for_overlaps(Xs, Rs, tol_overlap)
    overlap && error("Overlap found in the initial configuration!!\n"*message) 

    verbose && printstyled("\n Producing a jammed configuration of $N particles inside a $d-dimensional box of size $L:\n\n", bold=true)
    config_images = generate_system_images(d, L)
    #Initialization values
    sqrΓ=sqrΓ0; iter=0;  max_Si=1.0;  non_iso_count = 0; converged = false
    Nnr_iter = 0; 
    S⃗_cum = svectors(zeros(d, N-1), Val(d))
    distances = distances_between_centers(Xs)
    #Arrays to store the output
    constraints = Vector{Vector{ConstraintRef}}(undef, N-1); possible_neighs = Vector{Vector{Int64}}(undef, N-1)
    Γs_vs_t = zeros(max_iters+1); smax_vs_t = L*ones(max_iters+1); iso_vs_t = falses(max_iters+1); solve_ts = similar(Γs_vs_t); # +1 because it could be that an additional LP step might be needed after convergence to guarantee force balance

     #### Initialize linear model, using the optimizer passed as kwarg, as well as the attributes
     if add_bridges
        LP_model = Model(() -> optimizer) # create a model with the solver defined by 'optimizer'
    else # this is the choice by default
        LP_model = Model(() -> optimizer; add_bridges = false) # create a model with the solver defined by 'optimizer' without bridges to improve performance; it doesn't work with COSMO
        # LP_model = direct_model(optimizer) # this way of creating models is more efficient but does *not* Clp
    end

    set_optimizer_attributes(LP_model, solver_attributes...) # pass the optimizer attributes to the model
    
    set_string_names_on_creation(LP_model, false) # disable string names for better performance

    #######################################################################################
    # **** Here CALiPPSO begins, in order to obtain the jammed configuration ****
    # We're timing the full mainloop process and also storing the allocated memory using '@timed'.
    #######################################################################################
    exec_stats = @timed while !converged

        iter+=1
        if iter>max_iters 
            #=if the maximal number of iterations has been reached, print the current status of inflation factor (sqrΓ) and maximum displacement (max_Si).
            Then terminate the main loop. No errors are thrown. =#
            print_failed_max_iters(max_iters, sqrΓ, tol_Γ_convergence, max_Si, tol_S_convergence)
            break
        end

        # when to show info about the current LP instance, solution, and other variables to monitor the progress of CALiPPSO
        if verbose && any([iter<=initial_monitor,  iter%monitor_step==0, iter%interval_overlaps_check==0, iter<=initial_overlaps_check])
            verbose_LP_info = true
            println("Iteration: ", iter)
        else
            verbose_LP_info = false
        end

        ℓ, s_bound = bounds_and_cutoff(sqrΓ, Rs[i_max], ℓ0, d; thresholds=thresholds_bounds, sbound=sbound)
        s_update = ratio_sℓ_update*ℓ/sqrt(d)

        # obtain the optimal value of displacements, inflation factor, relevant constraints, and other quantities from the LP optimization
        S⃗, Γ, constraints, possible_neighs, t_solve = solve_LP_instance!(LP_model, Xs, Rs, ℓ, s_bound, config_images, distances; verbose_LP_info=verbose_LP_info)
        
        # track the force balance in the current iteration or not 
        if verbose  && (iter<=initial_monitor || iter%monitor_step==0)
            contacts_vecs_iter, forces_iter, contacts_lists_iter= network_of_contacts(Xs, constraints, possible_neighs, config_images, zero_force=zero_force)
            f_mismatch, small_fs = monitor_force_balance(contacts_vecs_iter, forces_iter)
            Γ_mismatch = monitor_Γ_constraint(forces_iter, Rs, contacts_lists_iter)
            # f_mismatch, small_fs = monitor_force_balance(Xs, constraints, possible_neighs, config_images; zero_force=zero_force)
        end

        ### **** UPDATE CONFIGURATION WITH THE VALUES FROM LP OPTIMIZATION ****
        sqrΓ = sqrt(Γ); Rs.*=sqrΓ; # update particles' size

        # update particles' position with optimal displacements
        for i in 1:N-1 #Since S⃗[:, N] = 0
            # Xs[i] += SVector{d}(S⃗[:, i]) 
            # It seems  that converting each column of S⃗ to an SVector takes more allocations than simply summing the elements without conversion
            Xs[i] += S⃗[:, i]
            S⃗_cum[i] += S⃗[:, i] # To monitor the total displacement of each particle
        end

        Γs_vs_t[iter] = sqrΓ; solve_ts[iter] = t_solve # store the optimal √Γ and solving time
        
        # obtain non-rattlers (important so only the displacements of these particles are considered)
        non_rattlers_iter, Nnr_iter, zs_iter = obtain_non_rattlers(constraints, possible_neighs, d; zero_force=zero_force)
        # update and store value of maximum displacement (max_Si) 
        (Nnr_iter>0) ? (max_Si = maximum(abs.(S⃗[:, non_rattlers_iter]))) : (max_Si = maximum(abs.(S⃗)))
        smax_vs_t[iter] = max_Si 

        # What follows are just checks for overlaps, verifying isostaticity, printing some output at given steps, etc.
        iso_cond = 0.5*sum(zs_iter[non_rattlers_iter]) == (d*(Nnr_iter-1)+1) # test whether the configuration is isostatic 
        iso_vs_t[iter] = iso_cond; # store such info 
    
        if verbose && (iter<=initial_monitor || iter%monitor_step==0)
            print_monitor_progress(max_Si, Rs, L, d, zs_iter, non_rattlers_iter, length.(constraints), iso_cond, Γ_mismatch, f_mismatch, small_fs)
        end

        # determine whether distances between centers should be updated
        update_distances!(S⃗_cum, s_update, Xs, distances, verbose=verbose_LP_info)

        # count how many NON-isostatic *consecutive* solutions have been produced and print info if currently NON-isostatic
        if iso_cond
            non_iso_count = 0
        else
            non_iso_count +=1
            
            if non_iso_warn && (!verbose || !(iter<=initial_monitor || iter%monitor_step==0)) # in case the non-isostatic package was obtained in an iteration for which no printed output was expected
                f_mismatch, small_fs = monitor_force_balance(Xs, constraints, possible_neighs, config_images; zero_force=zero_force)
                bound_si = bounds_and_cutoff(sqrΓ, Rs[i_max], ℓ0, d; thresholds=thresholds_bounds, sbound=sbound)[2]

                print_non_isostatic(d, zs_iter, non_rattlers_iter, iter, max_Si, bound_si, f_mismatch, small_fs)

            elseif verbose_LP_info
                print_non_isostatic(d, zs_iter, non_rattlers_iter)
            end
        end
        

        if non_iso_count>=non_iso_break
            printstyled("Aborting the jamming process because too many iterations (=$non_iso_break) yield NON-isostatic configurations!!\n", color=:red, bold=true)
            # converged = false
            break
        end

        # Check for overlaps during initial iterations, which might be the ones causing some troubles, or after a given interval
        # Note that these are temporary configurations.
        if iter%interval_overlaps_check==0 || iter<=initial_overlaps_check
            check_for_overlaps(Xs, Rs, S⃗, iter, possible_neighs; tol_overlap=tol_overlap)
        end

        if (sqrΓ-1>tol_Γ_convergence || max_Si>tol_S_convergence)
            empty!(LP_model) # so the model can be created afresh in the next iteration
        else
            converged=true
        end
    end 
    ###################
    # THIS ENDS THE MAIN LOOP OF CALiPPSO
    ###################
    t_exec = exec_stats.time # running time of the main loop
    bytes = exec_stats.bytes # memory allocated in the main loop

    # The rest of the function analyses and stores the output of the CALiPPSO process, construct the final packing, and performs some last checks.
    if converged # Check whether the final configuration corresponds to a jammed stated (i.e. convergence was achieved)
        jammed = true
        status_LP = termination_status(LP_model)
        print_converged(iter, status_LP, sqrΓ, max_Si, Rs, N, L, d, length.(constraints), Nnr_iter)
    else # or CALiPPSO was terminated because 'max_iters' was reached, or too many consecutive non-isostatic solutions were obtained.
        iter-=1
        jammed = false
    end

    # CONSTRUCT THE FINAL PACKING
    final_packing = PolyPacking(Xs, constraints, possible_neighs, Rs, config_images, jammed, tol_mechanical_equilibrium = tol_mechanical_equilibrium, zero_force=zero_force, verbose=verbose) 

    isostatic, Nc, Nnr = is_isostatic(final_packing) # test whether such final packing is isostatic
    force_mismatch = maximum(norm.(total_force(final_packing))) # test whether such final packing is in mechanical equilibrium

    ## if there is a (relatively) large force mismatch, perform a last LP optimizations so the force balance condition is met more accurately
    if force_mismatch>tol_mechanical_equilibrium 
        empty!(LP_model)
        ℓ, s_bound = bounds_and_cutoff(sqrΓ, Rs[i_max], ℓ0, d; thresholds=thresholds_bounds, sbound=sbound)
        
        fine_tune_stats = @timed begin
            iter+=1          
            t_solve, isostatic, Nc, Nnr  = fine_tune_forces!(final_packing, LP_model,force_mismatch, ℓ, config_images, distances; tol_mechanical_equilibrium=tol_mechanical_equilibrium,zero_force=zero_force)

            solve_ts[iter] = t_solve; Γs_vs_t[iter] = 1.0; smax_vs_t[iter] = 0.0; iso_vs_t[iter] = isostatic # store time, √Γ, and if isostatic of this last iteration.
        end
        t_exec += fine_tune_stats.time
        bytes += fine_tune_stats.bytes
    end
    empty!(LP_model)

    memory = bytes/(1024^3); # allocated memory in GB
    conv_info = convergence_info(converged, iter, t_exec, memory, solve_ts[1:iter]) # store info about the CALiPPSO process and termination status

    # when verbose output is on, print some extra info about the final packing.
    verbose && print_info_convergence(final_packing, isostatic, t_exec, memory)

    # check for overlaps in the FINAL packing
    check_for_overlaps(final_packing, iter, possible_neighs, jammed; tolerance=tol_overlap)

    return final_packing, conv_info, Γs_vs_t[1:iter], smax_vs_t[1:iter], iso_vs_t[1:iter]
end


function produce_jammed_configuration!(Xs::Matrix{T}, Rs::Vector{T}, L::T=1.0; kwargs...) where {T<:Float64}
    produce_jammed_configuration!(PeriodicVectors(Xs, L), Rs; kwargs...)
end



"""
    network_of_contacts(packing::AbstractPacking{d, T}, normalized::Bool=true) where {d, T<:Float64}

Obtain the list of contact indices (as ordered pairs, i.e. [i, j] with j>i), the corresponding contact vectors, and magnitudes of contact forces, from a given 'packing'.

This function is not used in the main CALiPPSO function (i.e. [`produce_jammed_configuration!`](@ref)), but 
is useful for extracting the system's micro-structural information once the jammed packings have been generated.
"""
function network_of_contacts(packing::AbstractPacking{d, T}, normalized::Bool=true; only_stable::Bool=true) where {d, T<:AbstractFloat}
    particles = packing.Particles
    N = length(particles)

    isostatic, Nc, Nnr = is_isostatic(packing; only_stable=only_stable)
    non_rattlers = get_non_rattlers(packing);
    (length(non_rattlers) == Nnr) || @warn "Error occurred identifying non_rattlers, and hence also isostaticity"

    isostatic || @warn "Analysing a NON-isostatic packing"

    contact_indices = Vector{Tuple{Int64,Int64}}(undef,Nc)
    contact_vectors = Vector{SVector{d, Float64}}(undef, Nc)
    forces_magnitudes = zeros(Nc)
    contact_count = 0
    for i in non_rattlers
        P = particles[i]
        neighs = P.neighbours[P.neighbours .>i] # list of neighbours of particle i, whose index is larger than i
        forces = P.forces[P.neighbours .>i] # corresponding list of forces magnitudes 
        cvecs = P.contact_vecs[P.neighbours .>i] # list of corresponding contact vectors
        for (index, j) in enumerate(neighs)
            contact_count +=1
            contact_indices[contact_count] = (i, j)
            contact_vectors[contact_count] = cvecs[index]
            forces_magnitudes[contact_count] = forces[index]
        end
    end

    #=  Note that given that forces are obtained as the dual variables, in general, mechanical equilibrium is *only guaranteed if using the UNORMALIZED contact vectors*.
    However, by rescaling the magnitudes and vectors consistently the condition should hold, even for polydisperse packings and *away* from jamming.
    However, the constructor defined for PolyPacking in terms of the constraints, already takes care of normalizing the contact vectors.
    So in general, with such type of packings, this function should be called with normalized=false; at least with packings for which the convergence criteria were met.
    =#
    if normalized
        norm_cvecs = norm.(contact_vectors)
        forces_magnitudes .= norm_cvecs.*forces_magnitudes
        contact_vectors .= contact_vectors./norm_cvecs
    end

    return contact_indices, contact_vectors, forces_magnitudes
end



function precompile_main_function(optimizer::MOI.AbstractOptimizer=default_optimizer, solver_attributes::Dict=default_solver_attributes; add_bridges::Bool=false)
    #= =====================================================================
    =====================================================================
    USING MAIN FUNCTIONS ON SMALL SYSTEM FOR COMPILATION/FIRST CALL PURPOSES
    =====================================================================
    ===================================================================== =#
    test_model = Model(() -> optimizer) # create a model with the solver defined by 'optimizer' without bridges to improve 
    set_optimizer_attributes(test_model, solver_attributes...)
    optim_name = solver_name(test_model); 
    
    printstyled("\t The following output just refers to simple first calls for precompiling needed functions\n\n", bold=true, color=:cyan)

    printstyled("\n\nUsing CALiPPSO (with *low accuracy*) on a small, predefined system in order to compile the needed functions.\n
    \t\t Solver used: ", optim_name, "\n\n", color=:cyan)
    Nt=30; rt=0.52; dt=3; 
    Lt=4.0;
    images_comp = generate_system_images(dt, Lt)
    tol_Γ = 1e-3; tol_S=1e-2; tol_f = 10*tol_S; 
    tol_ovlp = 1e-5 # I know, this is a huge value, but it's used only to make precompilation work with ANY solver (e.g. COSMO)

    Xs_comp = Matrix(transpose([
    1.3022036173223075	3.100681668734574	2.650145868235529
    1.1155120883337997	2.414092386379111	0.3717676414091988
    0.5930622140684934	0.42734603856532427	1.3262407197673065
    3.571441842561353	2.957136483108787	2.6002633302062668
    2.3529766061581237	0.48997057275911615	3.417191033324789
    1.7112056139610727	0.8347268069230269	2.229786164196563
    3.6817595377155348	0.6725943671221577	2.4334027050690317
    1.9126951936859324	1.5811607905044598	3.460749757890367
    2.10364317266594	2.310954222951657	2.6450251464477823
    3.2140912406591706	1.9574663303701785	3.5356829974762487
    1.5789002308497055	2.218608327103566	1.664307104091713
    0.720581832264064	3.921707081207914	3.2620078999221587
    2.459219015775952	3.503789658633769	0.9321778468174386
    3.6785107203364795	3.337897259842034	0.6462417289544637
    3.4240729350390406	2.2304592814837587	1.2859711802875369
    2.200649934148873	2.428279519213076	0.38501446881727297
    0.35487777941165266	2.5227027830749185	1.7324205587874006
    0.6937597057720284	1.675212825207569	1.1950801557221888
    2.5783261715853776	0.49417638526472807	1.2809483317712305
    0.25868616273463285	2.114656556778108	3.8413865080579885
    2.269006112530807	3.5797496686786596	2.8118840862720855
    0.4016449295907414	2.092677546492193	2.715279086250721
    0.936166141399041	0.19318442951467052	0.3547261200994525
    3.4700487156587823	1.163132895735302	1.0834877206251514
    3.307514775505691	3.9075523242046204	3.036230135532204
    1.1350771834782938	1.2830440239210912	0.11294699337022074
    1.7035741859408704	3.868891311445325	1.722953156549603
    3.896005446471116	3.1189426973423595	3.607500484893032320
    3.653047120545878	0.25044725175167404	3.9730036860708076
    0.6334501855938681	1.1831285025320382	3.1868918476405756]))
    cen_comp = PeriodicVectors(Xs_comp, Lt);
    printstyled("\nThe usual output of the tested optimizers has been disabled by default. Change the value of OutputFlag (or verbose, or the corresponding field of the solver of your choice) if you want such information to be printed out.\n\n", color=:magenta)

    dists_comp = distances_between_centers(cen_comp)
    Ds_comp, Γ_comp, cons_comp, neighs_comp, = solve_LP_instance!(test_model, cen_comp, rt, 1.5, 4*rt, images_comp, dists_comp)
    obtain_non_rattlers(cons_comp, neighs_comp, dt)

    MonoPacking(cen_comp, cons_comp, neighs_comp, rt, images_comp; verbose=false)
    empty!(test_model)
    
    Jpack_comp, =  produce_jammed_configuration!(cen_comp, rt, ℓ0=Lt, initial_monitor=0, tol_Γ_convergence= tol_Γ, tol_S_convergence=tol_S, tol_mechanical_equilibrium=tol_f, tol_overlap=tol_ovlp, verbose=false, optimizer=optimizer, solver_attributes=solver_attributes, add_bridges=add_bridges, max_iters=20 )
    
    network_of_contacts(Jpack_comp)

    # # Function for polydisperse packings. I'm using the same configuration, with 'Rs' a vector of equal elements
    Jpack_comp,  =  produce_jammed_configuration!(cen_comp, rt*ones(Nt), ℓ0=Lt, initial_monitor=0, tol_Γ_convergence= tol_Γ, tol_S_convergence=tol_S, tol_mechanical_equilibrium=tol_f, tol_overlap=tol_ovlp, verbose=false, optimizer=optimizer, solver_attributes=solver_attributes, max_iters=20, add_bridges=add_bridges)

    network_of_contacts(Jpack_comp)


    printstyled("\n\n Calling the main functions once again, to compile the second method.\n\n", color=:cyan)
    produce_jammed_configuration!(Xs_comp, rt, Lt; ℓ0=Lt, initial_monitor=0, monitor_step=0, verbose=false,
     tol_Γ_convergence= tol_Γ, tol_S_convergence=tol_S, tol_mechanical_equilibrium=tol_f, tol_overlap=tol_ovlp,
     optimizer=optimizer, solver_attributes=solver_attributes, add_bridges=add_bridges, max_iters=20 )

    produce_jammed_configuration!(Xs_comp, rt*ones(Nt), Lt; ℓ0=Lt, initial_monitor=0, tol_Γ_convergence= tol_Γ, tol_S_convergence=tol_S, tol_mechanical_equilibrium=tol_f, tol_overlap=tol_ovlp, verbose=false, optimizer=optimizer, solver_attributes=solver_attributes, add_bridges=add_bridges, max_iters=20)

    printstyled("\n________________________________________________________\n\tCompilation process finished! \t (with optimizer: ", optim_name, ")\n________________________________________________________\n\n\n\n\n", color=:cyan, bold=true)
end

end