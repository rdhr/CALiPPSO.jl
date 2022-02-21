############## 
### Testing CALiPPSO in the same system, with different solvers
############## 
include("../src/CALiPPSO.jl")
using .CALiPPSO  
using DelimitedFiles, Statistics


const N, d = 1024, 3 # size and dimensionality of the configuration  to be used
const L = 1.0 # size of each side of the system's volume

# Note that this value is x100 larger than the default one (that can be used with Gurobi). Some solvers need this smaller precision
const tol_Γ_convergence = 1e-12 
const tol_overlap = CALiPPSO.default_tol_overlap
const tol_zero_forces = CALiPPSO.default_tol_zero_forces
const default_tol_displacements = CALiPPSO.default_tol_displacements
const tol_optimality = CALiPPSO.default_tol_optimality
const max_threads = CALiPPSO.max_threads

const solvers = [:Gurobi, :HiGHS, :Clp, :GLPK, :Hypatia, :COSMO] # array of solver's names (or symbols)
using Gurobi, HiGHS, Clp, GLPK, Hypatia, COSMO


println("This script will test the following solvers: ",  map(x->x*", ", string.(solvers))...)
println("But before the test is executed, `produce_jammed_configuration` will be compiled with each of them.\n\n")

##############################
### Now let's define the arguments and attributes to be used for each solver, and create an empty model for initialization
##############################
# Xs_comp = CALiPPSO.Xs_comp; Lt = CALiPPSO.Lt; rt=CALiPPSO.rt

# Gurobi:
const grb_args = Gurobi.Env()
const grb_attributes = Dict("OutputFlag" => 0, "FeasibilityTol" => tol_optimality, "OptimalityTol" => tol_optimality, "Method" => 3, "Threads" => max_threads)
precompile_main_function(Gurobi, grb_attributes, grb_args)


# HiGHS
const highs_args = nothing
const highs_attributes = Dict("small_matrix_value"=>0.1*tol_optimality, "primal_feasibility_tolerance" => tol_optimality, "dual_feasibility_tolerance" => tol_optimality, "solver" => "ipm", "ipm_optimality_tolerance"=>tol_optimality,  "highs_max_threads" => max_threads, "parallel"=>"on", "output_flag"=>false)
printstyled("\n______________________________________________________________________________________\n", color=:yellow)
precompile_main_function(HiGHS, highs_attributes, highs_args)
# cen_comp = PeriodicVectors(Xs_comp, Lt)
# produce_jammed_configuration(cen_comp, rt, ℓ0=Lt, max_iters=1, verbose=false, solver=Gurobi)
# printstyled("Performing first call (of a single iteration) of main function on compilation system using HiGHS solver\n", color=:yellow)
# produce_jammed_configuration(cen_comp, rt, ℓ0=Lt, max_iters=1, verbose=false, solver=HiGHS, solver_attributes=highs_attributes, solver_args=highs_args);



# Clp
const clp_args = nothing
const clp_attributes = Dict("PrimalTolerance"=>tol_optimality, "DualTolerance" => tol_optimality, "LogLevel" => 0, "SolveType" => 5)
# cen_comp = PeriodicVectors(Xs_comp, Lt)
printstyled("\n______________________________________________________________________________________\n", color=:yellow)
precompile_main_function(Clp, clp_attributes, clp_args)
# printstyled("Performing first call (of a single iteration) of main function on compilation system using Clp solver\n", color=:yellow)
# produce_jammed_configuration(cen_comp, rt, ℓ0=Lt, max_iters=1, verbose=false, solver=Clp, solver_attributes=clp_attributes, solver_args=clp_args);


# GLPK
const glpk_args = (want_infeasibility_certificates=false, method=GLPK.MethodEnum(0) ) 
const glpk_attributes = Dict("msg_lev"=>GLPK.GLP_MSG_OFF, "tol_bnd"=>tol_optimality, "tol_dj"=>tol_optimality)
# cen_comp = PeriodicVectors(Xs_comp, Lt)
printstyled("\n______________________________________________________________________________________\n", color=:yellow)
precompile_main_function()
# printstyled("Performing first call (of a single iteration) of main function on compilation system using GLPK solver\n", color=:yellow)
# produce_jammed_configuration(cen_comp, rt, ℓ0=Lt, max_iters=1, verbose=false, solver=GLPK, solver_attributes=glpk_attributes, solver_args=glpk_args);


# Hypatia
const hypa_args = (verbose = false, tol_abs_opt = tol_optimality, tol_feas = 0.1*tol_overlap, tol_infeas = 0.1*tol_overlap)
const hypa_attributes = Dict()
# cen_comp = PeriodicVectors(Xs_comp, Lt)
printstyled("\n______________________________________________________________________________________\n", color=:yellow)
precompile_main_function(Hypatia, hypa_attributes, hypa_args)
# printstyled("Performing first call (of a single iteration) of main function on compilation system using Hypatia solver\n", color=:yellow)
# produce_jammed_configuration(cen_comp, rt, ℓ0=Lt, max_iters=1, verbose=false, solver=Hypatia, solver_attributes=hypa_attributes, solver_args=hypa_args);


# COSMO
const cosmo_args = nothing
const cosmo_attributes = Dict("verbose"=>false, "max_iter"=>30000, "eps_abs"=>1e3*tol_optimality, "adaptive_rho"=>false, "eps_prim_inf"=>1e3*tol_optimality, "eps_dual_inf"=>1e3tol_optimality)
# cen_comp = PeriodicVectors(Xs_comp, Lt)
printstyled("\n______________________________________________________________________________________\n", color=:yellow)
precompile_main_function(COSMO, cosmo_attributes, cosmo_args)
# printstyled("Performing first call (of a single iteration) of main function on compilation system using COSMO solver\n", color=:yellow)
# produce_jammed_configuration(cen_comp, rt, ℓ0=Lt, max_iters=1, verbose=false, solver=COSMO, solver_attributes=cosmo_attributes, solver_args=cosmo_args);


# # Tulip
# const tulip_args = nothing
# const tulip_attributes = Dict("OutputLevel"=>0)
# cen_comp = PeriodicVectors(Xs_comp, Lt)
# printstyled("\n______________________________________________________________________________________\n", color=:yellow)
# printstyled("Performing first call (of a single iteration) of main function on compilation system using Tulip solver\n", color=:yellow)
# produce_jammed_configuration(cen_comp, rt, ℓ0=Lt, max_iters=1, verbose=false, solver=:Tulip, solver_attributes=tulip_attributes, solver_args=tulip_args);


#### create arrays of these variables, for easier calling
const solvers_args = [grb_args, highs_args, clp_args, glpk_args, hypa_args, cosmo_args]
const solvers_attributes = [grb_attributes, highs_attributes,clp_attributes ,glpk_attributes, hypa_attributes, cosmo_attributes]


## Array of packings to perform later comparison
all_packings = Vector{MonoPacking{d, Float64}}(undef, length(solvers))

# The initial conditions are configurations compressed with Lubachevsky--Stillinger up to a (reduced) pressure p>>10
#### read the initial configuration on which the test is performed
σ0 = readdlm("Centers-after-LS--N-$(N)--d-$d.dat", skipstart=2)[1, 1] # initial diameter of the particles
r0_common = 0.5*σ0 # initial radius of particles
Xs_common = PeriodicVectors(Matrix(readdlm("Centers-after-LS--N-$(N)--d-$d.dat", skipstart=3)'), L) # initial particles' positions


printstyled("\n\n O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O\n", color=:blue, bold=true)
printstyled("Applying CALiPPSO for jamming a system of N = ", N, " in d = ", d, "\t Initial φ = ", packing_fraction(d, r0_common, N, L), "\n", bold=:true, color=:blue)
printstyled("O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O --- O\n", color=:blue, bold=true)

for (ns, solver) in enumerate(solvers)
    Xs0 = 1.0.*Xs_common; r0=1*r0_common
    
    if  Symbol(solver)==:COSMO 
        # continue
        overlap_tolerance = 1e3*tol_overlap
        zero_force = 1e-5
        # S_conv = 1e-9
        Γ_conv = tol_Γ_convergence
    elseif Symbol(solver)==:Hypatia
        # continue
        overlap_tolerance = tol_overlap
        zero_force = 1e-5
        # S_conv = 1e-9
        Γ_conv = 1e-10
    else
        # continue
        overlap_tolerance = tol_overlap
        zero_force = tol_zero_forces
        # S_conv = default_tol_displacements
        Γ_conv = tol_Γ_convergence
    end
    S_conv = default_tol_displacements

    printstyled("\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", color=:magenta, bold=true)
    printstyled("\t\t\tUSING SOLVER: ", Symbol(solver), "\n", bold=:true, color=:magenta)
    println("\tConvergence criteria: tolerance √Γ-1 = ",Γ_conv, "\t tolerance |sᵢ| = ", S_conv)
    printstyled("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", color=:magenta, bold=true)

    solver_args = solvers_args[ns]; solver_attrs = solvers_attributes[ns]

    @time jammed_packing, info_convergence, Γs_vs_t, smax_vs_t, iso_vs_t = produce_jammed_configuration(Xs0, r0; verbose=true, 
    solver=eval(solver), solver_args=solver_args, solver_attributes=solver_attrs, tol_Γ_convergence=Γ_conv, tol_S_convergence = S_conv,
    tol_overlap=overlap_tolerance, initial_monitor=30, zero_force=zero_force, max_iters=50) 

    println("_______________________________________________________________________________________\n\n")
    times = info_convergence.times_LP_optim
    println("Info about LP times")
    println(times)
    println("(min, avg ± std, max) LP times:\t", minimum(times), ", \t", mean(times), " ± ", std(times), ", \t", maximum(times))
    println("_______________________________________________________________________________________\n\n\n")

    
    network_of_contacts(jammed_packing)

    all_packings[ns] = jammed_packing

end

Δpacks = map(x-> difference_in_packings(all_packings[1], x), all_packings)
Rref = all_packings[1].R
ΔRs = Rref .- getfield.(all_packings, :R)
printstyled("Comparing packing differences (only considering stable particles) , with respect to: ", solvers[1], "\n", color=:green)
for (ns, solver) in enumerate(solvers)
    println("Solver: ", solver, "\t mean and std Δ = ", Δpacks[ns], "\t and of difference of R = ", ΔRs[ns])
end
