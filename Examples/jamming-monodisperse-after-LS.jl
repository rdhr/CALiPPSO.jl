include("../src/CALiPPSO.jl")
using .CALiPPSO  
using DelimitedFiles, Statistics

const N_LS = 1024 # size of configurations obtained after LS compression
const ds = [3, 4, 5] # dimensions to use
const L = 1.0 # size of each side of the system's volume

# Testing CALiPPSO in systems of different dimensions.
# The initial conditions are configurations compressed with Lubachevsky--Stillinger up to a (reduced) pressure p>>10
for d in ds
    σ0 = readdlm("Centers-after-LS--N-$(N_LS)--d-$d.dat", skipstart=2)[1, 1] # initial diameter of the particles
    r0 = 0.5*σ0 # initial radius of particles
    Xs0 = PeriodicVectors(Matrix(readdlm("Centers-after-LS--N-$(N_LS)--d-$d.dat", skipstart=3)'), L) # initial particles' positions

    printstyled("\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", color=:blue)
    printstyled("Using CALiPPSO to jam a system of N = ", N_LS, " in d = ", d, "\t Initial φ = ", packing_fraction(d, r0, N_LS, L), "\n", bold=:true, color=:blue)
    printstyled("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", color=:blue)


    @time jammed_packing, info_convergence, Γs_vs_t, smax_vs_t, iso_vs_t = produce_jammed_configuration(Xs0, r0; verbose=true, tol_Γ_convergence=1e-14)
    println("Evolution of convergence criteria:")
    sqrt_gms = max.(Γs_vs_t, 1) .-1
    println("\t √Γ-1 = \t", sqrt_gms)
    println("\t max sᵢ = \t", smax_vs_t)
    println("_______________________________________________________________________________________\n\n")
    times = info_convergence.times_LP_optim
    println("Info about LP times")
    println(times)
    println("(min, avg ± std, max) LP times:\t", minimum(times), ", \t", mean(times), " ± ", std(times), ", \t", maximum(times))
    println("_______________________________________________________________________________________\n\n\n")

    
    network_of_contacts(jammed_packing)



end