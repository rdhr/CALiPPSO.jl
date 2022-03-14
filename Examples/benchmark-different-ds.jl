using CALiPPSO
using DelimitedFiles, StatsBase, BenchmarkTools, HiGHS

const tol_optimality = CALiPPSO.default_tol_optimality

const highs_args = nothing
const highs_attributes = Dict("small_matrix_value"=>0.1*tol_optimality, "primal_feasibility_tolerance" => tol_optimality, "dual_feasibility_tolerance" => tol_optimality, "solver" => "ipm", "ipm_optimality_tolerance"=>tol_optimality,  "threads" => CALiPPSO.max_threads, "parallel"=>"on", "output_flag"=>false)
printstyled("\n______________________________________________________________________________________\n", color=:yellow)
precompile_main_function(HiGHS, highs_attributes, highs_args)


test_bench = @benchmarkable sin(1)
run(test_bench)


const N_LS = 1024 # size of configurations obtained after LS compression
const ds = [3, 4, 5] # dimensions to use
const L = 1.0 # size of each side of the system's volume

const io = IOContext(stdout)

# Testing CALiPPSO in systems of different dimensions.
# The initial conditions are configurations compressed with Lubachevsky--Stillinger up to a (reduced) pressure p>>10
for d in ds
    σ0 = readdlm("Centers-after-LS--N-$(N_LS)--d-$d.dat", skipstart=2)[1, 1] # initial diameter of the particles
    r0 = 0.5*σ0 # initial radius of particles
    Xs0 = PeriodicVectors(Matrix(readdlm("Centers-after-LS--N-$(N_LS)--d-$d.dat", skipstart=3)'), L) # initial particles' positions

    printstyled("\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", color=:blue)
    printstyled("Using CALiPPSO to jam a system of N = ", N_LS, " in d = ", d, "\t Initial φ = ", packing_fraction(d, r0, N_LS, L), "\n", bold=:true, color=:blue)
    printstyled("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", color=:blue)


    benchmark_iLP = @benchmarkable produce_jammed_configuration!(Xs, r; verbose=false, solver=HiGHS, solver_attributes=highs_attributes, solver_args=highs_args) setup=(Xs = copy($Xs0); r=copy($r0) )
    res_benchmark = run(benchmark_iLP, seconds = 3600, samples=20)

    printstyled("\n\n - - - o o o - - - o o o - - - o o o - - - o o o - - - o o o - - - o o o - - - o o o - - - \n\n", color=:green)
    println("Results of benchmarking CALiPPSO with 20 samples:\n\n")
    show(io, MIME("text/plain"), res_benchmark)
    printstyled("\n\n - - - o o o - - - o o o - - - o o o - - - o o o - - - o o o - - - o o o - - - o o o - - - \n\n", color=:green)

    # io2 = open("Benchmarks--d-$d-after-LS--changing-solver.txt", "w")
    io2 = open("Benchmarks--d-$d-after-LS--base.txt", "w")
    show(io2, MIME("text/plain"), res_benchmark)
    close(io2)
end
