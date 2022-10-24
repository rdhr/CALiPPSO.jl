#########################################################################################
### Some functions to print output and monitor the progress and termination status of CALiPPSO
#########################################################################################
#Memory usage in GB. It's not used by CALiPPSO, but it might be useful at some point to monitor the resources consumed during the execution of the main function
memory_usage() = round(parse(Int, split(read(`ps -p $(getpid()) -o rss`, String))[2]) / 1024^2, digits=3) 

"Print the state of CALiPPSO maximum displacement and inflation factor when maximum number of iterations is exceeded"
function print_failed_max_iters(max_iters::Int64, sqrΓ, tol_conv_Γ, max_Si, tol_S_conv; color::Symbol=:red)
    printstyled("Maximum number of iterations (=$max_iters) reached and failed to produce a jammed packing!!\n", bold=true, color=color)
    println("Convergence criterion for √Γ-1: ", tol_conv_Γ, "\t and for max s_i: ", tol_S_conv)
    printstyled("Final values before stopping: \t √Γ-1 = ", sqrΓ-1, ", \t max s_i = ", max_Si, "\n\n", bold=true)
end


"Print information about CALiPPSO progress and current status of the system whenever a non-isostatic configuration is created."
function print_non_isostatic(d::Int64, zs::Vector{Int64}, non_rattlers::Vector{Int64}, t::Int64, max_Si::T, bound_s::T, f_mismatch::T, small_fs::Vector{T}) where T<:AbstractFloat
    
    println("Iteration: ", t)
    printstyled("Non-isostatic packing!!!\n", color=:yellow, bold=true)
    
    Nc = 0.5*sum(zs[non_rattlers]); 
    N = length(zs); Nnr = length(non_rattlers)
    rattlers = setdiff(1:N, non_rattlers); Ncr = 0.5*sum(zs[rattlers])
    zs_rats = zs[rattlers]; inds_max_zs = sortperm(zs_rats, rev=true)[1:min(10, length(zs_rats))]
    iso_gap = d*(Nnr-1)+1 - Nc
    
    println("N_non-rattlers = ", Nnr, "\t N_c [stable particles] = ", Nc, "\t Isostaticity gap = ", iso_gap, "\t N_c [rattlers] = ", Ncr)
    println("Amount of rattlers = ", length(rattlers))
    println("(rattlers, z_rattlers) with contacts = ", [rattlers[inds_max_zs] zs_rats[inds_max_zs]])
    
    printstyled("max |sᵢ| (stable particles) = ", max_Si, "; \t bound on |sᵢ| = ", bound_s, "\t Difference of these quantities = ", max_Si - bound_s , "\n", bold=true)
    printstyled("----------------------------------\n", color=:yellow)
    println("Force mismatch = ", f_mismatch)
    println("Smallest forces = ", small_fs)
    
    printstyled("__________________________________________________________________________\n\n", color=:yellow)
end


function print_non_isostatic(d::Int64, zs::Vector{Int64}, non_rattlers::Vector{Int64})
    
    printstyled("Non-isostatic packing!!!\n", color=:yellow, bold=true)
    
    Nc = 0.5*sum(zs[non_rattlers])
    N = length(zs); Nnr = length(non_rattlers)
    rattlers = setdiff(1:N, non_rattlers); Ncr = 0.5*sum(zs[rattlers])
    zs_rats = zs[rattlers]; inds_max_zs = sortperm(zs_rats, rev=true)[1:min(10, length(zs_rats))]
    iso_gap = d*(Nnr-1)+1 - Nc
    
    println("N_non-rattlers = ", Nnr, "\t N_c [stable particles] = ", Nc, "\t Isostaticity gap = ", iso_gap, "\t N_c [rattlers] = ", Ncr)
    println("Amount of rattlers = ", length(rattlers))
    println("(rattlers, z_rattlers) with contacts = ", [rattlers[inds_max_zs] zs_rats[inds_max_zs]])
    printstyled("__________________________________________________________________________\n\n", color=:yellow)
end


"Print info about CALiPPSO progress, including new values of ``|s_i|``, density and ``R``, statistics of constraints and contacts, forces, etc."
function print_monitor_progress(max_Si::T, R::T, N::I, L::T, d::I, zs::Vector{I}, non_rattlers::Vector{I}, n_constr::Vector{I}, iso_cond::Bool, Γ_mismatch::T, f_mismatch::T, small_fs::Vector{T};
    dig_R::I=9, dig_z::I = Int(ceil(log10(N*d))), dig_S::I=13, color::Symbol=:cyan) where {I<:Int, T<:AbstractFloat}
     
     ϕ = packing_fraction(d, R, N, L)
     
     Nnr = length(non_rattlers)
     rattlers = setdiff(1:N, non_rattlers)
     zs_rat = zs[rattlers]

     # Statistics number of constraints
     max_constr = maximum(n_constr); avg_constr = round(mean(n_constr), digits=2); std_constr = round(std(n_constr), digits=2)

     # Parameters for printing output
     max_Si_print = round(max_Si, digits = dig_S)
     
    printstyled("\t Max |sᵢ| (stable particles) = ", max_Si_print ,"\n New  values of (φ, R) = ", round.([ϕ, R],  digits=dig_R), "\n", bold=true)

     if Nnr>0
        zs_non_rat = zs[non_rattlers]
        
        if iso_cond
            println("(Max, mean±std) constraints per particle per d: ", [max_constr, avg_constr, std_constr], "\t Isostatic: ", iso_cond, "\t (Mean, Max, total in rattlers) z = ", [round(mean(zs_non_rat), digits=dig_z), maximum(zs), sum(zs_rat) ] , "\t\tNon_rattlers= ", Nnr)
        else
            iso_gap = (d*(Nnr-1)+1) - 0.5*sum(zs_non_rat) 
            println("(Max, mean±std) constraints per particle per d: ", [max_constr, avg_constr, std_constr], "\t Isostatic: ", iso_cond, ";  isostaticity gap = ", iso_gap, "\t (Mean, Max, total in rattlers) z = ", [round(mean(zs_non_rat), digits=dig_z), maximum(zs), sum(zs_rat) ] , "\t\tNon-rattlers= ", Nnr)
        end

     else
        printstyled("All particles are rattlers!!\n", color=:yellow)
        println("Max constraints per particle per d: ", max_constr, "\t Isostatic: NA", "\tNon-rattlers= ", 0, "\tTotal z of rattlers = ", sum(zs_rat))
     end

    printstyled("----------------------------------\n", color=color)
    println("Mismatch in ∂ℒ/∂Γ = 0 constraint = ", Γ_mismatch)
    println("Force mismatch = ", f_mismatch)
    println("Sample of smallest forces = \t", small_fs)
    
    printstyled("__________________________________________________________________________\n\n", color=color)
end

"Print info about CALiPPSO progress, including new values of ``|s_i|``, density and some statistics of ``R``, statistics of constraints and contacts, forces, etc."
function print_monitor_progress(max_Si::T, Rs::Vector{T}, L::T, d::I, zs::Vector{I}, non_rattlers::Vector{I}, n_constr::Vector{I}, iso_cond::Bool, f_mismatch::T, small_fs::Vector{T};
     dig_R::I=9, dig_z::I = Int(ceil(log10(length(Rs)*d))), dig_S::I=13, color::Symbol=:cyan) where {I<:Int, T<:AbstractFloat}
     
     ϕ = packing_fraction(d, Rs, L)
     N = length(Rs)
     
     Nnr = length(non_rattlers)
     rattlers = setdiff(1:N, non_rattlers)
     zs_rat = zs[rattlers]

     # Statistics number of constraints
     max_constr = maximum(n_constr); avg_constr = round(mean(n_constr), digits=2); std_constr = round(std(n_constr), digits=2)

     # Parameters for printing output
     max_Si_print = round(max_Si, digits = dig_S)
     
     Rmin = minimum(Rs); Ravg = mean(Rs); Rmax = maximum(Rs)
     printstyled("\t Max |sᵢ| (stable particles) = ", max_Si_print ,"\n New  values of φ = ", round(ϕ, digits=dig_R)," and (Rₘᵢₙ, ⟨R⟩, Rₘₐₓ) = ", round.([Rmin, Ravg, Rmax],  digits=dig_R), "\n", bold=true)

     if Nnr>0
        zs_non_rat = zs[non_rattlers]
        
        if iso_cond
            println("(Max, mean±std) constraints per particle per d: ", [max_constr, avg_constr, std_constr], "\t Isostatic: ", iso_cond, "\t (Mean, Max, total in rattlers) z = ", [round(mean(zs_non_rat), digits=dig_z), maximum(zs), sum(zs_rat) ] , "\t\tNon_rattlers= ", Nnr)
        else
            iso_gap = (d*(Nnr-1)+1) - 0.5*sum(zs_non_rat) 
            println("(Max, mean±std) constraints per particle per d: ", [max_constr, avg_constr, std_constr], "\t Isostatic: ", iso_cond, ";  isostaticity gap = ", iso_gap, "\t (Mean, Max, total in rattlers) z = ", [round(mean(zs_non_rat), digits=dig_z), maximum(zs), sum(zs_rat) ] , "\t\tNon-rattlers= ", Nnr)
        end

     else
        printstyled("All particles are rattlers!!\n", color=:yellow)
        println("Max constraints per particle per d: ", max_constr, "\t Isostatic: NA", "\tNon-rattlers= ", 0, "\tTotal z of rattlers = ", sum(zs_rat))
     end

    printstyled("----------------------------------\n", color=color)
    println("Force mismatch = ", f_mismatch)
    println("Sample of smallest forces = \t", small_fs)
    
    printstyled("__________________________________________________________________________\n\n", color=color)
end

"Print that CALiPPSO has converged, and some basic info (final ϕ, number of non rattlers, etc.)."
function print_converged(t::I, status, sqrΓ::T, max_Si::T, R::T, N::I, L::T, d::I, n_constr::Vector{I}, Nnr::I; dig_R::I=8, dig_S::I=15, color::Symbol=:green) where {I<:Int, T<:AbstractFloat}
    ϕ = packing_fraction(d, R, N, L)

     # Statistics number of constraints
     max_constr = maximum(n_constr); avg_constr = round(mean(n_constr), digits=2); std_constr = round(std(n_constr), digits=2)

    printstyled("\n\n\tCALiPPSO converged! \t \t Status of last LP optimization: ", status, "\n", bold=true, color=color)
    printstyled("Iterations to convergence = ", t, ",\t√Γ-1 = ", sqrΓ-1, ",\t Max displacement = ", round(max_Si, digits=dig_S), ",\t (ϕ, R) = ",  round.([ϕ, R], digits=dig_R), "\n", bold=true, color=color)
    println("\tNon_rattlers= ", Nnr, "\t% of rattlers = ", round(100*(1-Nnr/N), digits=3),  "\t(Max, mean±std) constraints per particle per d: ", [max_constr, avg_constr, std_constr])
end

function print_converged(t::I, status, sqrΓ::T, max_Si::T, Rs::Vector{T}, N::I, L::T, d::I, n_constr::Vector{I}, Nnr::I; dig_R::I=8, dig_S::I=15, color::Symbol=:green) where {I<:Int, T<:AbstractFloat}
    ϕ = packing_fraction(d, Rs, L)
    # Statistics of radii
    Rmin = minimum(Rs); Ravg = mean(Rs); Rmax = maximum(Rs)

     # Statistics number of constraints
     max_constr = maximum(n_constr); avg_constr = round(mean(n_constr), digits=2); std_constr = round(std(n_constr), digits=2)

    printstyled("\n\n\tCALiPPSO converged! \t \t Status of last LP optimization: ", status, "\n", bold=true, color=color)
    printstyled("Iterations to convergence = ", t, ",\t√Γ-1 = ", sqrΓ-1, ",\t Max displacement = ", round(max_Si, digits=dig_S), ",\t Final  values of φ = ", round(ϕ, digits=dig_R), " and (Rₘᵢₙ, ⟨R⟩, Rₘₐₓ) = ", round.([Rmin, Ravg, Rmax],  digits=dig_R), "\n", bold=true, color=color)
    println("\tNon_rattlers= ", Nnr, "\t% of rattlers = ", round(100*(1-Nnr/N), digits=3),  "\t(Max, mean±std) constraints per particle per d: ", [max_constr, avg_constr, std_constr])
end


"""
    print_info_convergence(final_packing::MonoPacking, isostatic::Bool, time, memory; digs::Int64=2)

Print extra details about the configuration obtained once CALiPPSO converges.

The maximum mismatch in the force balance condition is computed (and printed), and also the 
isostaticity of the packing is assessed. Besides, some information about performance 
(execution time and allocated memory) is also printed out.

See also [`print_converged`](@ref), [`print_monitor_progress`](@ref).
"""
function print_info_convergence(final_packing::AbstractPacking{d, T}, isostatic::Bool, time, memory; digs::Int64=2) where {d, T<:AbstractFloat}
    
    isostatic ? (col_p = :green) : (col_p = :magenta)

    N = length(final_packing)
    non_rattlers = get_non_rattlers(final_packing); Nnr = length(non_rattlers)
    rattlers = setdiff(1:N, non_rattlers)
    zs = get_coordination_number(final_packing, only_stable=false)

    force_mismatch = maximum(norm.(total_force(final_packing)))
    minutes = round(time/60, digits=digs)
    memory_r = round(memory, digits=digs)

    if isostatic
        printstyled("Isostaticity achieved: ", isostatic, "\t Non-rattlers = ", Nnr, "\t% of rattlers = ", round(100*(1-Nnr/N), digits=3), "\tN_contacts [in stable particles] = ", 0.5*sum(zs[non_rattlers]),"\t z in rattlers = ", sum(zs[rattlers]), "\n", color=col_p)
        printstyled("Maximum force equilibrium mismatch = ", force_mismatch, "\n", color=col_p)
    else
        iso_gap = (d*(Nnr-1)+1) - 0.5*sum(zs)
        printstyled("Isostaticity achieved: ", isostatic, "  ; isostaticity gap = ", iso_gap, "\t Non-rattlers = ", Nnr, "\t% of rattlers = ", round(100*(1-Nnr/N), digits=3), "\tN_contacts [all particles] = ", 0.5*sum(zs),"\t z in rattlers = ", sum(zs[rattlers]), "\n", color=col_p)
        printstyled("Maximum force equilibrium mismatch = ", force_mismatch, "\n", color=col_p)
    end
    println("Time to finish = ", minutes, " minutes;\t Memory allocated (GB): ",  memory_r)
end
