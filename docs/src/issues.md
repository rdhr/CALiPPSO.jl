# Problem solving

SECTION UNDER CONSTRUCTION.

## Possible issues 

1. Modify `sbound` if isostaticity gap is negative for several steps due to `abs(s_i)=sbound*R`
2. If isostaticity gap is too big (and positive), possibly there is small precision in the value of dual variables (*e.g.* when using Tulip). This could be solved by looking for dual variable larger than a given threshold value, but still needs to be implemented.
   1. This might be the case even when such gap is not too big. For instance, using GLPK in the LS tests, the initial iterations throw a gap of about +6 (due to very small forces)
3. A *negative* isostaticity can be caused by the maximal displacement (of *stable* particles) being equal to the bound on $|s_{i,\mu}^\star|$. If this happens, try using a larger value of `sbound`. 
   1. In reality, what happens is that there is an *active* dual variable related to such bound on the displacements and *not* to a contacts. Thus, the system is nonetheless isostatic, but a contact has not being properly counted, since it is not included in the lists of `constraints` in the main loop.
4. The right convergence tolerance values, *i.e.* the right choice of `tol_Γ_convergence` and `tol_S_convergence` possibly depend on the solver
5. `Tulip.jl` give an isostatic config only when the status of the optimization is `ITERATION_LIMIT`, which also entails that forces are determined with very, very low accuracy. Notably, even though this happens, the force mismatch is very small.

### Setting the required precision

## Overlaps occur

### Due to numerical error

### REAL overlaps

When this occurs, it is likely that the value of `sbound` is too large with respect of the value of ``\ell`` used. So try using a larger value of `ℓ0` (if this occurred when ``\Gamma^\star`` was still relatively large).
   -  Try using a smaller value of `sbound`.
   -  Or try redefining the bounds with which ``\ell`` is adjusted as CALiPPSO makes progress.
   -  These three options can be set as keywords arguments of [`produce_jammed_configuration!`](@ref) as explained [in this section](@ref kwargs-control). 
   -  See the docstring of [`bounds_and_cutoff`](@ref CALiPPSO.bounds_and_cutoff) for more info.

## Convergence not attained

