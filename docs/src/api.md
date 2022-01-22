# Library

## Types

```@docs
PeriodicNumber
MonoParticle{d,T}
MonoPacking{d,T}
convergence_info
```

## Functions and Methods

### Main Exported Functions


```@docs
produce_jammed_configuration
network_of_contacts
check_for_overlaps
check_for_overlaps(::MonoPacking)
generate_random_configuration
total_force
packing_fraction
is_isostatic
get_coordination_number
get_non_rattlers
get_rattlers
```

### Other exported functions

```@docs
PeriodicVectors
volume_d_ball
norm(v::AbstractVector{<:PeriodicNumber})
```

### Constructors

```@docs
PeriodicVector
MonoPacking()
MonoPacking(::Vector, ::Vector, ::Vector, ::Real)
```

### Secondary functions

```@docs
CALiPPSO.solve_LP_instance
CALiPPSO.fine_tune_forces!
CALiPPSO.update_packing_forces!
CALiPPSO.add_non_overlapping_constraints!
CALiPPSO.bounds_and_cutoff
CALiPPSO.obtain_non_rattlers
CALiPPSO.MIC_vector
CALiPPSO.MIC_distance
```

### Functions for printing CALiPSSO info and progress

```@docs
CALiPPSO.print_monitor_progress
CALiPPSO.print_converged
CALiPPSO.print_info_convergence
CALiPPSO.print_failed_max_iters
CALiPPSO.print_non_isostatic
```

## Other functions

```

```@autodocs
Modules = [CALiPPSO]
Order = [:function]
```

```