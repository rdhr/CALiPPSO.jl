# Library

## Types

```@docs
PeriodicNumber
MonoParticle{d,T}
MonoPacking
convergence_info
```

## Functions and Methods

### Main Exported Functions


```@docs
produce_jammed_configuration
network_of_contacts
check_for_overlaps
packing_fraction
is_isostatic
get_coordination_number
total_force
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
```

### Secondary functions

```@docs
CALiPPSO.solve_LP_instance
CALiPPSO.fine_tune_forces!

```

### Functions for printing CALiPSSO info and progress

```@docs
CALiPPSO.print_monitor_progress
CALiPPSO.print_converged
```

## Other functions

```

```@autodocs
Modules = [CALiPPSO]
Order = [:function]
```

```