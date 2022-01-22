5. If the `verbose` option is set to `true`, the following information is printed (see [the dedicated section](@ref output-process) on how to interpret such info):
   1. The number of LP iteration and the state of the optimization thrown after `optimize!(LP_model)` finishes, *i.e.* the value of `termination_status(LP_model)`. This latter could be, *e.g.* `OPTIMAL`, `INFEASIBLE`, `TIME_LIMIT` , etc.
   2. The value of ``\sqrt{\Gamma^\star} -1`` and ``\max |s^\star_{i,\mu}|``, where ``i \in \text{non-rattlers}``.
   3. The time required to solve the LP problem.
   4. Few statistics about the number of constraints included in `LP_model`
   5. Information about whether the preliminary configuration is isostatic or not, as well as the isostaticity gap, ``N_c - N_{dof}``, in the latter case.
   6. Few statistics of the coordination number.
   7. Number of non-rattlers (or stable particles).
   8. Maximum mismatch in the force balance (per particle).
   9. Sample of the 10 smallest forces in the system. This information is also useful for evaluating if numerical issues are present.







defined few *composite types* or `struct`'s (aka *objects* in other languages). They are:

- `PeriodicNumber`
- two particle types (either `MonoParticle` or `Particle`)
- two packing types: `MonoPacking` (composed of `MonoParticle`'s) and `PolyPacking` (composed `Particle`'s), and used to model mono- or polydisperse systems, respectively.

The specific properties of each of these types are explained [below](#types-defined-in-this-package), but the idea is that each of them contains the necessary information of the object it represents. For instance:
- A `Particle` is assigned a position (as an `SVector` of `PeriodicNumber` elements), a radius, a set of neighbours, and the corresponding set of contact vectors and forces.
  - `MonoParticle` is essentially the same type of object, except that no radius is assigned to it (since it is stored as a field of the packing itself.)
- Analogously, a `PolyPacking` contains an array of `Particle` objects, and also include information about its isostaticity, mechanical stability, and whether it is jammed or not.
  - Similarly, `MonoPacking` has the same fields and an extra one to store the size of all `MonoParticles` it is composed of.