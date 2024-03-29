# CALiPPSO: A Linear Programming Algorithm for Jamming Hard Spheres

The CALiPPSO algorithm was introduced in [our paper](https://arxiv.org/abs/2203.05654) (we = *Artiaco, Díaz, Parisi, and Ricci-Tersenghi*). A [Julia](https://julialang.org/) implementation of this algorithm is available through the [CALiPPSO.jl package](https://github.com/rdhr/CALiPPSO.jl), developed by ourselves. This is the documentation of such Julia package. 

As we explain in our paper, CALiPPSO is an iterative Linear Programming algorithm to produce jammed packings of hard-spheres (HS) *in arbitrary* dimensions, $d$. Besides, it works for both *mono*- and *poly*disperse packings, as shown below

```@raw html
<p align="center" margin=0px>
  <img src="https://www.dropbox.com/s/u4718ftwoz1f99z/3dRCP.jpg?raw=1" alt="Monodisperse 3d" width=48%>
  <img src="https://www.dropbox.com/s/1s2t7so2nyih502/isostatic-2d-polydisperse.png?raw=1" alt="Polydisperse 2d" width=46%>

  <small> (Left: Monodisperse packing of 16k particles; coloured according to their number of contacts. Right: Polydisperse packing of 1024 disks, with radii from a log-normal distribution, and network of contacts drawn.</small>
</p>
```


Essentially, our method consists on a **C**hain of **A**pproximate **Li**near **P**rogramming for **P**acking **S**pherical **O**bjects.

The code contained in the [CALiPPSO.jl package](https://github.com/rdhr/CALiPPSO.jl) is written in pure [Julia](https://julialang.org/), and makes extensive use of  the (wonderful) optimization and modelling package [JuMP.jl](https://github.com/jump-dev/JuMP.jl), as well as other existing packages ([StaticArrays.jl](https://juliaarrays.github.io/StaticArrays.jl/stable/), [GLPK.jl](https://github.com/jump-dev/GLPK.jl), etc.). 

This package is licensed under the MIT license, so please feel free to use/modify/improve this code as better suits you. We only ask you to cite our work if you find it useful.

```
@article{artiacoHardsphereJammingLens2022,
  title = {Hard-Sphere Jamming through the Lens of Linear Optimization},
  author = {Artiaco, Claudia and D{\'i}az Hern{\'a}ndez Rojas, Rafael and Parisi, Giorgio and {Ricci-Tersenghi}, Federico},
  year = {2022},
  month = nov,
  journal = {Physical Review E},
  volume = {106},
  number = {5},
  pages = {055310},
  publisher = {{American Physical Society}},
  doi = {10.1103/PhysRevE.106.055310},
  url = {https://link.aps.org/doi/10.1103/PhysRevE.106.055310}
}
```

If you already know how CALiPPSO works (or know the theory behind our paper), you can skip to [Installation](@ref) or [Basic usage](@ref) sections (if you have already installed the package). Otherwise, the [Theory behind CALiPPSO](@ref) section provides the essentials to understand our method.


## Before reading this documentation

In the next sections we describe in some detail our implementation of the CALiPPSO algorithm. But before reading how our code works we suggest that if you're new to [Julia](https://julialang.org/) it might be useful to have [its documentation](https://docs.julialang.org/en/v1/) at hand, specially for understanding some terminology (although we tried to keep it to a minimum).
Besides, several of the functions we define rely on [JuMP](https://jump.dev/)'s functionality, so if you are unfamiliar with this package consider skimming through [its documentation](https://jump.dev/JuMP.jl/stable/) (specially the parts on [creating a model](https://jump.dev/JuMP.jl/stable/manual/models/#Create-a-model) and defining/accessing [constraints](https://jump.dev/JuMP.jl/stable/manual/constraints/)). 

## Contents

```@contents
Pages = ["theory.md", "installation.md", "basic_usage.md", "mainfunction.md", "changing_default.md", "tests.md", "types.md", "issues.md", "api.md"]
Depth = 4
```



