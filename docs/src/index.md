# CALiPPSO: A Linear Programming Algorithm for Jamming Hard Spheres


This is the accompanying code to the paper XXXX, by *Artiaco, Díaz, Parisi, and Ricci-Tersenghi*. There, we introduce and describe an iterative Linear Programming algorithm to produce jammed packings of hard-spheres (HS) *in arbitrary* dimensions, $d$.
Essentially, our method consists on a **C**hain of **A**pproximate **Li**near **P**rogramming for **P**acking **S**pherical **O**bjects.
The code contained here is written in pure [Julia](https://julialang.org/), and makes extensive use of  the (wonderful) optimization and modelling package [JuMP.jl](https://github.com/jump-dev/JuMP.jl). 

Feel free to use/modify/improve this code as better suits you. We only ask you to give credit to our work.

```
@article{...}
```
If you already know how CALiPPSO works (or know the theory behind our paper), you can safely skip to [Basic usage](@ref) section, where we describe our code in some detail. Otherwise, the [Theory behind CALiPPSO](@ref) section provides the essentials to understand our method.

## Disclaimer

## Contents

```@contents
Pages = ["theory.md", "usage.md", "basic_usage.md", "mainfunction.md", "changing_default.md", "tests.md", "types.md", "issues.md", "todos.md", "api.md"]
Depth = 4
```



