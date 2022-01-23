# Introduction: how CALiPPSO works and some terminology 

## Theory behind CALiPPSO

As we explain in the article, our approach works by mapping the exact (non-convex) problem of jamming hard spheres (HS), into a series of linear optimization problems (and thus readily solvable). To describe the basic idea of our method, let us assume that we are given a system of ``N`` (hyper-)spheres, with diameters ``\vec{\sigma}=\{\sigma_i\}_{i=1}^N`` and centres at positions ``\vec{\mathbf{r}}=\{\mathbf{r}_i\}_{i=1}^N``. (We use ``\mathbf{x}=(x_1, \dots, x_d)`` to denote a vector in ``d`` dimensions, and ``\vec{\bullet}=\{\bullet_i\}_{1=1}^N`` to write the set values of ``\bullet`` from every particle.) With this information, we formulate an instance of a *constrained* linear optimization problem, whose purpose is to maximise the system's density, but avoiding any overlap between the particles. More precisely, we aim to find the optimal displacements (``\vec{\mathbf{s}}^\star=\{\mathbf{s}_i^\star\}_{i=1}^N``) that allow to maximize an inflation factor of their diameter (``\Gamma``), without producing any overlap. 
Explicitly, for given ``\vec{\mathbf{r}}`` and ``\vec{\sigma}`` (*i.e.* these variables play the role of parameters), the exact optimization problem reads
```math
\tag{1} \begin{aligned}
\max & \quad \Gamma \\
  \text{s. t.} \qquad |\mathbf{r}_{i} + \mathbf{s}_i - (\mathbf{r}_j + \mathbf{s}_j)|^2  & \geq \Gamma \sigma_{ij}^2  \qquad \forall \ 1\leq i < j \leq N\, .
\end{aligned}
```
Instead, if second order terms of ``|\mathbf{s}|_i`` are neglected (for instance, if the initial configuration is already close to jamming), we can linearized the previous optimization problem and turn it into a Linear Programming (LP) one. Thus, with our algorithm we solve a LP instance that reads:
```math
\begin{aligned} \tag{2}
  \max & \quad \Gamma \\
  \text{s. t.} \qquad |\mathbf{r}_{ij}|^2 + 2 \mathbf{r}_{ij}\cdot \mathbf{s}_{ij}   & \geq \Gamma \sigma_{ij}^2  \qquad \forall \ 1\leq i < j \leq N\, ;
\end{aligned}
```
where ``\mathbf{r}_{ij} = \mathbf{r}_i - \mathbf{r}_j``, and analogously for ``\mathbf{s}_{ij}``, while ``\sigma_{ij}=(\sigma_i+\sigma_j)/2`` is the sum the pair radii. Thus, we look for the optimal values of ``(\vec{\mathbf{s}}^\star, \Gamma^\star)`` ─the design variables of our problem.

Once the LP optimization has been carried out, the positions and sizes are updated with these optimal values, following the rule ``(\vec{\mathbf{r}} \to \vec{\mathbf{r}} + \vec{\mathbf{s}^\star}, \sigma_i \to \sqrt{\Gamma^\star} \sigma_i)``, whence a new LP instance is created and solved. (The usage of the ``\sqrt{\Gamma}`` instead of ``\Gamma`` is a minor detail used to keep both the objective and the constraints linear in the design variables.) 
Updating a configuration using ``\Gamma^\star`` and ``\vec{\mathbf{s}}^\star`` and formulating a new LP instance leads to an iterative process that clearly approaches the jamming point, because at each step the density is increased without occurring into any overlaps.
This process continues until ``(\vec{\mathbf{s}}^\star, \Gamma^\star) = (\vec{\mathbf{0}}, 1)``, which we term *convergence condition*. For brevity, we call **Iterative Linear Programming (ILP)** such process of transforming a configuration using the chain of optimal solutions to the LP problems. Besides, we refer to one LP optimization as *one iteration* of the CALiPPSO algorithm.



It should be mentioned that we only deal with the *linearized* version of the exact non-overlapping constraints (in fact, the exact constraints are what causes the non-convexity of the original problem, Eqs. (1)). This is why our method is amenable to standard numerics, but also implies that several LP optimizations are required because, before convergence, particles are not really in contact. Nevertheless, it is rather intuitive that once CALiPPSO converges it produces a jammed packing, since particle cannot be displaced nor the system's density further increased (this is what the convergence conditions actually encode). A more detailed proof is given in our paper where we show that once the convergence criterion is met, the packing thus produced is a well defined jammed state. 
This means that:
1. Mechanical equilibrium holds for each particle;
2. The number of contacts, ``N_c``, exactly matches the *total* number of degrees of freedom (dof), ``N_{dof}``. This guarantees the global stability of the packing.
   - For systems with periodic boundary conditions, as we consider here, we have ``N_{dof}=d*(N-1)+1``, due to the ``d`` possible uniform translations of the system. The extra dof comes from the fact that ``\Gamma`` (or, equivalently, the density) is also a variable in our setup.
   - We henceforth refer to the condition of ``N_c=N_{dof}`` as *isostaticity* for short[^iso].
   - We mention in passing that the extra contact with respect fo the configurational degrees of freedom is very important because the onset of the critical properties of jamming are marked precisely by such additional contact. 


All in all, the idea behind CALiPPSO is that if the initial configuration is already close to its jamming point, then the error made by using a linear approximation of the exact optimization problem is not going to be very large. Importantly, the feasible set of the linearized problem (*i.e.* the set of solutions satisfying the constraints) is always contained in the feasible set of the original, non-convex problem. This implies that any solution to the LP problem of Eqs.(2), is also a valid configuration considering the *exact* non-overlapping constraints, Eqs. (1). However, as we describe [below](@ref The-initial-conditions), CALiPPSO is robust enough to produce jammed packings even if the initial configuration is *far* from jamming.


## Contact forces

Obtaining the contact forces in a HS configuration is in general not a trivial task, because the HS interaction potential is (*extremely*) singular: it is ``0`` if particles do not overlap, and ``\infty`` as soon as they do. The interaction energy is consequently always zero (in any valid configuration), and therefore forces cannot be obtained by taking the derivatives of a potential. In turn, they have to be determined by solving the force balance equations of the whole system. Fortunately, in our case the Karush─Kuhn─Tucker conditions are equivalent to this the same set of equations. 

The full proof is given in the paper, but the basic idea is that the *active* dual variables (or Lagrange multipliers in the physics parlor) associated to the linearized non-overlapping constraints, satisfy the same equations as the forces in an stable packing. Therefore, such dual variables must be the contact forces; see Eqs. (8) and (10) of our paper. (In some sense, this is just an application of the Lagrange equations for computing the generalized forces associated to some mechanical constraint.)

Recall that *active* dual variables are associated to constraints that are *saturated*, and therefore, once convergence ensues, they are always associated to real physical contacts.
- An interesting feature [not a bug ;)] is that the force balance condition actually holds at each CALiPPSO iteration and *not only* at jamming. That is, even if convergence has not been reached, the equations satisfied by the dual variables correspond to the ones of mechanical equilibrium. This means that after each LP optimization, the non-zero values of the dual variables are such that, given the current ``\vec{\mathbf{r}}``, they are in mechanical equilibrium, with respect to the *linear* forces (or *linearized constraints*).

In any case, the fact that there is a one-to-one mapping of forces and active dual variables solves, in a single step, the problem of unequivocally identifying real contacts and calculating the contact forces. 
- The first point is not trivial because gaps between nearly-touching particles in jammed packings have a non trivial distribution, that actually diverges as such separation ``\to 0``. Thus, it is likely to have many particles separate by very, very small distances but that, nonetheless, do not correspond to real contacts. By using the dual variables criterion, we can easily identify such contacts, usually without any problem associated to numerical accuracy.

### Rattlers

As with any other jamming algorithms, CALiPPSO also yields a (small) amount of non-stable particles. These are particles that have ``d`` contacts or less, and therefore cannot be mechanically stable, because at least one of its dof is not blocked. Such non-stable particles are usually called *rattlers* and should be excluded when analysing the stability of a packing. Therefore, the real number of degrees of freedom in our systems is ``N_{dof}=d(N_s-1)``, where ``N_s`` is the number of *stable* particles (*i.e.* those with at least ``d+1`` contacts). 

Despite the fact that rattlers are a small fraction of the total number of particles (always less than 3%), it is very important to identify them, both for properly computing the network of contacts at jamming, and for evaluating the performance of our CALiPPSO algorithm. This latter point will also be useful to understand possible numerical issues when running our ILP code, *even before convergence*. That is, from now on, we'll say that a particle is a *rattler* whenever it has at most ``d`` contacts, *independently if they are real or only linear contacts*. This will enable us to assess the isostaticity of a configuration even at intermediate steps of the [main CALiPPSO loop](@ref mainloop).

## The initial conditions

As mentioned above, CALiPPSO was developed under the assumption that the initial density, ``\varphi_0``, is already close to the jamming one, ``\varphi_J``. Yet, it is versatile enough to bring to their jamming point HS configurations even when ``\varphi_0 \simeq 0.3\varphi_J``, and the initial positions are randomly distributed (we provide some examples in [the examples section](@ref Some-examples-included)).

!!! note "POSSIBLE CAVEAT"
    Nonetheless, it should be kept in mind that when low density configurations are used as initial conditions, the packings thus obtained are somewhat atypical. That is, their value of ``\varphi_J`` is slightly, but detectably, smaller than the jamming density obtained through other methods (or by using “better” initial conditions). This indicates that it is likely that such configuration is a "high minimum" of the *free-energy landscape*. On the other hand, other important properties, such as isostaticity, critical distributions of contact forces and gaps, as well as a flat density of states, seem to be unaffected. But we have not carried out a systematic analysis in this scenario.

!!! warning "Non-isostaticity before convergence"
    When using a small value of ``\varphi_0`` it may happen that during intermediate iterations the solution of the LP instance of Eq. (2) yields a **non**-isostatic configuration, or one where mechanical equilibrium does *not* hold. This is an unusual situation, but it occurred in some of [our tests](@ref Some-examples-included) with random initial conditions, specially in ``d \geq 5`` and in ``d=2``. In high dimensions, this issue is likely caused by the fact that we included bounds on ``|\mathbf{s}_i| \ \forall 1\leq i\leq N`` and these extra constraints are not considered when performing the stability analysis. In turn, when ``d=2`` the problem is related to the fact that partial crystallization cannot be avoided, leading to an under-determined system of equations (see Sec. II.B of our paper). In any case, whenever ``d\geq 3``, the configuration obtained *once convergence is reached* should be isostatic and in mechanical equilibrium.

Now, as expected, CALiPPSO works best when the initial configuration is already highly compressed. For instance when the reduced pressure is at least ``p>100``.  (Recall that the reduced pressure is defined as ``p=\frac{\beta P}{\rho}``, with ``\beta`` the inverse temperature, ``P`` the “normal” pressure, and ``\rho=N/L^d`` the number density.) And, naturally, the higher the initial value of ``p`` the quicker CALiPPSO converges to a typical jammed state[^1]. 

How to produce a typical highly compressed configuration is itself a different problem, and we decided to use the Lubachevsky─Stillinger algorithm (with the code implementation available [here](https://cims.nyu.edu/~donev/Packing/C++/), based on [this work](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.74.041127); also accessible as [preprint](https://arxiv.org/abs/cond-mat/0608362)). But of course, you can use your favourite method. In any case, the [examples we provide](@ref Some-examples-included) where the initial condition is already close to ``\varphi_J`` were obtained with such compression method, and with ``p\geq 10^5``.


[^1]: In our paper we report results with ``p\geq 1000``, and extensive analyses have been carried out with ``p=10^7``.
[^iso]: Although, strictly speaking, isostaticity refer to the condition ``N_c=N_{dof}``. In that sense, our configuration are hyper-static, even if by a single contact. More precisely, our conditions have what is commonly called "a single state of self-stress".