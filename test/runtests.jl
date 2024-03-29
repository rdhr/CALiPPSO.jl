using CALiPPSO
using Test

# precompile_main_function()

Lt = 4.0;
Nt = 30;
rt = 0.52;
dt = 3;
tol_Γ = 1e-3;
tol_S = 1e-2;
tol_f = 10 * tol_S;

optimizer = CALiPPSO.GLPK.Optimizer(want_infeasibility_certificates=false)


Xs_comp = Matrix(transpose([
    1.3022036173223075 3.100681668734574 2.650145868235529
    1.1155120883337997 2.414092386379111 0.3717676414091988
    0.5930622140684934 0.42734603856532427 1.3262407197673065
    3.571441842561353 2.957136483108787 2.6002633302062668
    2.3529766061581237 0.48997057275911615 3.417191033324789
    1.7112056139610727 0.8347268069230269 2.229786164196563
    3.6817595377155348 0.6725943671221577 2.4334027050690317
    1.9126951936859324 1.5811607905044598 3.460749757890367
    2.10364317266594 2.310954222951657 2.6450251464477823
    3.2140912406591706 1.9574663303701785 3.5356829974762487
    1.5789002308497055 2.218608327103566 1.664307104091713
    0.720581832264064 3.921707081207914 3.2620078999221587
    2.459219015775952 3.503789658633769 0.9321778468174386
    3.6785107203364795 3.337897259842034 0.6462417289544637
    3.4240729350390406 2.2304592814837587 1.2859711802875369
    2.200649934148873 2.428279519213076 0.38501446881727297
    0.35487777941165266 2.5227027830749185 1.7324205587874006
    0.6937597057720284 1.675212825207569 1.1950801557221888
    2.5783261715853776 0.49417638526472807 1.2809483317712305
    0.25868616273463285 2.114656556778108 3.8413865080579885
    2.269006112530807 3.5797496686786596 2.8118840862720855
    0.4016449295907414 2.092677546492193 2.715279086250721
    0.936166141399041 0.19318442951467052 0.3547261200994525
    3.4700487156587823 1.163132895735302 1.0834877206251514
    3.307514775505691 3.9075523242046204 3.036230135532204
    1.1350771834782938 1.2830440239210912 0.11294699337022074
    1.7035741859408704 3.868891311445325 1.722953156549603
    3.896005446471116 3.1189426973423595 3.607500484893032320
    3.653047120545878 0.25044725175167404 3.9730036860708076
    0.6334501855938681 1.1831285025320382 3.1868918476405756]))
cen_comp = PeriodicVectors(Xs_comp, Lt);

@testset "Periodic Number tests" begin
    @test abs(PeriodicNumber(2.5, 2.0).value - 0.5) <= eps()
    @test abs(PeriodicNumber(1.3).value - 0.3)<=eps()
    @test CALiPPSO.value(cen_comp) == Xs_comp
    @test eltype(cen_comp) <: CALiPPSO.SVector
    @test length(cen_comp) == Nt
    @test all(≥(0), CALiPPSO.distances_between_centers(cen_comp))
end



Jpack_comp, conv_info_comp, Γs_comp, smax_comp, isos_comp = produce_jammed_configuration!(cen_comp, rt, ℓ0=Lt, initial_monitor=0, tol_Γ_convergence=tol_Γ, tol_S_convergence=tol_S, tol_mechanical_equilibrium=tol_f, verbose=false, add_bridges=true, optimizer=optimizer)



@testset "Packing compilation configuration (monodisperse method)" begin
    @test conv_info_comp.converged
    @test Jpack_comp.jammed
    @test Jpack_comp.mechanical_equilibrium
    @test length(Jpack_comp.Particles) == Nt
    @test Γs_comp[end] - 1 <= tol_Γ
    @test smax_comp[end] <= tol_S
    @test Jpack_comp.R ≥ rt
end

cen_comp = PeriodicVectors(Xs_comp, Lt); Rs_comp = rt*ones(Nt)
Jpack_comp, conv_info_comp, Γs_comp, smax_comp, isos_comp =  produce_jammed_configuration!(cen_comp, Rs_comp, ℓ0=Lt, initial_monitor=0, tol_Γ_convergence=tol_Γ, tol_S_convergence=tol_S, tol_mechanical_equilibrium=tol_f, verbose=false, optimizer=optimizer)

@testset "Packing compilation configuration (polydisperse method)" begin
    @test conv_info_comp.converged
    @test Jpack_comp.jammed
    @test Jpack_comp.mechanical_equilibrium
    @test length(Jpack_comp.Particles) == Nt
    @test Γs_comp[end] - 1 <= tol_Γ
    @test smax_comp[end] <= tol_S
    @test CALiPPSO.get_radii(Jpack_comp) == Rs_comp
end