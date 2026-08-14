# AUTO-GENERATED one-shot fixture — DO NOT EDIT.
# Baked ORIGINAL DG-Teulolsky (/data/haiyangw/claude/DG-Teulolsky @ 0a88967, its own env, julia 1.12.6)
# by scratchpad/capture_td_dg.jl on x86_64-linux-gnu. (a) Mode1D Schwarzschild ℓ=0 Price tail: ψ probe samples at r*=500 and scri + fit_tail exponents.
# EXACT CONFIG: Mode1DProblem(0; N=8, K=120, ρL=-200.0, R=150.0, s=1200.0, P=4, M=1, extract_rstar=[500.0, Inf]); solve_mode1d(tfinal=6000.0, dt=0.0866, record_every=100, momentum_ic=(amp=1/sqrt(2π*10^2), ρ0=30.0, w=10.0)); fit_tail window=(2500,6000)
const TD_MODE1D = (
    nsamp = 693,
    idx = [14, 69, 173, 346, 520, 693],
    tau = [112.58000000000148, 588.8799999999155, 1489.520000000275, 2987.7000000018234, 4494.540000001288, 5992.719999994969],
    psi_r500 = [-4.3494119340596235e-28, 0.013942226440503948, 6.888344221376623e-5, 7.151045842978891e-6, 2.0370656284914208e-6, 8.504845472262575e-7],
    psi_scri = [1.3205047100408938e-75, 5.733658169117945e-34, 0.0011044591759450844, 2.8719118004034496e-5, 8.474233375029043e-6, 4.008475059437962e-6],
    n_r500 = -3.067068266658501,
    n_scri = -2.9090262084732705,
    # problem-construction probes (grounds inverse_tortoise→Spacetime + verbatim C):
    probe_rho = [494.2246204961792, 1200.0],
    probe_k = [60, 120],
    C_3_17 = 0.0002982994266908252,
    C_5_60 = 1.7132745945921235e-8,
    C_scri = 0.0,
    H_5_100 = 0.821780794010727,
)
