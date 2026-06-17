"""
=============================================================================
 Kinetic ODE Simulation of a Branched Metabolic Pathway
 with Non-Competitive Allosteric Feedback Inhibition
=============================================================================

 Course    : Biotechnology — Systems Biology (BISB211605)
 Student   : Adelia Yusuf Ardhani
 University: Universitas Gadjah Mada — Faculty of Biology

 Description:
   This script simulates the time-course dynamics of a 4-reaction branched
   metabolic pathway using Ordinary Differential Equations (ODEs). The
   pathway converts external substrate X into a target product P through
   two intermediates A and B, while a competing branch (v4) diverts flux
   toward an unwanted byproduct.

   Product P exerts non-competitive allosteric feedback inhibition on the
   first committed step (v1), which is modeled using an extended
   Michaelis-Menten equation with an inhibition factor.

 Pathway Topology:
   X --[v1]--> A --[v2]--> B --[v3]--> P
                |                       |
               [v4]          (allosteric inhibition on v1)
                |
                v
            Byproduct

 Reactions:
   v1 : X → A   (Michaelis-Menten + non-competitive inhibition by P)
   v2 : A → B   (first-order)
   v3 : B → P   (first-order)
   v4 : A → Byproduct  (first-order, competing branch)

=============================================================================
"""

# =============================================================================
# 1. IMPORT LIBRARIES
# =============================================================================

import numpy as np                        # numerical computation
import matplotlib.pyplot as plt           # plotting
import matplotlib.gridspec as gridspec    # for multi-panel figure layout
from scipy.integrate import solve_ivp     # ODE solver
import os                                 # for creating output folder


# =============================================================================
# 2. PARAMETERS
# =============================================================================
# All concentrations are in mM, time is in hours (h), rates in mM/h.

PARAMS = {
    "V1max" : 5.0,   # Maximum velocity of reaction v1 (mM/h)
    "Km1"   : 2.0,   # Michaelis constant for v1 — concentration of X at
                     #   half-maximal rate (mM)
    "Ki"    : 3.0,   # Inhibition constant — concentration of P that reduces
                     #   v1 rate by 50% (mM)
    "X"     : 10.0,  # External substrate concentration, held constant (mM)
    "k2"    : 1.0,   # First-order rate constant for A → B (1/h)
    "k3"    : 0.8,   # First-order rate constant for B → P (1/h)
    "k4"    : 0.3,   # First-order rate constant for A → Byproduct (1/h)
}

# Simulation settings
T_START   = 0       # Start time (h)
T_END     = 48      # End time (h) — simulates a 48-hour fermentation
N_POINTS  = 1000    # Number of time points to record
Y0        = [0.0, 0.0, 0.0]   # Initial conditions: [A₀, B₀, P₀] = 0 mM


# =============================================================================
# 3. REACTION RATE FUNCTIONS
# =============================================================================

def rate_v1_inhibited(X, P, V1max, Km1, Ki):
    """
    Reaction v1: X → A
    Rate law: Michaelis-Menten with non-competitive inhibition by P.

    In non-competitive inhibition, the inhibitor (P) binds to an allosteric
    site on the enzyme — separate from the substrate-binding (active) site.
    This reduces the effective Vmax without changing Km.

    The inhibition factor (1 + [P]/Ki) multiplies the denominator, reducing
    the overall reaction velocity as [P] increases.
    When [P] = 0, the expression reduces to standard Michaelis-Menten.

    Formula:
              V1max · [X]
    v1 = ─────────────────────────
          (Km1 + [X]) · (1 + P/Ki)

    Parameters
    ----------
    X     : float  — external substrate concentration (mM), constant
    P     : float  — product P concentration (mM), the inhibitor
    V1max : float  — maximum velocity (mM/h)
    Km1   : float  — Michaelis constant (mM)
    Ki    : float  — inhibition constant (mM)

    Returns
    -------
    float : reaction rate v1 (mM/h)
    """
    return (V1max * X) / ((Km1 + X) * (1 + P / Ki))


def rate_v1_baseline(X, V1max, Km1):
    """
    Reaction v1: X → A — BASELINE (no inhibition).
    Standard Michaelis-Menten kinetics.

    Formula:
              V1max · [X]
    v1 = ─────────────────
              Km1 + [X]

    Parameters
    ----------
    X     : float  — external substrate concentration (mM)
    V1max : float  — maximum velocity (mM/h)
    Km1   : float  — Michaelis constant (mM)

    Returns
    -------
    float : reaction rate v1 (mM/h)
    """
    return (V1max * X) / (Km1 + X)


def rate_first_order(k, S):
    """
    Generic first-order reaction rate: v = k · [S]

    Used for reactions v2, v3, and v4 where the rate is proportional to
    the concentration of the substrate.

    Parameters
    ----------
    k : float — rate constant (1/h)
    S : float — substrate concentration (mM)

    Returns
    -------
    float : reaction rate (mM/h)
    """
    return k * S


# =============================================================================
# 4. ODE SYSTEM DEFINITIONS
# =============================================================================

def odes_inhibited(t, y, params):
    """
    ODE system WITH non-competitive allosteric feedback inhibition by P on v1.

    This is the biologically realistic model where the accumulation of product
    P slows down the entry of flux into the pathway via negative feedback.

    Mass balance equations (derived from the stoichiometric matrix S):
      d[A]/dt = v1 - v2 - v4     (A is produced by v1, consumed by v2 and v4)
      d[B]/dt = v2 - v3           (B is produced by v2, consumed by v3)
      d[P]/dt = v3                (P is only produced, no consumption/outflow)

    Parameters
    ----------
    t      : float       — current time (h), required by solve_ivp
    y      : list/array  — current state [A, B, P] (mM)
    params : dict        — dictionary of kinetic parameters

    Returns
    -------
    list : [dA/dt, dB/dt, dP/dt]
    """
    A, B, P = y  # unpack current concentrations

    # Calculate individual reaction rates at current state
    v1 = rate_v1_inhibited(params["X"], P, params["V1max"],
                           params["Km1"], params["Ki"])
    v2 = rate_first_order(params["k2"], A)
    v3 = rate_first_order(params["k3"], B)
    v4 = rate_first_order(params["k4"], A)

    # Mass balance ODEs
    dA = v1 - v2 - v4   # net change in A
    dB = v2 - v3         # net change in B
    dP = v3              # net change in P (only produced)

    return [dA, dB, dP]


def odes_baseline(t, y, params):
    """
    ODE system WITHOUT allosteric inhibition — baseline reference model.

    Identical to odes_inhibited except v1 uses standard Michaelis-Menten
    (no inhibition term). Used to compare how much the feedback inhibition
    suppresses product formation relative to an uninhibited system.

    Parameters
    ----------
    t      : float       — current time (h)
    y      : list/array  — current state [A, B, P] (mM)
    params : dict        — dictionary of kinetic parameters

    Returns
    -------
    list : [dA/dt, dB/dt, dP/dt]
    """
    A, B, P = y

    v1 = rate_v1_baseline(params["X"], params["V1max"], params["Km1"])
    v2 = rate_first_order(params["k2"], A)
    v3 = rate_first_order(params["k3"], B)
    v4 = rate_first_order(params["k4"], A)

    dA = v1 - v2 - v4
    dB = v2 - v3
    dP = v3

    return [dA, dB, dP]


# =============================================================================
# 5. RUN SIMULATION
# =============================================================================

def run_simulation(params, y0, t_start, t_end, n_points):
    """
    Solve both ODE systems (inhibited and baseline) over the defined
    time span using the Runge-Kutta RK45 method.

    RK45 (Dormand-Prince) is an explicit adaptive-step Runge-Kutta method
    of order 4(5). It automatically adjusts the internal step size to keep
    the local error below the specified tolerance, balancing accuracy and
    speed. It is the standard solver for non-stiff biological ODE systems.

    Parameters
    ----------
    params  : dict   — kinetic parameters
    y0      : list   — initial conditions [A₀, B₀, P₀]
    t_start : float  — start time (h)
    t_end   : float  — end time (h)
    n_points: int    — number of output time points

    Returns
    -------
    tuple : (sol_inhibited, sol_baseline) — scipy OdeSolution objects
    """
    t_span = (t_start, t_end)
    t_eval = np.linspace(t_start, t_end, n_points)

    print(f"Running simulation: t = {t_start} to {t_end} h  |  {n_points} points")
    print(f"Initial conditions: A={y0[0]}, B={y0[1]}, P={y0[2]} mM")
    print("-" * 60)

    # Simulate WITH inhibition
    sol_inh = solve_ivp(
        fun     = lambda t, y: odes_inhibited(t, y, params),
        t_span  = t_span,
        y0      = y0,
        t_eval  = t_eval,
        method  = "RK45",
        rtol    = 1e-8,   # relative tolerance
        atol    = 1e-10,  # absolute tolerance
    )

    # Simulate WITHOUT inhibition (baseline)
    sol_base = solve_ivp(
        fun     = lambda t, y: odes_baseline(t, y, params),
        t_span  = t_span,
        y0      = y0,
        t_eval  = t_eval,
        method  = "RK45",
        rtol    = 1e-8,
        atol    = 1e-10,
    )

    # Check solver success
    if sol_inh.success and sol_base.success:
        print("✓ Both ODE systems solved successfully.")
    else:
        print("✗ Solver error:", sol_inh.message or sol_base.message)

    return sol_inh, sol_base


# =============================================================================
# 6. COMPUTE FLUX TIME COURSES
# =============================================================================

def compute_fluxes(sol, params, inhibited=True):
    """
    Recompute the four reaction fluxes (v1–v4) at every recorded time point
    from the ODE solution. This is needed for plotting flux trajectories.

    Parameters
    ----------
    sol       : OdeSolution — output from solve_ivp
    params    : dict        — kinetic parameters
    inhibited : bool        — True = use inhibited v1 formula, False = baseline

    Returns
    -------
    tuple : (v1_t, v2_t, v3_t, v4_t) — numpy arrays of flux vs. time
    """
    A_t, B_t, P_t = sol.y   # shape: (3, n_points)

    if inhibited:
        v1_t = rate_v1_inhibited(params["X"], P_t,
                                  params["V1max"], params["Km1"], params["Ki"])
    else:
        v1_t = rate_v1_baseline(params["X"], params["V1max"], params["Km1"])
        v1_t = np.full_like(A_t, v1_t)  # constant, broadcast to array

    v2_t = rate_first_order(params["k2"], A_t)
    v3_t = rate_first_order(params["k3"], B_t)
    v4_t = rate_first_order(params["k4"], A_t)

    return v1_t, v2_t, v3_t, v4_t


# =============================================================================
# 7. PRINT SUMMARY STATISTICS
# =============================================================================

def print_summary(sol_inh, sol_base, params):
    """
    Print a summary of key simulation results to the terminal.

    Includes:
    - Final concentrations at t_end
    - Quasi-steady-state estimates (for A and B in baseline model)
    - Inhibition fold-change at t_end
    """
    t_end_val   = sol_inh.t[-1]
    A_i, B_i, P_i = sol_inh.y[:, -1]
    A_b, B_b, P_b = sol_base.y[:, -1]

    # Theoretical quasi-steady-state values (baseline, ignoring P growth)
    # At QSS: dA/dt = 0 → v1 = (k2 + k4)·A → A_qss = v1/(k2+k4)
    #         dB/dt = 0 → v2 = k3·B         → B_qss = k2·A/(k3)
    v1_const  = rate_v1_baseline(params["X"], params["V1max"], params["Km1"])
    A_qss     = v1_const / (params["k2"] + params["k4"])
    B_qss     = params["k2"] * A_qss / params["k3"]

    # Inhibition factor at final time point
    inh_factor = 1 + P_i / params["Ki"]
    v1_inh_end = rate_v1_inhibited(params["X"], P_i,
                                    params["V1max"], params["Km1"], params["Ki"])

    print("\n" + "=" * 60)
    print("  SIMULATION SUMMARY")
    print("=" * 60)
    print(f"  Simulation end time : {t_end_val:.0f} h")
    print()
    print(f"  {'Metabolite':<14} {'With Inhibition':>18} {'Baseline (no inh)':>20}")
    print(f"  {'-'*52}")
    print(f"  {'[A] (mM)':<14} {A_i:>18.4f} {A_b:>20.4f}")
    print(f"  {'[B] (mM)':<14} {B_i:>18.4f} {B_b:>20.4f}")
    print(f"  {'[P] (mM)':<14} {P_i:>18.4f} {P_b:>20.4f}")
    print()
    print(f"  Baseline QSS (A) : {A_qss:.4f} mM")
    print(f"  Baseline QSS (B) : {B_qss:.4f} mM")
    print()
    print(f"  v1 at t=0  (no P)  : {v1_const:.4f} mM/h")
    print(f"  v1 at t={t_end_val:.0f}h (inhibited) : {v1_inh_end:.4f} mM/h")
    print(f"  Inhibition fold    : {inh_factor:.2f}x slower")
    print(f"  P suppression      : {(1 - P_i/P_b)*100:.1f}% less P than baseline")
    print("=" * 60 + "\n")


# =============================================================================
# 8. PLOTTING
# =============================================================================

def plot_results(sol_inh, sol_base, params, output_path="figures/simulation_figure.png"):
    """
    Generate a 4-panel figure summarising the simulation results.

    Panel 1 (top-left)  : Metabolite dynamics WITH allosteric inhibition
    Panel 2 (top-right) : Metabolite dynamics WITHOUT inhibition (baseline)
    Panel 3 (bottom-left) : Reaction flux time courses (inhibited model)
    Panel 4 (bottom-right): Comparison of P accumulation (inhibited vs baseline)

    Parameters
    ----------
    sol_inh     : OdeSolution — inhibited model solution
    sol_base    : OdeSolution — baseline model solution
    params      : dict        — kinetic parameters (for flux recomputation)
    output_path : str         — file path to save the figure
    """
    # ── Color scheme ──────────────────────────────────────────────────────────
    C = {
        "A"    : "#2196F3",   # blue
        "B"    : "#FF5722",   # orange-red
        "P"    : "#4CAF50",   # green
        "v1"   : "#9C27B0",   # purple
        "v2"   : "#2196F3",   # blue
        "v3"   : "#4CAF50",   # green
        "v4"   : "#FF5722",   # orange-red
        "shade": "#BDBDBD",   # grey for shading
    }

    t     = sol_inh.t
    A_i   = sol_inh.y[0]
    B_i   = sol_inh.y[1]
    P_i   = sol_inh.y[2]

    # Recompute fluxes for the inhibited model
    v1_t, v2_t, v3_t, v4_t = compute_fluxes(sol_inh, params, inhibited=True)

    # ── Figure layout ─────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 10))
    fig.patch.set_facecolor("#FAFAFA")
    gs  = gridspec.GridSpec(2, 2, hspace=0.42, wspace=0.34)

    # ── Panel 1: Metabolite dynamics — inhibited ──────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(t, A_i,           color=C["A"], lw=2,   label="A (Intermediate 1)")
    ax1.plot(t, B_i,           color=C["B"], lw=2,   label="B (Intermediate 2)")
    ax1.plot(t, P_i,           color=C["P"], lw=2,   label="P (Target Product)")
    ax1.set_title("Panel 1 — Metabolite Dynamics\n(With Allosteric Feedback Inhibition)",
                  fontsize=10, fontweight="bold")
    ax1.set_xlabel("Time (h)")
    ax1.set_ylabel("Concentration (mM)")
    ax1.legend(fontsize=8, loc="upper left")
    ax1.grid(alpha=0.3)
    ax1.set_xlim(0, T_END)

    # ── Panel 2: Metabolite dynamics — baseline ───────────────────────────────
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(t, sol_base.y[0], color=C["A"], lw=2, ls="--", label="A (no inhibition)")
    ax2.plot(t, sol_base.y[1], color=C["B"], lw=2, ls="--", label="B (no inhibition)")
    ax2.plot(t, sol_base.y[2], color=C["P"], lw=2, ls="--", label="P (no inhibition)")
    ax2.set_title("Panel 2 — Metabolite Dynamics\n(Baseline: No Allosteric Inhibition)",
                  fontsize=10, fontweight="bold")
    ax2.set_xlabel("Time (h)")
    ax2.set_ylabel("Concentration (mM)")
    ax2.legend(fontsize=8, loc="upper left")
    ax2.grid(alpha=0.3)
    ax2.set_xlim(0, T_END)

    # ── Panel 3: Flux dynamics — inhibited ───────────────────────────────────
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(t, v1_t, color=C["v1"], lw=2,        label="v₁  X → A  (inhibited by P)")
    ax3.plot(t, v2_t, color=C["v2"], lw=2,        label="v₂  A → B")
    ax3.plot(t, v3_t, color=C["v3"], lw=2,        label="v₃  B → P")
    ax3.plot(t, v4_t, color=C["v4"], lw=2, ls=":", label="v₄  A → Byproduct")
    ax3.set_title("Panel 3 — Reaction Flux Over Time\n(Inhibited Model)",
                  fontsize=10, fontweight="bold")
    ax3.set_xlabel("Time (h)")
    ax3.set_ylabel("Flux (mM/h)")
    ax3.legend(fontsize=8)
    ax3.grid(alpha=0.3)
    ax3.set_xlim(0, T_END)

    # ── Panel 4: P accumulation comparison ───────────────────────────────────
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.plot(t, P_i,           color=C["P"], lw=2.5,        label="P — With Inhibition")
    ax4.plot(t, sol_base.y[2], color=C["P"], lw=2.5, ls="--", label="P — Baseline (no inhibition)")
    ax4.fill_between(t, P_i, sol_base.y[2],
                     alpha=0.12, color=C["shade"],
                     label="Difference (suppressed by feedback)")
    ax4.set_title("Panel 4 — Product P: Effect of Allosteric Feedback\n"
                  "(Inhibited vs. Baseline Comparison)",
                  fontsize=10, fontweight="bold")
    ax4.set_xlabel("Time (h)")
    ax4.set_ylabel("P Concentration (mM)")
    ax4.legend(fontsize=8)
    ax4.grid(alpha=0.3)
    ax4.set_xlim(0, T_END)

    # ── Super title ───────────────────────────────────────────────────────────
    fig.suptitle(
        "Kinetic ODE Simulation — Branched Metabolic Pathway with Allosteric Feedback\n"
        r"$v_1$ : Michaelis-Menten + Non-competitive Inhibition   |   "
        r"$v_2, v_3, v_4$ : First-order Kinetics",
        fontsize=11, y=1.02
    )

    # ── Save ──────────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else ".", exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"✓ Figure saved to: {output_path}")
    plt.close()


# =============================================================================
# 9. MAIN ENTRY POINT
# =============================================================================

def main():
    """
    Main function — runs the full simulation pipeline:
      1. Solve ODEs (inhibited and baseline)
      2. Print summary statistics to terminal
      3. Save 4-panel figure to output file
    """
    print("\n" + "=" * 60)
    print("  Metabolic Pathway Kinetic ODE Simulation")
    print("  Course: Biotechnology — Systems Biology (BISB211605)")
    print("=" * 60 + "\n")

    # Run ODE solver
    sol_inh, sol_base = run_simulation(
        params   = PARAMS,
        y0       = Y0,
        t_start  = T_START,
        t_end    = T_END,
        n_points = N_POINTS,
    )

    # Print results to terminal
    print_summary(sol_inh, sol_base, PARAMS)

    # Save figure
    plot_results(sol_inh, sol_base, PARAMS, output_path="figures/simulation_figure.png")

    print("Done. Check the 'figures/' folder for the output plot.\n")


# Standard Python idiom: only run main() when this file is executed directly,
# not when it is imported as a module by another script.
if __name__ == "__main__":
    main()
