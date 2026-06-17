# Kinetic ODE Simulation of a Branched Metabolic Pathway
### with Non-Competitive Allosteric Feedback Inhibition

> **Course:** Biotechnology — Systems Biology (BISB211605)  
> **Name:** Adelia Yusuf Ardhani.  
> **Institution:** Faculty of Biology, Universitas Gadjah Mada  
> **Semester:** Even Semester 2025/2026

---

## Overview

This repository contains a Python simulation of a **4-reaction branched metabolic pathway** modeled using **Ordinary Differential Equations (ODEs)**. The simulation explores how a target product P accumulates over a 48-hour fermentation cycle, with and without allosteric feedback inhibition.

The pathway represents a common scenario in **metabolic engineering**: an external substrate X is converted through two intermediates (A and B) into a target product P, while a competing branch reaction (v₄) drains flux toward an unwanted byproduct.

---

## Pathway Topology

```
                     Allosteric Inhibition (non-competitive)
        ┌──────────────────────────────────────────────────────┐
        │                                                      │
        ▼                                                      │
  X ──[v₁]──► A ──[v₂]──► B ──[v₃]──► P (Target Product)
                   │
                  [v₄]
                   │
                   ▼
              Byproduct (unwanted)
```

| Reaction | Description | Rate Law |
|----------|-------------|----------|
| **v₁** | X → A (entry step) | Michaelis-Menten + non-competitive inhibition by P |
| **v₂** | A → B | First-order: v₂ = k₂ · [A] |
| **v₃** | B → P | First-order: v₃ = k₃ · [B] |
| **v₄** | A → Byproduct | First-order: v₄ = k₄ · [A] |

**Key features:**
- Substrate X is maintained at a **constant external concentration** (not a dynamic variable)
- Internal pools A, B, and P are **dynamic** (change over time)
- Product P exerts **non-competitive allosteric feedback inhibition** on v₁, reducing effective Vmax as [P] increases

---

## Scientific Background

### Michaelis-Menten Kinetics

Enzyme-catalyzed reactions often follow Michaelis-Menten kinetics, where the reaction rate saturates at high substrate concentrations:

$$v = \frac{V_{max} \cdot [S]}{K_m + [S]}$$

- **V_max** — maximum reaction velocity (when all enzyme active sites are occupied)
- **K_m** — Michaelis constant (substrate concentration at half-maximal rate; reflects enzyme-substrate affinity)

### Non-Competitive Allosteric Inhibition

In non-competitive inhibition, the inhibitor (product P) binds to an **allosteric site** — a site separate from the substrate-binding active site. This does not block substrate binding (K_m is unchanged) but reduces the effective V_max by altering enzyme conformation.

The rate equation for v₁ with non-competitive inhibition is:

$$v_1 = \frac{V_{1,max} \cdot [X]}{(K_{m1} + [X]) \cdot \left(1 + \dfrac{[P]}{K_i}\right)}$$

The term $(1 + [P]/K_i)$ is the **inhibition factor**:
- When [P] = 0 → inhibition factor = 1 → standard Michaelis-Menten
- When [P] = K_i → inhibition factor = 2 → v₁ is halved
- As [P] → ∞ → inhibition factor → ∞ → v₁ approaches zero

This creates a **negative feedback loop**: more P → slower v₁ → less A → less flux through v₂/v₃ → slower P production. The system is self-limiting.

### Stoichiometric Matrix

The mass balance ODEs are derived directly from the **stoichiometric matrix S**, where rows represent internal metabolites and columns represent reactions:

$$S = \begin{pmatrix} & v_1 & v_2 & v_3 & v_4 \\ A & +1 & -1 & 0 & -1 \\ B & 0 & +1 & -1 & 0 \\ P & 0 & 0 & +1 & 0 \end{pmatrix}$$

Reading each row gives the ODE:

$$\frac{d[A]}{dt} = v_1 - v_2 - v_4$$

$$\frac{d[B]}{dt} = v_2 - v_3$$

$$\frac{d[P]}{dt} = v_3$$

Note: d[P]/dt has no negative term because no reaction consumes P — it accumulates monotonically, which drives the feedback inhibition.

---

## Parameters

| Parameter | Symbol | Value | Unit | Description |
|-----------|--------|-------|------|-------------|
| Max velocity of v₁ | V₁_max | 5.0 | mM/h | Maximum rate when enzyme is fully saturated |
| Michaelis constant | K_m1 | 2.0 | mM | Substrate [X] at half-maximal v₁ |
| Inhibition constant | K_i | 3.0 | mM | [P] at which v₁ is reduced by 50% |
| External substrate | X | 10.0 | mM | Held constant (not a dynamic variable) |
| Rate constant v₂ | k₂ | 1.0 | 1/h | First-order rate for A → B |
| Rate constant v₃ | k₃ | 0.8 | 1/h | First-order rate for B → P |
| Rate constant v₄ | k₄ | 0.3 | 1/h | First-order rate for A → Byproduct |

**Initial conditions:** A(0) = B(0) = P(0) = 0 mM  
**Simulation duration:** 48 hours

---

## Repository Structure

```
.
├── simulation.py        # Main simulation script (all code lives here)
├── requirements.txt     # Python package dependencies
├── README.md            # This file
└── figures/             # Output folder — created automatically when you run the script
    └── simulation_figure.png
```

---

## Installation and Usage

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/<your-repo-name>.git
cd <your-repo-name>
```

### 2. (Recommended) Create a virtual environment

```bash
python -m venv venv
source venv/bin/activate        # on macOS/Linux
venv\Scripts\activate           # on Windows
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Run the simulation

```bash
python simulation.py
```

The script will:
1. Solve both ODE systems (inhibited and baseline) over 48 hours
2. Print a summary table of results to the terminal
3. Save a 4-panel figure to `figures/simulation_figure.png`

---

## Output

### Terminal Output (example)

```
============================================================
  Metabolic Pathway Kinetic ODE Simulation
  Course: Biotechnology — Systems Biology (BISB211605)
============================================================

Running simulation: t = 0 to 48 h  |  1000 points
Initial conditions: A=0.0, B=0.0, P=0.0 mM
------------------------------------------------------------
✓ Both ODE systems solved successfully.

============================================================
  SIMULATION SUMMARY
============================================================
  Metabolite        With Inhibition    Baseline (no inh)
  ----------------------------------------------------
  [A] (mM)                   0.3078               3.2051
  [B] (mM)                   0.3898               4.0064
  [P] (mM)                  28.4793             147.3743

  Baseline QSS (A) : 3.2051 mM
  Baseline QSS (B) : 4.0064 mM

  v1 at t=0  (no P)  : 4.1667 mM/h
  v1 at t=48h (inhibited) : 0.3971 mM/h
  Inhibition fold    : 10.49x slower
  P suppression      : 80.7% less P than baseline
============================================================
```

### Figure Output

The script generates a 4-panel figure saved to `figures/simulation_figure.png`:

| Panel | Content |
|-------|---------|
| **Panel 1** (top-left) | Metabolite dynamics A, B, P over 48h **with** allosteric inhibition |
| **Panel 2** (top-right) | Metabolite dynamics A, B, P over 48h **without** inhibition (baseline) |
| **Panel 3** (bottom-left) | Reaction flux time courses v₁–v₄ (inhibited model) |
| **Panel 4** (bottom-right) | Comparison of P accumulation: inhibited vs. baseline |

---

## Key Results and Interpretation

### 1. Negative Feedback Strongly Suppresses Product Accumulation

At t = 48h, the inhibited model produces only **28.5 mM** of P compared to **147.4 mM** in the baseline — an **80.7% reduction**. This demonstrates that the allosteric feedback loop is a powerful regulatory mechanism that significantly limits product titer in the bioreactor.

### 2. The Entry Flux v₁ Declines Dramatically Over Time

In the inhibited model, v₁ drops from **4.17 mM/h** at t = 0 (when [P] = 0) to only **0.40 mM/h** at t = 48h — a **10.5-fold reduction**. This is the direct consequence of the inhibition factor $(1 + [P]/K_i)$ growing as P accumulates.

### 3. Intermediates A and B Reach Quasi-Steady State

In the inhibited model, [A] and [B] quickly reach a quasi-steady state (around 0.31 mM and 0.39 mM respectively), whereas in the baseline model they stabilize at much higher values (3.21 mM and 4.01 mM). The difference arises because in the inhibited model, reduced v₁ input limits how much flux flows through the entire pathway.

### 4. P Never Reaches True Steady State

Because there is no consumption or export reaction for P in this model (d[P]/dt = v₃ > 0 always), P accumulates indefinitely. In a real bioreactor, this would be addressed by product export, downstream reactions, or in situ product removal (ISPR).

---

## Metabolic Engineering Implications

From a practical standpoint, the feedback inhibition by P limits production efficiency. Potential engineering strategies to overcome this include:

- **Enzyme engineering (directed evolution):** Mutate the enzyme catalyzing v₁ to reduce its sensitivity to P (increase K_i), maintaining high flux even at elevated [P]
- **In situ product removal (ISPR):** Continuously extract P from the bioreactor (e.g., via membrane separation or two-phase systems) to keep [P] below the inhibitory threshold
- **Overexpression of v₃:** Increase expression of the enzyme catalyzing B → P to pull flux forward, reducing [B] buildup
- **Knockout of v₄:** Delete or repress the gene encoding the enzyme responsible for v₄ to eliminate byproduct formation and redirect all flux from A toward B and P

---

## References

- Klipp, E., et al. (2016). *Systems Biology: A Textbook* (2nd ed.). Wiley-Blackwell.
- Nielsen, J., & Keasling, J. D. (2016). Engineering Cellular Metabolism. *Cell*, 164(6), 1185–1197. https://doi.org/10.1016/j.cell.2016.02.004
- Orth, J. D., Thiele, I., & Palsson, B. Ø. (2010). What is flux balance analysis? *Nature Biotechnology*, 28(3), 245–248. https://doi.org/10.1038/nbt.1614
- SciPy Documentation — `scipy.integrate.solve_ivp`: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
- Course GitHub reference: https://github.com/lab-biotek-bio-ugm/S1_BISB211605_Biotechnology

---

## License

This project is submitted as part of academic coursework at Universitas Gadjah Mada. For educational use only.
