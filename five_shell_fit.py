"""Fit a five-shell vesicle SAXS model with thickness symmetry only.

Research-team requirements:
    n = 5 shells
    thickness1 == thickness5
    thickness2 == thickness4
    thickness3 is independent
    all five shell SLDs are independent

"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from bumps.fitproblem import FitProblem
from bumps.fitters import fit
from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
from sasmodels.data import Data1D

from vesicle import analyze_vesicle

# Original settings from three shell model
HERE = Path(__file__).resolve().parent
DATA_FILE = HERE / "DL_0_100_T20.0C_abs.dat"
Q_MIN, Q_MAX = 0.005, 0.7  # Å^-1
RHO_WATER_SLD = 9.41       # 10^-6 Å^-2, as supplied by research team


def load_saxs_data(path):
    """Loads the experimental SAXs Data
    Return q, I, dI after checking the research team's 3-column format"""
    if not path.is_file():
        raise FileNotFoundError(f"Missing SAXS data file: {path}")

    raw = np.loadtxt(path, ndmin=2) # loads measurements into an array
    if raw.shape[1] < 3:
        raise ValueError("SAXS data must have at least 3 columns: q, I, dI.")

    # Select rows with the first column in the applicable range
    selected = raw[(raw[:, 0] >= Q_MIN) & (raw[:, 0] <= Q_MAX), :3]
    if len(selected) == 0:
        raise ValueError("No SAXS data points within the requested q range.")
    if not np.all(np.isfinite(selected)):
        raise ValueError("Selected SAXS data contain non-finite entries.")
    q, intensity, uncertainty = selected.T
    if np.any(q <= 0) or np.any(uncertainty <= 0):
        raise ValueError("q and intensity uncertainty dI must be positive.")
    return q, intensity, uncertainty


def build_model():
    """Construct the five-shell model and tie only the shell thicknesses."""
    model = Model(
        load_model("core_multi_shell"),
        scale=0.04,
        background=0.0005,
        radius=305.0,
        radius_pd=0.22,
        radius_pd_type="schulz",
        sld_core=RHO_WATER_SLD,
        sld_solvent=RHO_WATER_SLD,
        n=5,
        # Total shell thickness = 10+5+9+5+10 = 39 Å
        thickness1=10.0,
        thickness2=5.0,
        thickness3=9.0,
        thickness4=5.0,
        thickness5=10.0,
        # All five SLDS must be independent
        sld1=11.0,
        sld2=10.2,
        sld3=8.1,
        sld4=10.1,
        sld5=11.1,
    )

    # Initial fitting bounds for optimizer
    model.scale.range(0.03, 0.06)
    model.background.range(0.0004, 0.0006)
    model.radius.range(300.0, 700.0)
    model.radius_pd.range(0.20, 0.30)

    # Only three independent thicknesses.
    model.thickness1.range(5.0, 15.0)
    model.thickness2.range(3.0, 15.0)
    model.thickness3.range(5.0, 30.0)

    # Enforcing the thickness symmetry constraint
    model.thickness5.equals(model.thickness1)
    model.thickness4.equals(model.thickness2)

    # Five independent SLDs: set ranges for all five
    model.sld1.range(9.8, 12.0)
    model.sld2.range(8.0, 12.0)
    model.sld3.range(7.0, 9.0)
    model.sld4.range(8.0, 12.0)
    model.sld5.range(9.8, 12.0)
    return model


def main():
    # Load the experimental data
    q, intensity, uncertainty = load_saxs_data(DATA_FILE)
    print(f"Loaded {len(q)} SAXS points from {DATA_FILE.name}.")

    # Convert into a data object compatible with sasmodels
    data = Data1D(x=q, y=intensity, dy=uncertainty)

    # Build the five shell model
    model = build_model()

    # Performing optimization
    experiment = Experiment(data=data, model=model)

    # Create the optimizaton problem
    problem = FitProblem(experiment)

    # ExperimentalData + model > Fit problem

    print("Running five-shell SAXS fit...")
    result = fit(problem, method="amoeba", steps=3000)
    problem.setp(result.x)  # Make the best-fit parameter vector explicit
    experiment.update()

    fitted_I = np.asarray(experiment.theory())
    chisq = problem.chisq()  # Bumps' reported normalized chi-square

    print("\nINDEPENDENT FITTED PARAMETERS")
    for name, value in zip(problem.labels(), problem.getp()):
        print(f"{name:25s} {value:.8g}")
    print(f"{'chi-square':25s} {chisq:.8g}")

 
    thicknesses = [getattr(model, f"thickness{i}").value for i in range(1, 6)]
    slds = [getattr(model, f"sld{i}").value for i in range(1, 6)]
    print("\nALL SHELL VALUES")
    for i, (thickness, sld) in enumerate(zip(thicknesses, slds), start=1):
        print(f"shell {i}: thickness = {thickness:.6f} Å, SLD = {sld:.6f} (10^-6 Å^-2)")

    # Reuse the geometry calculator
    geo = analyze_vesicle(model.radius.value, thicknesses)

    # Stop execution if geometry check fails
    if not (geo["symmetry_valid"] and geo["volume_valid"]):
        raise RuntimeError("Fitted geometry failed the checks. Inspect the parameter symmetry.")

    # Save results
    output_base = HERE / (DATA_FILE.stem.replace(".", "p") + "_five_shell")
    np.savetxt(
        str(output_base) + "_fit.csv",
        np.column_stack([q, intensity, uncertainty, fitted_I,
                         (intensity - fitted_I) / uncertainty]),
        delimiter=",",
        header="q_A^-1,I_observed,dI,I_fitted,standardized_residual",
        comments="",
    )

    fig, ax = plt.subplots(figsize=(9, 7))
    positive = (intensity > 0) & (fitted_I > 0)
    ax.scatter(q[positive], intensity[positive], s=16, c="black",
               label="Experimental data", zorder=2)
    ax.plot(q[positive], fitted_I[positive], linewidth=2.0, c="red",
            label="Five-shell fit", zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$q\;(\mathrm{\AA}^{-1})$")
    ax.set_ylabel(r"$I(q)$")
    ax.set_title("Five-shell vesicle SAXS fit (thickness symmetry enforced)")
    ax.text(0.04, 0.96, f"Bumps chi-square = {chisq:.4g}",
            transform=ax.transAxes, va="top")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(output_base) + "_fit.png", dpi=300)
    plt.close(fig)
    print(f"\nSaved: {output_base}_fit.png")
    print(f"Saved: {output_base}_fit.csv")

if __name__ == "__main__":
    main()
