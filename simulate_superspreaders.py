"""Monte Carlo simulator for the superspreader epidemic project.

The implementation follows Fujie and Odagaki's spatial SIR model closely:
individuals live in a two-dimensional spatial domain, infection probability
depends on distance, and superspreaders are represented either by stronger
local infectiousness or by a longer contact radius.

Running this file generates the figures used by main.tex:

    python3 simulate_superspreaders.py

The reference paper averages over 1000 Monte Carlo trials. This script also
uses MC_RUNS=1000 by default; lower it only for quick exploratory runs.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


R0_DISTANCE = 1.0
L = 10.0 * R0_DISTANCE
W0 = 1.0
DEFAULT_GAMMA = 1.0
MAX_STEPS = 60
SEED = 20250618

PERCOLATION_DISTANCE = 0.5 * L
PAPER_TRIALS = 1000
MC_RUNS = int(os.environ.get("MC_RUNS", str(PAPER_TRIALS)))

# Six-day binned SARS Singapore counts digitized for the qualitative
# comparison figure.
SARS_SINGAPORE_6DAY_COUNTS = np.array(
    [0, 0, 5, 10, 20, 50, 17, 16, 40, 27, 12, 9] + [0] * 14,
    dtype=float,
)


@dataclass
class SimulationResult:
    percolated: bool
    new_cases: list[int]
    front: list[float]
    secondary: np.ndarray
    final_attack_rate: float
    peak_time: int
    peak_size: int


def individuals_from_density(density_scaled: float) -> int:
    """Convert rho*pi*r0^2 into N for L=10*r0."""
    return max(2, int(round(density_scaled * L * L / (np.pi * R0_DISTANCE**2))))


def density_scaled_from_individuals(n: int) -> float:
    """Convert N into rho*pi*r0^2 for L=10*r0."""
    return n * np.pi * R0_DISTANCE**2 / (L * L)


def infection_probability(distance: np.ndarray, is_super: bool, model: str) -> np.ndarray:
    if model == "strong":
        if is_super:
            return np.where(distance <= R0_DISTANCE, W0, 0.0)
        return np.where(
            distance <= R0_DISTANCE,
            W0 * (1.0 - distance / R0_DISTANCE) ** 2,
            0.0,
        )

    if model == "hub":
        cutoff = np.sqrt(6.0) * R0_DISTANCE if is_super else R0_DISTANCE
        return np.where(distance <= cutoff, W0 * (1.0 - distance / cutoff) ** 2, 0.0)

    raise ValueError(f"unknown model: {model}")


def run_simulation(
    n: int,
    lambda_super: float,
    model: str,
    rng: np.random.Generator,
    gamma: float = DEFAULT_GAMMA,
    max_steps: int = MAX_STEPS,
    wrap_y: bool = True,
) -> SimulationResult:
    positions = rng.random((n, 2)) * L
    positions[0] = np.array([0.5 * L, 0.0])

    is_super = rng.random(n) < lambda_super
    states = np.zeros(n, dtype=np.int8)  # 0=S, 1=I, 2=R
    states[0] = 1

    secondary = np.zeros(n, dtype=np.int32)
    new_cases: list[int] = []
    front: list[float] = []
    percolated = False

    for _step in range(max_steps):
        infected = np.flatnonzero(states == 1)
        if infected.size == 0:
            break

        new_count = 0
        for infector in infected:
            susceptible = np.flatnonzero(states == 0)
            if susceptible.size == 0:
                break

            dx = np.abs(positions[susceptible, 0] - positions[infector, 0])
            dx = np.minimum(dx, L - dx)
            dy = np.abs(positions[susceptible, 1] - positions[infector, 1])
            if wrap_y:
                dy = np.minimum(dy, L - dy)
            distance = np.sqrt(dx * dx + dy * dy)

            probs = infection_probability(distance, bool(is_super[infector]), model)
            infected_now = susceptible[rng.random(susceptible.size) < probs]

            # Targets infected earlier in this same step are no longer susceptible.
            infected_now = infected_now[states[infected_now] == 0]
            if infected_now.size:
                states[infected_now] = 1
                secondary[infector] += infected_now.size
                new_count += int(infected_now.size)

        recover_draws = rng.random(infected.size)
        states[infected[recover_draws < gamma]] = 2

        ever_infected = np.flatnonzero(states > 0)
        if ever_infected.size:
            dx0 = np.abs(positions[ever_infected, 0] - positions[0, 0])
            dx0 = np.minimum(dx0, L - dx0)
            dy0 = np.abs(positions[ever_infected, 1] - positions[0, 1])
            if wrap_y:
                dy0 = np.minimum(dy0, L - dy0)
            direct_distance = np.sqrt(dx0 * dx0 + dy0 * dy0)
            front.append(float(direct_distance.max()))
            if direct_distance.max() >= PERCOLATION_DISTANCE:
                percolated = True
        else:
            front.append(0.0)

        new_cases.append(new_count)

    total_infected = int(np.count_nonzero(states > 0))
    peak_size = max(new_cases) if new_cases else 0
    peak_time = int(np.argmax(new_cases)) if new_cases else 0
    return SimulationResult(
        percolated=percolated,
        new_cases=new_cases,
        front=front,
        secondary=secondary,
        final_attack_rate=total_infected / n,
        peak_time=peak_time,
        peak_size=peak_size,
    )


def mean_padded(series: list[list[float]]) -> tuple[np.ndarray, np.ndarray]:
    max_len = max((len(s) for s in series), default=0)
    if max_len == 0:
        return np.array([]), np.array([])
    arr = np.full((len(series), max_len), np.nan)
    for i, values in enumerate(series):
        arr[i, : len(values)] = values
    return np.nanmean(arr, axis=0), np.nanstd(arr, axis=0) / np.sqrt(len(series))


def mean_zero_padded(series: list[list[float]], length: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Average time series after padding missing tail values with zero."""

    max_len = length if length is not None else max((len(s) for s in series), default=0)
    if max_len == 0:
        return np.array([]), np.array([])
    arr = np.zeros((len(series), max_len))
    for i, values in enumerate(series):
        clipped = values[:max_len]
        arr[i, : len(clipped)] = clipped
    return arr.mean(axis=0), arr.std(axis=0) / np.sqrt(len(series))


def theoretical_critical(lambda_values: np.ndarray, rc: float) -> np.ndarray:
    return rc / (lambda_values + (1.0 - lambda_values) / 6.0)


def interpolate_half_probability(densities: np.ndarray, probabilities: np.ndarray) -> float:
    if probabilities[0] >= 0.5:
        return float(densities[0])
    if probabilities[-1] <= 0.5:
        return float("nan")
    for left, right, p_left, p_right in zip(
        densities[:-1], densities[1:], probabilities[:-1], probabilities[1:]
    ):
        if (p_left - 0.5) * (p_right - 0.5) <= 0 and p_left != p_right:
            weight = (0.5 - p_left) / (p_right - p_left)
            return float(left + weight * (right - left))
    return float("nan")


def save_percolation_figures(out_dir: Path, rng: np.random.Generator) -> dict[str, np.ndarray]:
    n_values = np.arange(150, 901, 50, dtype=int)
    densities = np.array([density_scaled_from_individuals(int(n)) for n in n_values])
    lambdas = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], dtype=float)
    runs = MC_RUNS
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(lambdas)))
    criticals: dict[str, list[float]] = {"strong": [], "hub": []}

    for model in ("strong", "hub"):
        fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=160)
        for color, lam in zip(colors, lambdas):
            probs = []
            for n in n_values:
                hits = sum(
                    run_simulation(int(n), float(lam), model, rng).percolated
                    for _ in range(runs)
                )
                probs.append(hits / runs)

            probs_arr = np.maximum.accumulate(np.array(probs))
            criticals[model].append(interpolate_half_probability(densities, probs_arr))
            ax.plot(densities, probs_arr, marker="o", lw=1.8, ms=4, color=color, label=f"{lam:.1f}")

        ax.set_xlabel(r"scaled density $\rho\pi r_0^2$")
        ax.set_ylabel("percolation probability")
        ax.set_xlim(left=0)
        ax.set_ylim(-0.04, 1.04)
        ax.grid(alpha=0.28)
        ax.legend(title=r"$\lambda$", ncol=3, fontsize=8, title_fontsize=9)
        ax.set_title("Strong infectiousness model" if model == "strong" else "Hub model")
        fig.tight_layout()
        fig.savefig(out_dir / f"fig_{model}_percolation.png")
        plt.close(fig)

    lam_grid = np.linspace(0, 1, 200)
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=160)
    ax.plot(lam_grid, theoretical_critical(lam_grid, 4.5), color="#2563eb", lw=2.2, label="strong theory")
    ax.plot(lam_grid, theoretical_critical(lam_grid, 3.2), color="#dc2626", lw=2.2, label="hub theory")
    ax.scatter(lambdas, criticals["strong"], color="#2563eb", marker="o", label="strong simulation")
    ax.scatter(lambdas, criticals["hub"], color="#dc2626", marker="s", label="hub simulation")
    ax.set_xlabel(r"superspreader fraction $\lambda$")
    ax.set_ylabel(r"critical density $\rho_c\pi r_0^2$")
    ax.set_ylim(0, 28)
    ax.grid(alpha=0.28)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_critical_density.png")
    plt.close(fig)

    return {"densities": densities, "lambdas": lambdas}


def save_dynamic_figures(out_dir: Path, rng: np.random.Generator) -> None:
    n = individuals_from_density(20.0)
    lambdas = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    runs = MC_RUNS

    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=160)
    colors = plt.cm.plasma(np.linspace(0.08, 0.92, len(lambdas)))
    for color, lam in zip(colors, lambdas):
        fronts = [
            run_simulation(n, float(lam), "strong", rng).front
            for _ in range(runs)
        ]
        mean_front, err_front = mean_padded(fronts)
        x = np.arange(mean_front.size)
        ax.plot(x, mean_front, color=color, lw=1.9, label=f"{lam:.1f}")
        ax.fill_between(x, mean_front - err_front, mean_front + err_front, color=color, alpha=0.12)
    ax.set_xlabel("time step")
    ax.set_ylabel(r"front distance $r_f/r_0$")
    ax.grid(alpha=0.28)
    ax.legend(title=r"$\lambda$", ncol=3, fontsize=8, title_fontsize=9)
    ax.set_title("Propagation front in the strong infectiousness model")
    fig.tight_layout()
    fig.savefig(out_dir / "fig_propagation_front_strong.png")
    plt.close(fig)

    velocities = {"strong": [], "hub": []}
    velocity_err = {"strong": [], "hub": []}
    for model in ("strong", "hub"):
        for lam in lambdas:
            samples = []
            for _ in range(runs):
                front = np.array(run_simulation(n, float(lam), model, rng).front)
                if front.size >= 3:
                    t = np.arange(0, min(8, front.size))
                    slope = np.polyfit(t, front[t], deg=1)[0]
                    samples.append(max(0.0, float(slope)))
            if samples:
                velocities[model].append(float(np.mean(samples)))
                velocity_err[model].append(float(np.std(samples) / np.sqrt(len(samples))))
            else:
                velocities[model].append(0.0)
                velocity_err[model].append(0.0)

    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=160)
    ax.errorbar(lambdas, velocities["strong"], yerr=velocity_err["strong"], marker="o", lw=2, capsize=3, label="strong")
    ax.errorbar(lambdas, velocities["hub"], yerr=velocity_err["hub"], marker="s", lw=2, capsize=3, label="hub")
    ax.set_xlabel(r"superspreader fraction $\lambda$")
    ax.set_ylabel(r"propagation velocity $(r_0/\mathrm{step})$")
    ax.grid(alpha=0.28)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "fig_velocity.png")
    plt.close(fig)

    curve_specs = [("none", "strong", 0.0), ("strong", "strong", 0.2), ("hub", "hub", 0.2)]
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=160)
    for label, model, lam in curve_specs:
        curves = [
            run_simulation(n, lam, model, rng).new_cases
            for _ in range(runs)
        ]
        mean_curve, err_curve = mean_zero_padded(curves)
        x = np.arange(mean_curve.size)
        ax.plot(x, mean_curve, lw=2, label=label)
        ax.fill_between(x, mean_curve - err_curve, mean_curve + err_curve, alpha=0.13)
    ax.set_xlabel("time step")
    ax.set_ylabel("new infections")
    ax.grid(alpha=0.28)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "fig_epidemic_curves.png")
    plt.close(fig)


def save_sars_comparison_figure(out_dir: Path, rng: np.random.Generator) -> None:
    """Experiment 4: qualitative comparison with SARS Singapore data."""

    runs = MC_RUNS
    n = individuals_from_density(15.0)
    max_steps = SARS_SINGAPORE_6DAY_COUNTS.size
    specs = [
        ("no superspreaders", "strong", 0.0, "#22c7d6", "^"),
        ("strong, lambda=0.4", "strong", 0.4, "#dc2626", "o"),
        ("hub, lambda=0.4", "hub", 0.4, "#2563eb", "s"),
    ]

    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=160)
    x = np.arange(max_steps)

    # Observed SARS data: each model time step corresponds to six days.
    ax.bar(
        x,
        SARS_SINGAPORE_6DAY_COUNTS,
        width=0.82,
        color="#f59e0b",
        alpha=0.72,
        label="SARS Singapore data",
    )

    for label, model, lam, color, marker in specs:
        # Use wrap_y=False so the simulated outbreak has a bottom-to-top
        # direction, matching the qualitative setup in the reference comparison.
        curves = [
            run_simulation(
                n,
                lam,
                model,
                rng,
                max_steps=max_steps,
                wrap_y=False,
            ).new_cases
            for _ in range(runs)
        ]
        mean_curve, err_curve = mean_zero_padded(curves, length=max_steps)
        ax.plot(x, mean_curve, color=color, lw=1.8, marker=marker, ms=3.5, label=label)
        ax.fill_between(
            x,
            mean_curve - err_curve,
            mean_curve + err_curve,
            color=color,
            alpha=0.11,
            linewidth=0,
        )

    ax.set_xlabel("time step (1 step = 6 days)")
    ax.set_ylabel("new patients / infections")
    ax.set_xlim(-0.5, 25.5)
    ax.set_ylim(0, 80)
    ax.grid(axis="y", alpha=0.28)
    ax.legend(fontsize=8)
    ax.set_title("SARS Singapore epidemic curve comparison")
    fig.tight_layout()
    fig.savefig(out_dir / "fig_sars_epidemic_comparison.png")
    plt.close(fig)


def save_secondary_figures(out_dir: Path, rng: np.random.Generator) -> None:
    """Experiment 3: secondary-infection count distributions."""

    runs = MC_RUNS
    n = individuals_from_density(15.0)

    # Collect every individual's number of secondary infections across runs.
    secondary_no_super = []
    secondary_strong = []
    secondary_hub = []
    for _ in range(runs):
        secondary_no_super.extend(run_simulation(n, 0.0, "strong", rng).secondary.tolist())
        secondary_strong.extend(run_simulation(n, 0.2, "strong", rng).secondary.tolist())
        secondary_hub.extend(run_simulation(n, 0.2, "hub", rng).secondary.tolist())

    max_k = 24
    bins = np.arange(max_k + 2) - 0.5

    # Baseline distribution without superspreaders.
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=160)
    ax.hist(secondary_no_super, bins=bins, density=True, color="#64748b", edgecolor="white")
    ax.set_xlim(-0.5, max_k + 0.5)
    ax.set_xlabel("number of secondary infections")
    ax.set_ylabel("probability")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 1.0)
    ax.grid(axis="y", alpha=0.28)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_link_distribution_no_superspreaders.png")
    plt.close(fig)

    # Superspreader distributions. Log scale makes the long tail visible.
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=160)
    ax.hist(secondary_strong, bins=bins, density=True, histtype="step", lw=2.2, label="strong")
    ax.hist(secondary_hub, bins=bins, density=True, histtype="step", lw=2.2, label="hub")
    ax.set_xlim(-0.5, max_k + 0.5)
    ax.set_xlabel("number of secondary infections")
    ax.set_ylabel("probability")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 1.0)
    ax.grid(axis="y", alpha=0.28)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "fig_link_distribution_superspreaders.png")
    plt.close(fig)


def main() -> None:
    """Generate all figures used in the report."""

    out_dir = Path("figures")
    out_dir.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)

    # The four experiment blocks are independent and can be assigned separately.
    save_percolation_figures(out_dir, rng)
    save_dynamic_figures(out_dir, rng)
    save_secondary_figures(out_dir, rng)
    save_sars_comparison_figure(out_dir, rng)
    print(f"Generated figures in {out_dir.resolve()}")


if __name__ == "__main__":
    main()
