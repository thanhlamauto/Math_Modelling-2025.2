"""Animate infection-route simulations for the superspreader project.

The output is a single GIF with three panels, matching the qualitative setup
of Figs. 9--11 in Fujie and Odagaki:

    python3 visualize_superspreaders_gif.py

The generated file is:

    figures/superspreader_route_simulation.gif
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from simulate_superspreaders import (
    DEFAULT_GAMMA,
    L,
    R0_DISTANCE,
    W0,
    SEED,
    individuals_from_density,
    infection_probability,
)


DENSITY_SCALED = 15.0
LAMBDA_SUPER = 0.2
MAX_STEPS = 45
FPS = 4


@dataclass
class RouteResult:
    title: str
    model: str
    lambda_super: float
    positions: np.ndarray
    is_super: np.ndarray
    infection_time: np.ndarray
    state_history: list[np.ndarray]
    edges: list[tuple[int, int, int]]
    new_cases: list[int]
    wrap_x: bool
    wrap_y: bool

    @property
    def final_infected(self) -> int:
        return int(np.count_nonzero(self.infection_time >= 0))

    @property
    def final_step(self) -> int:
        if not self.new_cases:
            return 0
        return len(self.new_cases)

    def state_at(self, frame: int) -> np.ndarray:
        frame = min(frame, len(self.state_history) - 1)
        return self.state_history[frame]


def wrapped_distance_to(
    positions: np.ndarray,
    source: int,
    targets: np.ndarray,
    wrap_x: bool = True,
    wrap_y: bool = True,
) -> np.ndarray:
    """Compute distance from one source to many targets."""

    dx = np.abs(positions[targets, 0] - positions[source, 0])
    dy = np.abs(positions[targets, 1] - positions[source, 1])
    if wrap_x:
        dx = np.minimum(dx, L - dx)
    if wrap_y:
        dy = np.minimum(dy, L - dy)
    return np.sqrt(dx * dx + dy * dy)


def run_route_simulation(
    title: str,
    n: int,
    lambda_super: float,
    model: str,
    seed: int,
    gamma: float = DEFAULT_GAMMA,
    max_steps: int = MAX_STEPS,
    wrap_x: bool = True,
    wrap_y: bool = True,
) -> RouteResult:
    """Run one simulation and record who infected whom at each time step."""

    rng = np.random.default_rng(seed)
    positions = rng.random((n, 2)) * L
    positions[0] = np.array([0.5 * L, 0.0])

    is_super = rng.random(n) < lambda_super
    states = np.zeros(n, dtype=np.int8)
    states[0] = 1

    infection_time = np.full(n, -1, dtype=np.int16)
    infection_time[0] = 0

    state_history = [states.copy()]
    edges: list[tuple[int, int, int]] = []
    new_cases: list[int] = []

    for step in range(max_steps):
        infected = np.flatnonzero(states == 1)
        if infected.size == 0:
            break

        new_count = 0
        for infector in infected:
            susceptible = np.flatnonzero(states == 0)
            if susceptible.size == 0:
                break

            distance = wrapped_distance_to(positions, infector, susceptible, wrap_x, wrap_y)
            probs = infection_probability(distance, bool(is_super[infector]), model)
            infected_now = susceptible[rng.random(susceptible.size) < probs]
            infected_now = infected_now[states[infected_now] == 0]

            if infected_now.size:
                states[infected_now] = 1
                infection_time[infected_now] = step + 1
                new_count += int(infected_now.size)
                edges.extend((int(infector), int(target), step + 1) for target in infected_now)

        recover_draws = rng.random(infected.size)
        states[infected[recover_draws < gamma]] = 2
        new_cases.append(new_count)
        state_history.append(states.copy())

    return RouteResult(
        title=title,
        model=model,
        lambda_super=lambda_super,
        positions=positions,
        is_super=is_super,
        infection_time=infection_time,
        state_history=state_history,
        edges=edges,
        new_cases=new_cases,
        wrap_x=wrap_x,
        wrap_y=wrap_y,
    )


def find_representative_run(
    title: str,
    n: int,
    lambda_super: float,
    model: str,
    base_seed: int,
    minimum_infected: int,
    maximum_infected: int | None = None,
    attempts: int = 80,
) -> RouteResult:
    """Search deterministic seeds until the animation has a readable route."""

    best: RouteResult | None = None
    for offset in range(attempts):
        result = run_route_simulation(title, n, lambda_super, model, base_seed + offset)
        if best is None or result.final_infected > best.final_infected:
            best = result
        if result.final_infected < minimum_infected:
            continue
        if maximum_infected is not None and result.final_infected > maximum_infected:
            continue
        return result

    if best is None:
        raise RuntimeError("no route simulation was generated")
    return best


def draw_route_panel(ax: plt.Axes, result: RouteResult, frame: int) -> None:
    """Draw one panel of the current animation frame."""

    ax.clear()
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    visible_edges = [(i, j, t) for i, j, t in result.edges if t <= frame]
    if visible_edges:
        starts = np.array([result.positions[i] for i, _, _ in visible_edges])
        ends = np.array([result.positions[j] for _, j, _ in visible_edges])
        delta = ends - starts
        if result.wrap_x:
            delta[:, 0] = np.where(delta[:, 0] > 0.5 * L, delta[:, 0] - L, delta[:, 0])
            delta[:, 0] = np.where(delta[:, 0] < -0.5 * L, delta[:, 0] + L, delta[:, 0])
        if result.wrap_y:
            delta[:, 1] = np.where(delta[:, 1] > 0.5 * L, delta[:, 1] - L, delta[:, 1])
            delta[:, 1] = np.where(delta[:, 1] < -0.5 * L, delta[:, 1] + L, delta[:, 1])
        ax.quiver(
            starts[:, 0],
            starts[:, 1],
            delta[:, 0],
            delta[:, 1],
            angles="xy",
            scale_units="xy",
            scale=1,
            color="black",
            alpha=0.52,
            width=0.0021,
            headwidth=3.8,
            headlength=5.0,
            headaxislength=4.3,
            zorder=1,
        )

    states = result.state_at(frame)
    susceptible = states == 0
    active_infected = states == 1
    recovered = states == 2
    ever_infected = result.infection_time >= 0

    susceptible_normal = susceptible & ~result.is_super
    susceptible_super = susceptible & result.is_super
    active_normal = active_infected & ~result.is_super
    active_super = active_infected & result.is_super
    recovered_normal = recovered & ~result.is_super
    recovered_super = recovered & result.is_super

    ax.scatter(
        result.positions[susceptible_normal, 0],
        result.positions[susceptible_normal, 1],
        s=13,
        facecolors="white",
        edgecolors="black",
        linewidths=0.55,
        zorder=2,
    )
    ax.scatter(
        result.positions[susceptible_super, 0],
        result.positions[susceptible_super, 1],
        s=15,
        facecolors="black",
        edgecolors="black",
        linewidths=0.45,
        zorder=3,
    )
    ax.scatter(
        result.positions[recovered_normal, 0],
        result.positions[recovered_normal, 1],
        s=15,
        facecolors="white",
        edgecolors="#0b45ff",
        linewidths=0.8,
        zorder=4,
    )
    ax.scatter(
        result.positions[recovered_super, 0],
        result.positions[recovered_super, 1],
        s=18,
        facecolors="#0b45ff",
        edgecolors="#0b45ff",
        linewidths=0.6,
        zorder=5,
    )
    ax.scatter(
        result.positions[active_normal, 0],
        result.positions[active_normal, 1],
        s=20,
        facecolors="white",
        edgecolors="#d71920",
        linewidths=1.0,
        zorder=6,
    )
    ax.scatter(
        result.positions[active_super, 0],
        result.positions[active_super, 1],
        s=23,
        facecolors="#d71920",
        edgecolors="#d71920",
        linewidths=0.7,
        zorder=7,
    )

    new_cases = result.new_cases[frame - 1] if 0 < frame <= len(result.new_cases) else 0
    ax.set_title(
        f"{result.title}  (model={result.model}, lambda={result.lambda_super:g})\n"
        f"t={frame:02d}, new={new_cases}, I={int(active_infected.sum())}, "
        f"R={int(recovered.sum())}, ever={int(ever_infected.sum())}/{len(states)}",
        fontsize=10,
        pad=8,
    )


def make_gif(out_path: Path) -> None:
    """Generate the comparison animation."""

    out_path.parent.mkdir(exist_ok=True)
    n = individuals_from_density(DENSITY_SCALED)

    panels = [
        find_representative_run(
            "Strong infectiousness",
            n,
            LAMBDA_SUPER,
            "strong",
            SEED + 100,
            minimum_infected=180,
        ),
        find_representative_run(
            "Hub model",
            n,
            LAMBDA_SUPER,
            "hub",
            SEED + 500,
            minimum_infected=260,
        ),
        find_representative_run(
            "No superspreaders",
            n,
            0.0,
            "strong",
            SEED + 900,
            minimum_infected=35,
            maximum_infected=180,
        ),
    ]

    total_frames = max(panel.final_step for panel in panels) + 10
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.85), dpi=130)

    legend_handles = [
        Line2D([0], [0], color="black", lw=1.2, marker=">", markersize=5, label="route of infection"),
        Line2D([0], [0], marker="o", color="black", markerfacecolor="black", lw=0, label="S superspreader"),
        Line2D([0], [0], marker="o", color="black", markerfacecolor="white", lw=0, label="S normal"),
        Line2D([0], [0], marker="o", color="#d71920", markerfacecolor="#d71920", lw=0, label="I superspreader"),
        Line2D([0], [0], marker="o", color="#d71920", markerfacecolor="white", lw=0, label="I normal"),
        Line2D([0], [0], marker="o", color="#0b45ff", markerfacecolor="#0b45ff", lw=0, label="R superspreader"),
        Line2D([0], [0], marker="o", color="#0b45ff", markerfacecolor="white", lw=0, label="R normal"),
    ]
    setting_text = (
        rf"Settings: $N={n}$, $L=10r_0$, $r_0={R0_DISTANCE:g}$, "
        rf"$w_0={W0:g}$, $\gamma={DEFAULT_GAMMA:g}$, "
        rf"$\rho\pi r_0^2={DENSITY_SCALED:g}$, max steps={MAX_STEPS}, "
        "periodic boundaries"
    )
    update_text = (
        r"Update rule: infected-at-start transmit first, then recover with probability "
        rf"$\gamma$; newly infected transmit from the next Monte Carlo step"
    )

    def update(frame: int) -> list[plt.Artist]:
        for ax, panel in zip(axes, panels):
            draw_route_panel(ax, panel, min(frame, panel.final_step))
        fig.suptitle(
            "Route of infection animation with recovery",
            y=0.992,
            fontsize=13,
            fontweight="bold",
        )
        fig.text(0.5, 0.942, setting_text, ha="center", va="top", fontsize=9)
        fig.text(0.5, 0.914, update_text, ha="center", va="top", fontsize=8.2)
        fig.legend(
            handles=legend_handles,
            loc="lower center",
            ncol=7,
            frameon=False,
            bbox_to_anchor=(0.5, 0.005),
            fontsize=8.2,
        )
        fig.tight_layout(rect=(0.015, 0.09, 0.985, 0.875))
        return []

    animation = FuncAnimation(fig, update, frames=range(total_frames), interval=1000 / FPS)
    animation.save(out_path, writer=PillowWriter(fps=FPS))
    plt.close(fig)


def export_gif_frames(gif_path: Path, frames_dir: Path) -> None:
    """Export GIF frames for LaTeX Beamer's animate package."""

    frames_dir.mkdir(parents=True, exist_ok=True)
    for old_frame in frames_dir.glob("frame_*.png"):
        old_frame.unlink()

    with Image.open(gif_path) as gif:
        for frame_index in range(gif.n_frames):
            gif.seek(frame_index)
            frame = gif.convert("RGBA")
            frame.save(frames_dir / f"frame_{frame_index}.png")
            if frame_index == gif.n_frames - 1:
                frame.save(gif_path.parent / "superspreader_route_final.png")


def main() -> None:
    gif_path = Path("figures/superspreader_route_simulation.gif")
    make_gif(gif_path)
    export_gif_frames(gif_path, Path("figures/superspreader_route_frames"))


if __name__ == "__main__":
    main()
