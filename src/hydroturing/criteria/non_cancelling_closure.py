"""Block-wise water closure without cancellation between output intervals."""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import (
    FAIL, PASS, CriterionResult, criterion, make_window, reported_states, segments,
)
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


@criterion("non_cancelling_closure")
def non_cancelling_closure(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Score created and lost water separately in each prescribed block.

    Positive residual means unexplained loss. Negative residual means creation.
    Every declared water store is counted, including optional groundwater and
    channel storage. The first scored interval starts from the last spinup state.
    """
    allowed = {"threshold", "absolute_tolerance_mm", "step_deadband_mm", "segment_column"}
    unknown = set(params) - allowed
    if unknown:
        raise ValueError(f"non_cancelling_closure: unknown parameters {sorted(unknown)}")
    threshold = float(params.get("threshold", 0.05))
    floor = float(params.get("absolute_tolerance_mm", 0.01))
    deadband = float(params.get("step_deadband_mm", 0.001))
    if not np.isfinite([threshold, floor, deadband]).all() or min(threshold, floor, deadband) < 0:
        raise ValueError("non_cancelling_closure tolerances must be finite and nonnegative")
    if threshold == floor == 0:
        raise ValueError("non_cancelling_closure needs a positive block allowance")

    w = make_window(run, probe)
    states = reported_states(w, probe)
    required = ["pr", "evspsbl", "mrro", *states]
    missing = [name for name in required if name not in w.table]
    if missing:
        raise ValueError(f"non_cancelling_closure needs model result columns {missing}")
    if "pr" not in w.forcing:
        raise ValueError("non_cancelling_closure needs supplied precipitation")

    forcing_pr = np.asarray(w.forcing["pr"], dtype=float)
    names = ["evspsbl", "mrro", *states]
    if "gwex" in w.table:
        names.append("gwex")
    values = [forcing_pr, *(np.asarray(w.table[name], dtype=float) for name in names)]
    initial = w.storage_initial(states)
    if not np.isfinite(initial) or not all(np.isfinite(value).all() for value in values):
        return CriterionResult(
            name="non_cancelling_closure", status=FAIL, threshold=1.0,
            message="non-finite water budget data", diagnostics={"states": list(states)},
        )
    if (forcing_pr < 0).any():
        raise ValueError("non_cancelling_closure needs nonnegative precipitation")

    p = forcing_pr * w.dt_days
    et = np.asarray(w.table["evspsbl"], dtype=float) * w.dt_days
    q = np.asarray(w.table["mrro"], dtype=float) * w.dt_days
    x = np.asarray(w.table["gwex"], dtype=float) * w.dt_days if "gwex" in w.table else np.zeros(len(p))
    storage = w.storage(states)
    residual = p + x - et - q - np.diff(storage, prepend=initial)
    if not np.isfinite(residual).all():
        return CriterionResult(
            name="non_cancelling_closure", status=FAIL, threshold=1.0,
            message="unrepresentable interval residual", diagnostics={"states": list(states)},
        )

    blocks = []
    for label, start, stop in segments(w, params.get("segment_column", "_regime")):
        r = residual[start:stop]
        u = float((p[start:stop] + np.maximum(x[start:stop], 0)).sum())
        v = float((et[start:stop] + q[start:stop] + np.maximum(-x[start:stop], 0)).sum())
        denominator = max(u, v)
        allowance = max(threshold * denominator, floor)
        loss = float(np.maximum(r - deadband, 0).sum())
        creation = float(np.maximum(-r - deadband, 0).sum())
        if not np.isfinite([u, v, allowance, loss, creation]).all():
            return CriterionResult(
                name="non_cancelling_closure", status=FAIL, threshold=1.0,
                message=f"unrepresentable budget in block {label}",
            )
        blocks.append({
            "label": label, "start": start, "stop": stop,
            "supplied_mm": u, "exported_mm": v, "denominator_mm": denominator,
            "allowed_mm": allowance, "created_mm": creation, "lost_mm": loss,
            "signed_residual_mm": float(r.sum()), "absolute_residual_mm": float(np.abs(r).sum()),
            "max_interval_residual_mm": float(np.abs(r).max()),
            "allowance_ratio": max(loss, creation) / allowance,
        })
    if not blocks:
        return CriterionResult(
            name="non_cancelling_closure", status=FAIL, threshold=1.0,
            message="no scored wetting or drying blocks",
        )
    worst = max(blocks, key=lambda block: block["allowance_ratio"])
    failed = sum(block["allowance_ratio"] > 1 for block in blocks)
    return CriterionResult(
        name="non_cancelling_closure", status=FAIL if failed else PASS,
        value=worst["allowance_ratio"], threshold=1.0,
        message=(f"worst block {worst['label']}: created {worst['created_mm']:.6g} mm, "
                 f"lost {worst['lost_mm']:.6g} mm; allowed {worst['allowed_mm']:.6g} mm "
                 f"({failed}/{len(blocks)} blocks fail)"),
        diagnostics={"states": list(states), "step_deadband_mm": deadband,
                     "relative_tolerance": threshold, "absolute_tolerance_mm": floor,
                     "n_failed_blocks": failed, "worst_block": worst, "blocks": blocks},
    )
