"""Warm wetting, drying, rewetting and recovery blocks for interval closure."""

from __future__ import annotations

import numpy as np
import pandas as pd

SPINUP_DAYS = 365
PERIOD_DAYS = 365
STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 40.0,
}


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    n = SPINUP_DAYS + PERIOD_DAYS
    rain = np.zeros(n)
    pet = np.empty(n)
    regimes = np.full(n, "spinup", dtype=object)

    # Familiar warm weather fills the stores without a snow confounder.
    # Keep spinup dry so the deliberately broken reference's first two wet
    # rows are unambiguously inside the scored record.
    rain[:SPINUP_DAYS] = 0.0
    pet[:SPINUP_DAYS] = rng.uniform(1.2, 3.0, SPINUP_DAYS)

    # Every wetting block begins with two consecutive wet days. This is the
    # witness that lets a one-day storage offset cancel by the block end.
    # Keep each wetting pulse short enough that the 4 mm witness remains above
    # the block's 5% allowance while still giving the model long dry/recovery
    # stretches in which antecedent storage supplies ET and runoff.
    pattern = [("wetting", 2), ("drying", 28), ("rewetting", 2), ("recovery", 11)]
    pos = SPINUP_DAYS
    while pos < n:
        for label, duration in pattern:
            if pos >= n:
                break
            stop = min(pos + duration, n)
            regimes[pos:stop] = label
            length = stop - pos
            if label in {"wetting", "rewetting"}:
                rain[pos:stop] = rng.uniform(35.0, 55.0, length)
                pet[pos:stop] = rng.uniform(0.5, 1.5, length)
            else:
                pet[pos:stop] = rng.uniform(1.0, 2.5, length)
            pos = stop

    time = pd.date_range("2000-01-01", periods=n, freq="D")
    forcing = pd.DataFrame({
        "time": time.strftime("%Y-%m-%d"),
        "pr": np.round(rain, 6),
        "tas": np.full(n, 15.0),
        "pet": np.round(pet, 6),
        "_regime": regimes,
    })
    return forcing, dict(STATIC)
