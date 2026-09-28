import pandas as pd

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.harness import build_case, load_generator
from hydroturing.runner import get_runner
from hydroturing.seeds import gate_seeds


PROBE_ID = "mass/non-cancelling-wet-dry-closure"


def test_generator_is_deterministic_and_labels_contiguous_blocks():
    probe = registry.find_probe(PROBE_ID)
    generator = load_generator(probe)
    first, static_first = generator.generate(20260903)
    second, static_second = generator.generate(20260903)
    pd.testing.assert_frame_equal(first, second)
    assert static_first == static_second
    assert first["_regime"].iloc[0] == "spinup"
    assert set(first["_regime"]) == {"spinup", "wetting", "drying", "rewetting", "recovery"}


def test_reference_models_pass_and_storage_reset_fails(tmp_path):
    probe = registry.find_probe(PROBE_ID)
    case = build_case(probe, gate_seeds(PROBE_ID, 5)[0])
    for name in ("reference_bucket", "flex_lumped", "flex_topo", "sacsma_snow17"):
        model = registry.find_model(name)
        run = get_runner(model).run(model, probe, case, tmp_path / name)
        result = get("non_cancelling_closure")(run, probe, {
            "threshold": 0.05,
            "absolute_tolerance_mm": 0.01,
            "step_deadband_mm": 0.001,
            "segment_column": "_regime",
        })
        assert result.passed, f"{name}: {result.message}"

    model = registry.find_model("reference_event_storage_reset")
    run = get_runner(model).run(model, probe, case, tmp_path / model.name)
    result = get("non_cancelling_closure")(run, probe, {
        "threshold": 0.05,
        "absolute_tolerance_mm": 0.01,
        "step_deadband_mm": 0.001,
        "segment_column": "_regime",
    })
    assert not result.passed
    assert any(block["created_mm"] > 0 for block in result.diagnostics["blocks"])
    assert any(block["lost_mm"] > 0 for block in result.diagnostics["blocks"])
