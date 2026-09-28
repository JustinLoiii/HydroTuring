# Non-cancelling water conservation through wetting-drying transitions

## The idea

This probe tests whether the reported water fluxes and reported water stores
describe the same finite-volume update at each output interval. A signed error
can disappear when a model creates water in one wet interval and removes the
same amount in the next interval. The probe therefore scores the interval
residuals separately inside short wetting, drying, rewetting and recovery
blocks.

For each interval, with fluxes in mm/day and `dt` in days, the residual is

```text
r = (pr + gwex - evspsbl - mrro) * dt - (S_end - S_start)
```

`gwex` is included only when the model reports it. `S` includes every water
state the model reports, including `mrso`, `gw`, `canopy`, `snw` and
`channel`. `dis`, snowmelt and internal groundwater-surface-water exchange
are not additional external terms.

Within each contiguous `_regime` block, positive and negative residuals are
accumulated independently after a 0.001 mm per-interval output deadband. Each
total must be no larger than `max(0.05 * D, 0.01 mm)`, where `D` is the larger
of the block's supplied and exported water throughput.

## Case construction

The deterministic generator uses a dry spinup followed by repeated warm,
non-freezing blocks:

```text
wetting → drying → rewetting → recovery
```

The two-day wetting blocks contain consecutive rainfall, while the dry blocks
test storage-supported evaporation and drainage. Rainfall intensity, PET and
the number of repeated blocks vary with the seed. No data files are committed.

## Gate references

`reference_bucket`, `flex_lumped`, `flex_topo` and `sacsma_snow17` are required
to pass. `reference_event_storage_reset` adds 4 mm to `mrso` after the first
wet interval and removes it on the next contiguous wet interval, without
changing fluxes or the final state. Whole-event signed closure cancels this
error, while `non_cancelling_closure` fails it.

## Scope and limitation

The probe catches temporary errors visible in the aggregate reported stores.
Errors that cancel entirely inside an unreported internal layer cannot be
detected without additional model fields. It does not certify root-zone
physiology, layer-level Richards variables or equality of hydrographs between
different model structures.

## Reproduce locally

```bash
ht gate --probe mass/non-cancelling-wet-dry-closure
```
