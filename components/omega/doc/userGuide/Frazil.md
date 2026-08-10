(omega-user-frazil)=

# Frazil

This page describes user-facing configuration for frazil tendencies in Omega.
Frazil physics represents the formation and melt of frazil ice within
the ocean water column. It impacts local layer pseudo-thickness (i.e. mass),
temperature, and salinity tendencies. The vertical sum of the frazil energy,
mass of water and mass of salt are passed to the coupler (if coupled) or
discarded (in ocean standalone mode).

## Configuration overview

Frazil behavior is controlled by one enable switch in Tendencies and one
Frazil configuration block:

```yaml
Omega:
  Tendencies:
    FrazilTendencyEnable: true

  Frazil:
    FrazilType: basic
    MassLimit: 0.1
    DepthLimit: -1.0
    ConservationCheck: false
```

- `Tendencies.FrazilTendencyEnable`
  - Enables/disables application of frazil tendency contributions.
- `Frazil.FrazilType`
  - Selects frazil option.
  - Current supported option: `basic`.
- `Frazil.MassLimit`
  - limits per-layer frazil mass/thickness tendency magnitude (applied
   to formation and melt)
- `Frazil.DepthLimit`
  - Limits depth range where frazil is computed.
  - Negative values mean no depth limit.
- `Frazil.ConservationCheck`
  - Enables a column-level diagnostic conservation check with logging.

## Available frazil option

Omega currently supports one active frazil pathway:

- `basic`
  - freezing is based on the formation of fresh solid ice, to which
  salt is added (similar the mpas-ocean implementation).

This pathway contributes to:

- pseudo-thickness tendency
- temperature tracer tendency
- salinity tracer tendency

## Notes

- Frazil tendencies are applied through the `Tendencies` tracer-step workflow.
- The frazil tendency hook assumes tracer names include `Temperature` and
  `Salinity`.
- For implementation and algorithm details, see
  [Developer Frazil Guide](../devGuide/Frazil.md).
