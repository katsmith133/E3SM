(omega-dev-frazil)=

# Frazil

This page describes frazil design and implementation details in Omega,
for the current basic pathway. A TEOS-10 option will be added soon.

## Purpose and coupling points

Frazil computes phase-change-related tendencies that modify:

- pseudo-thickness tendency
- temperature tracer tendency
- salinity tracer tendency

The tendency hook-up is implemented through `FrazilOnCell` in the tracer
tendency phase, where frazil contributions are added to accumulated
`PseudoThicknessTend` and `TracerTend` arrays.

## Data flow and call sequence

1. `Tendencies.computeTracerTendenciesOnly` checks
   `Tendencies.FrazilTendencyEnable`.
2. If enabled, `FrazilOnCell.operator()` retrieves the default Frazil object
   and zeros frazil tendency and accumulator arrays.
3. `FrazilOnCell` extracts `Temperature` and `Salinity` tracer subviews, then calls
   `Frazil.computeFrazil(CT, SA, PressureMid, PseudoThickness)`.
4. `Frazil.computeFrazil` dispatches using FrazilType. The only current supported
   option is basic. A TEOS-10 option will be added soon.
5. Returned frazil tendencies are added into `PseudoThicknessTend` and
   `TracerTend` for active cell layers.

## Configuration coupling

Frazil behavior is configured with:

- `Omega.Tendencies.FrazilTendencyEnable`
  - Global switch for applying frazil tendency terms.
- `Omega.Frazil.FrazilType`
  - Implementation choice. Only current supported value: `basic`.
- `Omega.Frazil.MassLimit`
  - Per-layer mass and thickness limiter used by basic formation and melt.
- `Omega.Frazil.DepthLimit`
  - Optional depth cutoff for frazil activity. Negative means no cutoff.
- `Omega.Frazil.ConservationCheck`
  - Optional post-compute column conservation diagnostic logging.

## Basic pathway summary

- Freezing-point checks are based on conservative temperature and absolute
  salinity with pressure-dependent freezing temperature.
- Vertical accumulation order is bottom-to-top within each active column.
- Frazil tendencies are not time-step scaled inside `Frazil`; they are
  accumulated as tendency contributions.
- Surface salt redistribution is applied at the top active layer.
- Column accumulators are converted to coupler units at the end of the
  per-column loop.

### Basic pathway (`FrazilType: basic`)

- Uses simplified energetics: the energy of the super-cooled water sets the
amount of solid ice formed (used constant latent heat of fusion of fresh ice).
Salt is added based on a constant bulk salinity `IceRefSal` (default) or a
manual toggle (for now) using the local salinity and the frazil porosity.
Melting of existing frazil is set by the amount of pure ice that can be melted
by the warm layer, and the enthalpy of melt water at the local freezing point
is added to the temperature tendency.
- Computes local layer tendencies (`HTend`, `TTend`, `STend`) and updates
  accumulated frazil stores.
- Applies surface salt redistribution adjustment and converts accumulators to
  coupler units at the end of the column loop.
- Warnings: 1) basic frazil formation/melt is does not conserve energy.
  2) Using porosity to set the salt content includes a redistribution of excess
  salt in the surface layer, which can be very significant.

## Existing ctest coverage

The existing frazil test driver is in
`components/omega/test/ocn/FrazilTest.cpp` and covers:

- basic frazil formation in cold and warm single-layer states
- mixed warm/cold column behavior with sign checks for branch switching
- depth-limit behavior ensuring excluded deep layers have zero frazil tendency

## Extensibility

The frazil choice dispatch and configuration plumbing are retained so future
frazil options can be added without changing Tendencies call sites.

## Related pages

- User-facing options: [User Frazil Guide](../userGuide/Frazil.md)
- Tendency container and hook-up: [Tendencies](Tendencies.md)
