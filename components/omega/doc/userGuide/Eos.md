(omega-user-eos)=

## Eos

The `Eos` class provides a container for the equation of state in OMEGA, and member functors that are the linear and TEOS-10 equations of state.
Upon creation of an `Eos` instance, these functors are initialized and the array for the specific volume is allocated.
The user-configurable options are: `EosType` (either `Linear` or `teos-10`), as well as the parameters needed for the linear equation of state.

```yaml
Eos:
   EosType : Linear
   Linear:
      DRhoDT: -0.2
      DRhoDS : 0.8
      RhoT0S0 : 1000.0
```
