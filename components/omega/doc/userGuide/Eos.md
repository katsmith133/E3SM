(omega-user-eos)=

## Eos

The `Eos` class provides a container for the equation of state in OMEGA. Upon creation of an
`Eos` instance, functors are initialized and arrays for the specific volume (`SpecVol`) and 
displaced specific volume (`SpecVolDisp`) are allocated.

The user-configurable options are: `EosType` (either `Linear` or `Teos-10`), as well as the parameters needed for the linear equation of state.
```yaml
Eos:
   EosType : Linear
   Linear:
      DRhoDT: -0.2
      DRhoDS : 0.8
      RhoT0S0 : 1000.0
```
The `Linear` `EosType` computes specific volume using the linear expansion coefficients specified in the YAML file, and `Teos-10` `EosType` computes specific volume using the TEOS-10 75 term expansion from 
[Roquet et al. 2015](https://www.sciencedirect.com/science/article/pii/S1463500315000566).

Note: `SpecVolDisp` can not be calculated if the `Linear` `EosType` is selected, it will produce a run time error.