(omega-dev-eos) =

## Equation of State

Omega includes an `Eos` class that provides functions that compute `SpecVol` and `SpecVolDisplaced`. 
Current EOS options are a linear EOS or an EOS computed using the TEOS-10 75 term expansion from 
Roquet et al. 2015. `SpecVol` is provided for either of these EOS options, while `SpecVolDisplaced` 
is only available for the TEOS-10 option. Trying to get `SpecVolDisplaced` while using the linear EOS
type will throw a run time error.

## Eos type

An enumeration listing all implemented schemes is provided. It needs to be extended every time an 
EOS is added. It is used to identify which EOS method is to be used at run time.
```c++
enum class EosType { Linear, Teos10Poly75t };
```

## Initialization

An instance of the `Eos` class requires a [`HorzMesh`](#omega-dev-horz-mesh), so the mesh class 
and all of its dependencies need to be initialized before the `Eos` class can be. The static method:
```c++
OMEGA::Eos::init();
```
initializes the default `Eos`. A pointer to it can be retrieved at any time using:
```c++
OMEGA::Eos* DefEos = OMEGA::Eos::getDefault();
```
The create method:
```c++ 
OMEGA::Eos::create(Name, Mesh, NVertLevels);
```
allocates the `SpecVol` and `SpecVolDisplaced` arrays for the mesh and vertical level dimensions.

## Computation of Eos

To compute `SpecVol` for a particular set of temperature, salinity, and pressure arrays, do 
```c++
Eos.computeSpecVol(SpecVol, ConsrvTemp, AbsSalinity, Pressure);
```
`SpecVolDisplaced` is calculated using local temperature and salinity values, but a pressure 
value at KDisp. To compute `SpecVolDisplaced` for a particular set of temperature, salinity, 
and pressure arrays and displaced vertical index level, do
```c++
Eos.computeSpecVolDisp(SpecVol, ConsrvTemp, AbsSalinity, Pressure, KDisp);
```
where KDisp is the vertical `k` index of the level you want to displace the specific volume to. 
For example, to displace the specific volume to the surface, set `KDisp = 0` (or `minLevelCell` 
when that becomes available).

## Removal of Eos

To erase a specific named Eos instance use `erase`
```c++
OMEGA::Eos::erase(Name);
```
To clear all instances do:
```c++
OMEGA::Eos::clear();
```