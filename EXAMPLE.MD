# E2EDNA - Tinker Implementation 

## EXAMPLE
Here is an example of executing this code on a minimal system.
### System
For speed and simplicity, we will look at a tiny hairpin, GCTTTGC, and asses it's binding with UTP(-4) in PBS at 310K.

### Setup
We will use our default settings, and track the distances between bases which are expected to pair (1-7 and 2-6) as our reaction coordinates. 
The specific atoms to follow in these bases must be determined manually by users.
```
# Simulation parameters
params['reaction coordinates'] = [[26,213],[53,185]]
params['force field'] = 'AMOEBA'
params['target'] = 'UTP-4' # 
params['minimization gradrms'] = 0.05
params['equilibration time'] = 0.5 # equilibration time in nanoseconds
params['sampling time'] = 10 # sampling time in nanoseconds
params['time step'] = 2.0 # in fs
params['num charges'] = 4 # number of positive counterions (Na+) to add to simulation box to counter the analyte.
params['NaCl concentration'] = 163 # in mmol
params['box offset'] = 5 # in angstroms 
params['print step'] = 1 # printout step in ps 
params['heavy hydrogen'] = False 
params['outisde secondary structure'] = False 

# analyte placement - only relevant if 'simulation type' is 'binding'
params['analyte position'] = 'random' # 'random' - will be at least 2 nm from the aptamer in one of 8 cubic face directions, 'manual' places the aptamer center at the given coordinates relative to origin

# AMOEBA parameters
params['vdW radius'] = 12 # van der Waals radius
params['polarization version'] = 'OPT4' # OPT3 or OPT4 - 4 is slower but more accurate
params['polar eps'] = 1e-5 # induced dipole convergence parameter default 1e-6
params['polar predict'] = True # whether to run with polar-predict or not
```

### Free Aptamer Run

This first run returns the trajectory of the specified reaction coordinates above, in this case corresponding to the expected hydrogen bonding bases.
We compute the free energy of the system along these coordinates by histogramming the reaction coordinate trajectories, and taking their negative logarithm.

_Reaction Coordinates Histograms and Free Energies_

This is a trivial example, and we can see that, as expected, the paired structure is stable. 

This run will also return a 3D structure, labelled 'representativeStructure.xyz' which lies in the minima of all specified reaction coordinate free energy surfaces.
This will be the initial aptamer configuration for an eventual binding simulation. 

### Binding evaluation

_Representative Aptamer Structure Complexed with UTP-4 Analyte_

_Reaction Coordinate Histograms and Free Energies