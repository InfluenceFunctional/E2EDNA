# E2EDNA - Tinker Implementation

Michael Kilgour, Tao Liu, Lena Simine
mjakilgour gmail com

## INSTALLATION
* Get the repository
```
it clone git@github.com:InfluenceFunctional/e2edna.git WHAT_YOU_WANT_IT_TO_BE_CALLED
```
* Setup python environment with the following packages
```
biopython==1.78
matplotlib==3.3.2
MDAnalysis==1.0.0
numpy==1.19.2
pandas==1.1.3
ParmEd==3.2.0
Pillow==7.2.0
pyparsing==2.4.7
scipy==1.4.1
seqfold==0.7.7
```
* Install appropriate MacroMoleculeBuilder for your system from [SimTK here](https://simtk.org/projects/rnatoolbox) - old and new versions have been tested and seem to work
* Install Tinker from [here](https://github.com/tinkertools) or [here](https://dasher.wustl.edu/tinker/). Tinker9 with GPU support is recommended, but somewhat more involved to install. Tinker8 binaries can be easily downloaded.
* Create working directory where runs will take place
* Set paths in main.py
```
* minimize, dynamic, pdbxyz, xyzedit, archive --> all set to the Tinker executables
* mmb --> set to the MacroMoleculeBuilder executable
* workdir --> the working directory created in the previous step
* analyte --> the analyte structure in tinker .xyz format
* various .key files --> parameters may be specific to your analyte if custom force-field parameters are required. Adjust keyfiles accordingly, leaving the keywords on top as-is.
```
	
## RUNNING A JOB
* set 'run num' to zero for a fresh run, or if you want to pick-up on a previous run, put its number here 
  - using argparse, one can directly set this via command line from sub.sh rather than editing aptamerBinding.py itself
* Simulation parameters
  - note that AMOEBA and AMBER are possible, though AMBER is a work-in-progress in some ways
  - minimization gradrms - the tightness of the minimization routine - 0.1 is good. 1 may be OK but it's on the high side
  - equilibration time - in ns
  - samling time - in ns
  - time step - in fs. If we 1-2ish without heavy-hydrogen. up to 3 or maybe 3.5 with heavy-hydrogen 
* adjust Keyfiles: All of the other tinker parameters are set in the keyfiles. 
  - there are also .in files which are basically tinker scripts, and .dat files which define the MMB folding protocol
* set the DNA sequence as a string to the 'sequence' variable at the bottom of the file. It will use seqfold and MMB to fold it (note, seqfold is pretty conservative, and will generally only recommend high-probability secondary structure, potentially leaving long structureless strings. Also I'm not 100% certain it's pairList function is working, but I have a better one I can sub in when I find
time to check it.
* the script will send tinker output files under the workdir/outfiles directory
* simulation parameters and outputs including energy traces are saved to 'bindingOutputs.npy', which can be reloaded as
  - outputs = np.load('bindingOutputs.npy',allow_pickle=True)
  - outputs = outputs.item()
* Visualizing outputs
  - note that both VMD and pymol give crazy results if you have the periodic box information left in your .xyz or .arc files. For .xyz you can remove manually. For .arc trajectories, there is a function in binder which removes all periodic information.
  - one may also plot the energy and analyte-sequence center-of-mass distance with a commented script directly after the output save command

