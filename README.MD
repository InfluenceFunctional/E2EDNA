# E2EDNA - Tinker Implementation

Michael Kilgour, Tao Liu, Lena Simine
mjakilgour gmail com

* If using our code, please consult and cite the following publication: **J. Chem. Inf. Model. 2021, 61, 9, 4139–4144** (https://doi.org/10.1021/acs.jcim.1c00696)

* For an updated version (E2EDNA 2.0) implemented using OpenMM, please go to: https://github.com/siminegroup/E2EDNA2

## INSTALLATION
* Get the repository
```
git clone git@github.com:InfluenceFunctional/e2edna.git WHAT_YOU_WANT_IT_TO_BE_CALLED
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
* Set paths in main.py. Note, there is a toggle for running on a 'local' machine vs 'cluster', with distinct paths.
```
* params['minimize', 'dynamic', 'pdbxyz', 'xyzedit', 'archive'] --> all set to the Tinker executables on your system
* params['mmb'] --> set to the MacroMoleculeBuilder executable on your system
* params['workdir'] --> the working directory created in the previous step
* params['analyte xyz'] --> path to the analyte structure in tinker .xyz format
* various .key files --> parameters may be specific to your analyte if custom force-field parameters are required. Adjust keyfiles accordingly, leaving the keywords on top as-is.
```
	
## RUNNING A JOB
Before running a job, ensure installation, analyte parameterization and simulation parameters have been set.
1. Set 'run num' in main.py to zero for a fresh run, which will create a new run directory in your workdir. If you want to pick-up on a previous run, put its number here. 
   - Using argparse, one can directly set this via command line from a submission script rather than editing main.py itself.
2. Set simulation parameters in main.py, and select the mode to run in as params['simulation type']. 'free aptamer' runs the aptamer by itself, 'binding' runs a representative structure of the aptamer complexed with the analyte, and 'analysis' redoes reaction coordinate analysis given pre-existing trajectories.
4. Set the DNA sequence as a FASTA string to the 'sequence' variable at the bottom of main.py.
5. Run the code, e.g., for run_num = 0
   ```
   python main.py --run_num=0
   ```
   Inputs and outputs will go to /workdir/run0. In this directory, /outfiles contains outputs for debugging the various processes within the script.

   In general, one would run first a 'free aptamer' simulation to retrieve a representaitive 3D structure and confirm the secondary structure stability. 
   Then, one would change the operating mode to 'binding' to assess aptamer-analyte binding.
6. Simulation parameters and outputs including energy traces are saved to 'e2ednaOutputs.npy', which can be reloaded as
   - outputs = np.load('e2ednaOutputs.npy',allow_pickle=True).item()
