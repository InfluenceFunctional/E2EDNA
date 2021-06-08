### README FOR APTAMERBINDING CODE ###
Michael Kilgour
mjakilgour@gmail.com

INSTALLATION
==> if it's not already, connect your github account https://docs.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent
==> in your desired directory from which to run (I think this should work! Not an expert but it's also not crazy complicated so).
	: git clone git@github.com:InfluenceFunctional/tinkertest.git WHAT_YOU_WANT_IT_TO_BE_CALLED
==> in the aptamerBinding.py script, after elif params['device'] == 'cluster': there are a number of paths to set. 
	>> set all tinker8 and tinker9 (GPU-only) paths to their appropriate places (for Simine group members, this is already done)
	>> set macromolecule builder (MMB) path (for Simine group members this is already done)
	>> choose a workdir - I suggest making a new directory in your personal scratch, all runs and outpts will be placed here, with a new run directory made for each simulation
	>> Python environment: 
		-> I like to set the python version first - shouldn't matter a lot but I'm a little superstitious 
			: module load python/3.7
		-> make and activate a new python environment 
			: virtualenv --no-download ~/YOUR-ENV-NAME
			: source ~/YOUR-ENV-NAME/bin/activate
		-> update pip
			: pip install -no-index --upgrade pip
		-> install the requisite packages. Compute Canada has their own wheels for certain programs, so we add --no-index so we get the optimal CC version. Some of these are dependencies for MDAnalysis so they may be redundant.
			: pip install numpy --no-index
			: pip install MDAnalysis --no-index
			: pip install matplotlib
			: pip install tqdm
			: pip install pandas --no-index
			: pip install Pillow --no-index
			: pip install scipy --no-index
			: pip install seqfold
			: pip install parmed
			

==> Other paths: 
	>> set: analyte xyz (xyz for your binding analyte)
		-> for UTP-4 - this is pre-loaded
		-> UTP-3 is deprecated (bad parameters)
	>> in your workdir, must have a directory '/params' which contains the relevant FF parameters - most importantly amoebabio18.prm
	>> all others should be automatic
	
==> Submission setup:
	>> in the directory which contains aptamerBinding.py etc. 
		-> if you have modified any files, delete them
		-> get the newest version of the code
			: git pull
	>> sbatch sub.sh will run the code
	>> ensure settings in sub.sh are appropriate (email, memory, walltime). Beluga only has v100 GPUs, whereas Cedar has v100l's. This must be represented in the resource request.
	>> set the path to your virtual environment for this project
	
	
RUNNING A JOB
	>> set 'run num' to zero for a fresh run, or if you want to pick-up on a previous run, put its number here
		-> using argparse, one can directly set this via command line from sub.sh rather than editing aptamerBinding.py itself
	>> Simulation parameters
		-> note that AMOEBA and AMBER are possible, though AMBER is a work-in-progress in some ways
		-> minimization gradrms - the tightness of the minimization routine - 0.1 is good. 1 may be OK but it's on the high side
		-> equilibration time - in ns
		-> samling time - in ns
		-> time step - in fs. If we 1-2ish without heavy-hydrogen. up to 3 or maybe 3.5 with heavy-hydrogen
	>> adjust Keyfiles: All of the other tinker parameters are set in the keyfiles. 
		-> there are also .in files which are basically tinker scripts, and .dat files which define the MMB folding protocol
	>> set the DNA sequence as a string to the 'sequence' variable at the bottom of the file. It will use seqfold and MMB to fold it (note, seqfold is pretty conservative, and will generally only recommend high-probability secondary structure, potentially leaving long structureless strings. Also I'm not 100% certain it's pairList function is working, but I have a better one I can sub in when I find
time to check it.
	>> the script will send tinker output files under the workdir/outfiles directory
	>> simulation parameters and outputs including energy traces are saved to 'bindingOutputs.npy', which can be reloaded as
		: outputs = np.load('bindingOutputs.npy',allow_pickle=True)
		: outputs = outputs.item() 
	>> Visualizing outputs
		-> note that both VMD and pymol give crazy results if you have the periodic box information left in your .xyz or .arc files. For .xyz you can remove manually. For .arc trajectories, there is a function in binder which removes all periodic information.
		-> one may also plot the energy and analyte-sequence center-of-mass distance with a commented script directly after the output save command


Some notes from LS:
1. the name of the workdir  in aptamerBinding.py must end with 'run' , I tried 'runs' and the job failed; it can start with anything: testrun, nextrun, etc.. are ok
