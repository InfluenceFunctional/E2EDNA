# import statements
from shutil import copyfile
import re
from seqfold import dg, fold
from utils import *


'''
script which accepts a DNA sequence and an analyte, and returns the binding affinity

To-Do:
==> implement AMBER mode
    -> testing
    -> fix analysis script
    -> fix charge of peptides
==> check pairlist
==> implement nupack with saltiness
==> add docker
==> impelement docker recommendations
==> allow for more than 1000 frames and/or kill last frame grabber
'''


params = {}
params['device'] = 'cluster' # 'local' or 'cluster'
params['explicit run enumeration'] = False # if this is True, the next run be fresh, in directory 'run%d'%run_num, if false, regular behaviour. Note: only use this on fresh runs

if params['device'] == 'cluster':
    params['run num'] = get_input()
elif params['device'] == 'local':
    params['run num'] = 0 # manual setting, for 0, do a fresh run, for != 0, pickup on a previous run.

# Simulation parameters
params['force field'] = 'AMOEBA' # AMOEBA or AMBER (WiP)
params['target'] = 'UTP-4' # 'UTP-4' setup only for AMOEBA, 'peptide' setup only for AMBER
params['minimization gradrms'] = 0.05 # flag that minimization is converged 0.1-0.01 is an ok range usually
params['equilibration time'] = 0.1 # equilibration time in nanoseconds # don't make this or sampling time divisible by 2xdt
params['sampling time'] = 10 # sampling time in nanoseconds
params['time step'] = 2.0 # in fs
params['num charges'] = 4 # number of positive charges (Na+) to add to simulation box
params['NaCl concentration'] = 163 # concentration of NaCl in mmol ~ approximating for other heteroions (PBS)
params['box offset'] = 5 # in angstroms box size will be order of longest edge + 2 * vdW radius + 2 * this offset
params['print step'] = 10 # printout step in ps - if using 'getFinalFrame' you want there to be more than 2 and less than 1000 total frames in any trajectory
params['heavy hydrogen'] = False # True or False - yes means hydrogens have extra mass, allowing longer time-steps up to 3-4 fs
params['outisde secondary structure'] = False # use a secondary structure manually coded as a list in a relevant file

# AMBER params
#???

# Tinker parameters
params['vdW radius'] = 12 # van der Waals radius
params['polarization version'] = 'OPT4' # OPT3 or OPT4 - 4 is slower but more accurate
params['polar eps'] = 0.00001 # induced dipole convergence parameter default 1e-6
params['polar predict'] = True # whether to run with polar-predict or not


# paths
if params['device'] == 'local':
    params['minimize path'] = 'minimize'
    params['dynamic path'] = 'dynamic'
    params['pdbxyz path'] = 'pdbxyz'
    params['xyzedit path'] = 'xyzedit'
    params['archive path'] = 'archive'
    params['workdir'] = 'C:/Users\mikem\Desktop/tinkerruns'
    params['mmb'] = 'C:/Users/mikem/Desktop/Installer.2_14.Windows/MMB.2_14.exe'
    params['mmb params'] = 'lib/MMB/parameters.csv'
    params['mmb template'] = 'lib/MMB/commands.template.dat'
elif params['device'] == 'cluster':
    params['minimize path'] = 'sh ~/projects/def-simine/programs/tinker9/bin/cc70/gpu-m/minimize9.sh'
    params['dynamic path'] = 'sh ~/projects/def-simine/programs/tinker9/bin/cc70/gpu-m/dynamic9.sh'
    params['pdbxyz path'] = '~/projects/def-simine/programs/tinker8/pdbxyz'
    params['xyzedit path'] = '~/projects/def-simine/programs/tinker8/xyzedit'
    params['archive path'] = '~/projects/def-simine/programs/tinker8/archive'
    params['workdir'] = '/home/kilgourm/scratch/simruns' # specify your working directory here
    params['mmb'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64/MMB.2_14.Linux64'
    params['mmb params'] = 'lib/MMB/parameters.csv'
    params['mmb template'] = 'lib/MMB/commands.template.dat'


# find analyte-specific files
if params['target'] == 'UTP-4':
    params['analyte xyz'] = 'lib/UTP-4/UTP-4.xyz'

    # waters
    params['waterbox'] = 'lib/water/waterbox.in'
    params['water'] = 'lib/water/water.xyz'

    # copy Tinker keyfiles
    params['min key'] = 'lib/keyfiles/UTP-4/minimize.key'
    params['equil key'] = 'lib/keyfiles/UTP-4/equilibrate.key'
    params['dyn key'] = 'lib/keyfiles/UTP-4/dynamics.key'
    params['combined params key'] = 'lib/keyfiles/UTP-4/params_combined.key'

    # copy tinker utilities
    params['grablastframe'] = 'lib/infiles/grablastframe.in'
    params['movesomething'] = 'lib/infiles/movesomething.in'
    params['killWater'] = 'lib/infiles/killWater.in'
    params['addIons'] = 'lib/infiles/addIons.in'
    params['addSodium'] = 'lib/infiles/addSodium.in'
    params['addChloride'] = 'lib/infiles/addChloride.in'
    params['origin'] = 'lib/infiles/origin.in'
    params['mmb template'] = 'lib/commands.template.dat'
elif params['target'] == 'peptide':
    params['analyte xyz'] = 'lib/peptide/peptide.xyz'

    # waters
    params['waterbox'] = 'lib/water/waterbox2.in'
    params['water'] = 'lib/water/water2.xyz'

    # copy Tinker keyfiles
    params['min key'] = 'lib/keyfiles/peptide/minimize.key'
    params['equil key'] = 'lib/keyfiles/peptide/equilibrate.key'
    params['dyn key'] = 'lib/keyfiles/peptide/dynamics.key'
    params['combined params key'] = 'lib/keyfiles/peptide/params_combined.key'

    # copy tinker utilities
    params['grablastframe'] = 'lib/infiles/grablastframe.in'
    params['movesomething'] = 'lib/infiles/movesomething.in'
    params['killWater'] = 'lib/infiles/killWater2.in'
    params['addIons'] = 'lib/infiles/addIons2.in'
    params['addSodium'] = 'lib/infiles/addSodium2.in'
    params['addChloride'] = 'lib/infiles/addChloride2.in'
    params['origin'] = 'lib/infiles/origin.in'
    params['mmb template'] = 'lib/commands.template.dat'


class binder():
    def __init__(self,sequence,params):
        self.params = params
        self.sequence = sequence

        if self.params['target'] == 'UTP-4':
            self.MM_params = '../params/amoebabio18.prm'
        elif self.params['target'] == 'peptide':
            self.MM_params = '../params/amber99.prm'

        self.setup()
        self.getCheckpoint()


    def setup(self):
        '''
        setup working directory
        copy in relevant xyz and keyfiles
        move to relevant directory
        :return:
        '''
        if (self.params['explicit run enumeration'] == True) or (self.params['run num'] == 0):
            if self.params['run num'] == 0:
                self.makeNewWorkingDirectory() # make new workdir in sequence
            else:
                self.workDir = self.params['workdir'] + '/run%d'%self.params['run num'] # make a new workdir with the given index
                os.mkdir(self.workDir)

            os.mkdir(self.workDir + '/outfiles')
            # copy structure files
            copyfile(self.params['analyte xyz'], self.workDir + '/analyte.xyz')
            copyfile(self.params['waterbox'], self.workDir + '/waterbox.in')
            copyfile(self.params['water'], self.workDir + '/water.xyz')

            # copy relevant keyfiles
            copyfile(self.params['min key'], self.workDir + '/minimize.key')
            copyfile(self.params['equil key'], self.workDir + '/equilibrate.key')
            copyfile(self.params['dyn key'], self.workDir + '/dynamics.key')
            copyfile(self.params['combined params key'], self.workDir + '/params_combined.key')

            # append keyfile parameters
            replaceText(self.workDir + '/minimize.key', "VDWCUTOFF", str(self.params['vdW radius']))
            replaceText(self.workDir + '/dynamics.key', "VDWCUTOFF", str(self.params['vdW radius']))
            replaceText(self.workDir + '/equilibrate.key', "VDWCUTOFF", str(self.params['vdW radius']))

            replaceText(self.workDir + '/minimize.key', "OPTVER", str(self.params['polarization version']))
            replaceText(self.workDir + '/dynamics.key', "OPTVER", str(self.params['polarization version']))
            replaceText(self.workDir + '/equilibrate.key', "OPTVER", str(self.params['polarization version']))

            replaceText(self.workDir + '/minimize.key', "POLAR-EPS", str(self.params['polar eps']))
            replaceText(self.workDir + '/dynamics.key', "POLAR-EPS", str(self.params['polar eps']))
            replaceText(self.workDir + '/equilibrate.key', "POLAR-EPS", str(self.params['polar eps']))

            if self.params['heavy hydrogen'] == False:
                replaceText(self.workDir + '/minimize.key', "heavy-hydrogen", "#heavy-hydrogen")
                replaceText(self.workDir + '/dynamics.key', "heavy-hydrogen", "#heavy-hydrogen")
                replaceText(self.workDir + '/equilibrate.key', "heavy-hydrogen", "#heavy-hydrogen")

            if self.params['polar predict'] == False:
                replaceText(self.workDir + '/minimize.key', "polar-predict", "#polar-predict")
                replaceText(self.workDir + '/dynamics.key', "polar-predict", "#polar-predict")
                replaceText(self.workDir + '/equilibrate.key', "polar-predict", "#polar-predict")

            # copy MMB params
            copyfile(self.params['mmb params'], self.workDir + '/parameters.csv')

            #copy tinker utilities
            copyfile(params['grablastframe'], self.workDir + '/grablastframe.in')
            copyfile(params['movesomething'], self.workDir + '/movesomething.in')
            copyfile(params['killWater'], self.workDir + '/killWater.in')
            copyfile(params['addIons'], self.workDir + '/addIons.in') # not necessary anymore
            copyfile(params['addSodium'], self.workDir + '/addSodium.in')
            copyfile(params['addChloride'], self.workDir + '/addChloride.in')
            copyfile(params['origin'], self.workDir + '/origin.in')
            copyfile(params['mmb template'], self.workDir + '/commands.template.dat')

            if params['outisde secondary structure']:
                copyfile('fullPairList.npy', self.workDir + '/fullPairList.npy')


        else:
            self.workDir = self.params['workdir'] + '/' + 'run%d' %self.params['run num']


        # move to working dir
        os.chdir(self.workDir)

        # save params to outputs
        outputs = {}
        outputs['params'] = params
        np.save('bindingOutputs', outputs)  # unpack with load then outputs.item()


    def makeNewWorkingDirectory(self):    # make working directory
        '''
        make a new working directory
        non-overlapping previous entries
        :return:
        '''
        workdirs = glob.glob(self.params['workdir'] + '/' + 'run*') # check for prior working directories
        if len(workdirs) > 0:
            prev_runs = []
            for i in range(len(workdirs)):
                prev_runs.append(int(workdirs[i].split('run')[-1]))

            prev_max = max(prev_runs)
            self.workDir = self.params['workdir'] + '/' + 'run%d' %(prev_max + 1)
            os.mkdir(self.workDir)
            print('Starting Fresh Run %d' %(prev_max + 1))
        else:
            self.workDir = self.params['workdir'] + '/' + 'run1'
            os.mkdir(self.workDir)


    def getCheckpoint(self):
        '''
        identify if any work has been done in this directory, and if so, where to pick up
        :return:
        '''
        try:
            f = open('checkpoint.txt','r')#if it does, see how long it is
            text = f.read()
            f.close()
            self.checkpoints = len([m.start() for m in re.finditer('\n', text)])

            if self.checkpoints == 1:
                self.checkpoint = 'initialized'
            elif self.checkpoints == 2:
                self.checkpoint = 'something else'
                '''
                and so on
                '''

            print("Resuming Run #%d"%int(self.params['run num']))


        except:
            f = open('checkpoint.txt','w')
            f.write('Started!\n')
            f.close()
            self.checkpoints = 1
            self.checkpoint = 'initialized'


    def run(self):
        '''
        run the binding simulation end-to-end
        consult checkpoints to not repeat prior outputs
        :return:
        '''
        if self.checkpoints < 2: # if we have initialized
            # get secondary structure
            print("Get Secondary Structure")
            if os.path.exists("pre_fold.pdb"): # if we have a pre-folded structure, do nothing
                pass
            elif os.path.exists('fullPairList.npy'):
                self.pairList = np.asarray(np.load('fullPairList.npy',allow_pickle=True).item()['list'])
            else:
                self.ssString, self.pairList = self.getSecondaryStructure(sequence)
            writeCheckpoint("Got Secondary Structure")


        if self.checkpoints < 3: # if we have the secondary structure
            # fold it!
            print("Folding Sequence")
            if os.path.exists("pre_fold.pdb"):
                copyfile('pre_fold.pdb','sequence.pdb')
                #os.system(self.params['pdbxyz path'] + ' sequence.pdb -k params_combined.key > outfiles/pdbconversion.out')
            else:
                #self.foldSequence(sequence, self.pairList)
                pass
            writeCheckpoint("Folded Sequence")


        if self.checkpoints < 4: # if we've folded the sequence
            # add the analyte!
            print('Add Analyte')
            #self.addAnalyte('sequence.xyz', 'analyte.xyz','params_combined.key')
            #os.rename('sequence.xyz','complex.xyz')
            writeCheckpoint('Added Analyte')


        if self.checkpoints < 5: # if we added the analyte
            print('Soak Complex')
            copyfile('complex.xyz','complex_pre_soak.xyz')
            complex_xrange, complex_yrange, complex_zrange = self.soak('complex.xyz','params_combined.key')
            replaceText('minimize.key','XX','%.1f'%complex_xrange)
            replaceText('minimize.key','YY','%.1f'%complex_yrange)
            replaceText('minimize.key','ZZ','%.1f'%complex_zrange)
            replaceText('equilibrate.key','XX','%.1f'%complex_xrange)
            replaceText('equilibrate.key','YY','%.1f'%complex_yrange)
            replaceText('equilibrate.key','ZZ','%.1f'%complex_zrange)
            replaceText('dynamics.key','XX','%.1f'%complex_xrange)
            replaceText('dynamics.key','YY','%.1f'%complex_yrange)
            replaceText('dynamics.key','ZZ','%.1f'%complex_zrange)
            os.rename('complex.xyz','complex_soaked.xyz')
            writeCheckpoint('Complex Soaked')


        if self.checkpoints < 6: # if we have soaked it
            # neutralize the analyte
            print('Neutralize Complex')
            copyfile('complex_soaked.xyz','complex_to_neutralize.xyz')
            self.neutralize('complex_to_neutralize.xyz','params_combined.key')
            os.rename('complex_to_neutralize.xyz','complex_neutralized.xyz')
            writeCheckpoint('Complex Neutralized')


        if self.checkpoints < 7: # if we have neutralized
            print('Minimize Combined Complex')
            copyfile('complex_neutralized.xyz','complex_to_min.xyz')
            self.minimize('complex_to_min.xyz','minimize.key') # minimize the combined structure
            os.rename('complex_to_min.xyz','complex_minimized.xyz')
            writeCheckpoint('Complex Minimized')


        if self.checkpoints < 8: # if we've minimized the combined structure
            print('Equilibrate Complex')
            copyfile('complex_minimized.xyz','complex_to_equil.xyz')
            self.equilibrate('complex_to_equil.xyz','equilibrate.key')
            killPeriodicInfo('complex_to_equil.arc')
            os.rename('complex_to_equil.xyz','complex_equil.xyz')
            os.rename('complex_to_equil.arc','complex_equil.arc')
            coarsenArcfile(self.MM_params, 'complex_equil.arc', 100)
            writeCheckpoint('Complex Equilibrated')


        if True: # always run sampling and binding analysis if everything else is done - self.checkpoints < 9: # if we've equilibrated the combined structure
            print('Run Dynamics')
            if os.path.exists('complex_sampled.arc'): # if we already sampled, resample
                pass
            else:
                copyfile('complex_equil.xyz','complex_to_sample.xyz')
            self.sampling('complex_to_sample.xyz','dynamics.key') # sample binding
            copyfile('complex_to_sample.xyz','complex_sampled.xyz')
            copyfile('complex_to_sample.arc','complex_sampled.arc')
            killPeriodicInfo('complex_sampled.arc')
            coarsenArcfile(self.MM_params, 'complex_sampled.arc', 1000)
            writeCheckpoint('Complex Sampled')

            # in the same step compute the binding
            print('Compute Binding')
            self.comProfile, self.bindingScore = [1,1] #evaluateBinding('complex_sampled.arc') # output a binding score
            writeCheckpoint('Got Binding Score')

            return self.comProfile, self.bindingScore

#=======================================================
#=======================================================
#=======================================================

    def getSecondaryStructure(self, sequence):
        '''
        get the secondary structure for a given sequence
        using seqfold here - identical features are available using nupack, though results are sometimes different
        :param sequence:
        :return: a dot-bracket string and list of paired bases (assuming single-strand DNA aptamer)
        '''
        temperature = 37.0
        dg(sequence, temp=temperature) # get energy of the structure
        #print(round(sum(s.e for s in structs), 2)) # predicted energy of the final structure

        structs = fold(sequence) # identify structural features
        desc = ["."] * len(sequence)
        pairList = []
        for s in structs:
            pairList.append(s.ij[0])
            pairList[-1]# list of bound pairs indexed from 1
            if len(s.ij) == 1:
                i, j = s.ij[0]
                desc[i] = "("
                desc[j] = ")"

        ssString = "".join(desc)
        pairList = np.asarray(pairList) + 1

        return ssString, pairList


    def foldSequence(self, sequence, pairList):
        '''
        generate a coarsely folded structure for the aptamer
        :param sequence: the ssDNA sequence to be folded
        :param pairList: list of binding base pairs
        :return:
        '''
        # write pair list as forces to the MMB command file
        comFile = 'commands.fold.dat'  # name of command file
        copyfile('commands.template.dat', comFile)  # make command file
        replaceText(comFile, 'SEQUENCE', sequence)


        baseString = '#baseInteraction A IND WatsonCrick A IND2 WatsonCrick Cis'
        lineNum = findLine(comFile, baseString)  # find the line number to start enumerating base pairs

        for i in range(len(pairList)):
            filledString = 'baseInteraction A {} WatsonCrick A {} WatsonCrick Cis'.format(pairList[i,0], pairList[i,1])
            addLine(comFile, filledString, lineNum + 1)

        # run fold
        os.system(self.params['mmb'] + ' -c ' + comFile + ' > outfiles/fold.out')
        os.rename('frame.pdb','sequence.pdb')
        os.system(self.params['pdbxyz path'] + ' sequence.pdb -k params_combined.key > outfiles/pdbconversion.out')


    def soak(self, structure, key):
        '''
        soak sample structure in solvent
        :param structure: structure to be solvated
        :param water: solvent
        :return: save wet .xyz file

        0.032 waters per cubic angstrom (slightly less than STP)
        nwaters = 0.032 * size**3
        '''
        # set params
        waterDensity = 0.032 # waters per cubic angstrom
        offset = self.params['box offset'] # box spacing from nearest atom in angstroms

        # get the size of the box to be made
        copyfile(structure, './to_be_wetted.arc')
        u = mda.Universe('to_be_wetted.arc')
        coords = u.atoms.positions
        xrange = np.ptp(coords[:,0])
        yrange = np.ptp(coords[:,1])
        zrange = np.ptp(coords[:,2])

        #''' more realistic box size (provided sequence is larger than analyte)
        maxsize = max([xrange,yrange,zrange])

        #xrange = max([xrange,maxsize/2])
        #yrange = max([yrange,maxsize/2])
        #zrange = max([zrange,maxsize/2])
        #'''

        xrange = np.ceil(xrange + 2*offset + 2*self.params['vdW radius']) # may also need an EWALD offset
        yrange = np.ceil(yrange + 2*offset + 2*self.params['vdW radius'])
        zrange = np.ceil(zrange + 2*offset + 2*self.params['vdW radius'])

        #xrange = max([xrange,yrange,zrange]) # make the box a cube
        #yrange = max([xrange,yrange,zrange])
        #zrange = max([xrange,yrange,zrange])

        numWaters = int(xrange*yrange*zrange*waterDensity)

        # edit waterbox.in with these parameters
        replaceText('waterbox.in','NWATERS','%d'%numWaters)
        replaceText('waterbox.in','XX','%.1f'%xrange)
        replaceText('waterbox.in','YY','%.1f'%yrange)
        replaceText('waterbox.in','ZZ','%.1f'%zrange)

        if os.path.exists('water.xyz_2'): # delete any prior waterboxes
            os.remove('water.xyz_2')

        os.system(self.params['xyzedit path'] +' < waterbox.in > outfiles/make_' + structure + '_water_box.out') # make box of waters

        # run via keyfile
        os.system(self.params['xyzedit path'] + ' ' + structure + ' -k ' + key + ' 12 0 > outfiles/' + structure + '_centered.out') # make sure structure is centered
        removetwo(structure)
        os.system(self.params['xyzedit path'] + ' ' + structure + ' -k ' + key + ' 20 ' + 'water.xyz_2 > outfiles/' + structure + 'wet.out') # place structure in box of waters
        removetwo(structure)
        box_data = copyLine('water.xyz_2',2) # copy periodic box data from the waterbox
        addLine(structure, box_data, 2) # re-add the periodic box data to the .xyz file

        # reset waterbox
        replaceText('waterbox.in','%d'%numWaters,'NWATERS')
        replaceText('waterbox.in','%.1f'%(xrange),'XX')
        replaceText('waterbox.in','%.1f'%(yrange),'YY')
        replaceText('waterbox.in','%.1f'%(zrange),'ZZ')

        return xrange, yrange, zrange


    def neutralize(self, structure, key):
        '''
        neutralize system charge with sodium ions
        number of positive charges added is in params

        also add appropriate NaCl concentration

        note - the old addIons2 works for amber - new addSodium and addChloride need amber codes appended if we want it to work
        :param structure:
        :param key:
        :return:
        '''

        # need all the atom numbers for the non-solvent parts - do this by finding the first water
        max_solute_index = findTXYZEndSoluteEnd(structure, self.MM_params)

        # convert from mmol to molecules / cubic angstrom
        # derivation for STP - will be slightly off for different temperatures
        densityScale = 0.032 / 55.5 # (water molecules per cubic angstrom divided by molar (pure water) - gives molecules / molar / A^3)

        saltDensity = densityScale * 1e-3 * self.params['NaCl concentration'] # scale concentration in mmol

        # compute number of molecules in the total box
        copyfile(structure, './to_be_neutralized.arc')
        u = mda.Universe('to_be_neutralized.arc')
        coords = u.atoms.positions
        xrange = np.ptp(coords[:,0])
        yrange = np.ptp(coords[:,1])
        zrange = np.ptp(coords[:,2])

        numSalts = int(saltDensity * xrange * yrange * zrange) # per-angstrom density * volume

        # compute number of sodiums to add
        nNa = int(numSalts + self.params['num charges']) # add the charges necessary to neutralize the analyte

        # compute number of chlorides to add

        nCl = int(numSalts)

        replaceText('addSodium.in','NUMIONS',str(nNa))
        replaceText('addSodium.in','SOLNUM',str(max_solute_index))
        replaceText('addChloride.in','NUMIONS',str(nCl))
        replaceText('addChloride.in','SOLNUM',str(max_solute_index))
        os.system(self.params['xyzedit path'] + ' ' + structure + ' -k ' + key + ' 21 < addChloride.in > outfiles/addIons.out')  # add 3 counterions - make sure you have a box already
        removetwo(structure)
        os.system(self.params['xyzedit path'] + ' ' + structure + ' -k ' + key + ' 21 < addSodium.in > outfiles/addIons.out')  # add 3 counterions - make sure you have a box already
        removetwo(structure)


    def addAnalyte(self, sequence, analyte, key):
        '''
        remove water, then add the analyte to the dna water box, then add water back
        :param analyte:
        :return:
        '''
        #replaceText('killWater.in','aaa.xyz',sequence) # edit killwater
        #os.system(self.params['xyzedit path'] + ' < killWater.in > outfiles/kill_waters_from_sequence.out') # remove waters
        #replaceText('killWater.in',sequence,'aaa.xyz') # restore file
        #removetwo(sequence)
        # move the analyte to the side of the sequence
        #os.system('xyzedit.exe ' + analyte + ' ' + self.MM_params + ' 11 20,20,20 0 > outfiles/shift_analyte.out') # command line won't work
        np.random.seed()
        copyfile(sequence, './pre_comb_seq.arc')
        u = mda.Universe('pre_comb_seq.arc')
        coords = u.atoms.positions
        seqCog = u.atoms.center_of_geometry()
        moveDist = np.ptp(coords)/2 + 10 # move it at least out of the sequence + 1 nm
        moveComp = np.sqrt(moveDist**2/3) # the components of the shift
        rands = np.random.randint(0,2,size=3)
        rands = (rands - 0.5) * 2
        moveComp = np.array((moveComp,moveComp,moveComp))
        moveComp = moveComp * rands # randomize the directional components of the move - the analyte can then approach from any of 8 directions
        xMove = moveComp[0] + seqCog[0]
        yMove = moveComp[1] + seqCog[1]
        zMove = moveComp[2] + seqCog[2]

        os.system(self.params['xyzedit path'] + ' < origin.in > outfiles/analyte_origin.out') # move analyte to the origin (not exactly center of mass) - the central ring nitrogen
        removetwo(analyte)
        replaceText('movesomething.in', 'aaa.xyz', analyte)
        replaceText('movesomething.in', 'XVEC', '%.1f'%xMove) # adjust input file
        replaceText('movesomething.in', 'YVEC', '%.1f'%yMove) # adjust input file
        replaceText('movesomething.in', 'ZVEC', '%.1f'%zMove) # adjust input file
        os.system(self.params['xyzedit path'] + ' < movesomething.in > outfiles/shift_analyte.out')
        #replaceText('movesomething.in', '%.1f'%xMove, 'XVEC')
        #replaceText('movesomething.in', '%.1f'%yMove, 'YVEC')
        #replaceText('movesomething.in', '%.1f'%zMove, 'ZVEC')
        #replaceText('movesomething.in', analyte, 'aaa.xyz')
        removetwo(analyte)
        os.system(self.params['xyzedit path'] + ' ' + sequence + ' -k ' + key + ' 18 ' + analyte + ' 0 > outfiles/add_analyte.out') # add analyte
        removetwo(sequence)


    def minimize(self, structure, keyfile):
        '''
        quickly minimize a given structure
        :param structure:
        :return:
        '''
        gradrms = self.params['minimization gradrms']
        os.system(self.params['minimize path'] + ' ' + structure + ' -k ' + keyfile + ' ' + str(gradrms)+ ' > outfiles/minimize_' + structure + '.out') # call tinker to minimize the structure
        removetwo(structure)


    def equilibrate(self, structure, keyfile):
        '''
        equilibrate given structure
        :param structure:
        :param keyfile:
        :return:
        '''
        simTime = self.params['equilibration time'] # time in ns
        dt = self.params['time step'] # timestep in fs
        timeSteps = int(simTime * 1e6 / dt)
        printAt = self.params['print step'] # printout spacing in ps #max(1.0,round(simTime * 1e3 / 143,1)) # print times in ps, always save a certain number of steps, one decimal point
        os.system(self.params['dynamic path'] + ' ' + structure + ' -k ' + keyfile + ' {} {} {} 2 310 >> outfiles/'.format(timeSteps,dt,printAt) + structure + '_equil_sequence.out')
        finished = 0
        reruns = 0
        max_reruns = 10
        finished, missingFrames = checkIfDynamicsFinished(structure, (simTime * 1e3) // printAt)
        while (finished == 0) and (reruns < max_reruns):
            finished, missingFrames = checkIfDynamicsFinished(structure, (simTime*1e3)//printAt)
            missingTime = (missingFrames // (1e3 * printAt))
            missingSteps = int(missingTime * 1e6 / dt)
            print('Rerunning Equilibration #%d' %reruns)
            os.system(self.params['dynamic path'] + ' ' + structure + ' -k ' + keyfile + ' {} {} {} 2 310 >> outfiles/'.format(missingSteps, dt, printAt) + structure + '_equil_sequence.out')
            reruns += 1

        #getFinalFrame(self.params['archive path'], structure)


    def sampling(self, structure, keyfile):
        '''
        run a trajectory two see if the molecules in the structure file bind
        :param structure:
        :param keyfile:
        :return:
        '''
        simTime = self.params['sampling time'] # time in ns
        dt = self.params['time step']# timestep in fs
        timeSteps = int(simTime * 1e6 / dt)
        printAt = self.params['print step'] # printout spacing in ps #max(1.0,round(simTime * 1e3 / 143,1)) # print times in ps, always save a certain number of steps, one decimal point
        os.system(self.params['dynamic path'] + ' ' + structure + ' -k ' + keyfile + ' {} {} {} 2 310 >> outfiles/'.format(timeSteps,dt,printAt) + structure + '_sampling.out')
        finished = 0
        reruns = 0
        max_reruns = 10
        finished, missingFrames = checkIfDynamicsFinished(structure, (simTime * 1e3) // printAt)
        while (finished == 0) and (reruns < max_reruns):
            finished, missingFrames = checkIfDynamicsFinished(structure, (simTime*1e3)//printAt)
            missingTime = (missingFrames // (1e3 * printAt))
            missingSteps = int(missingTime * 1e6 / dt)
            print('Rerunning Sampling #%d' %reruns)
            os.system(self.params['dynamic path'] + ' ' + structure + ' -k ' + keyfile + ' {} {} {} 2 310 >> outfiles/'.format(missingSteps, dt, printAt) + structure + '_sampling.out')
            reruns += 1

        #getFinalFrame(self.params['archive path'], structure)


'''
==============================================================
'''

if __name__ == '__main__':
    sequence =  'TATGCATGTGGGCGACGCAGTGCCCGTGGGATTTACTTGCAC' #'TATGCATGTGGGCGACGCAGTGCCCGTGGGATTTACTTGCAC' # stem truncated aptamer # #'TGCATGTGCGTGGGATTTACTTGCA' truncated central hairpin# 'TAGATCCGCATGAGGCTCGATCTGCATGTGGGCGACGCAGTGCCCGTGGGATTTACTTGCAC' # full aptamer #
    binder = binder(sequence, params)
    comProfile, bindingScore = binder.run() # retrieve binding score and center-of-mass time-series

    #os.chdir('C:/Users\mikem\Desktop/tinkerruns\clusterTests/fullRuns/run36')
    #comProfile, bindingScore = evaluateBinding('complex_sampling.arc') # normally done inside the run loop
    timeEq, potEq, kinEq = getTinkerEnergy('outfiles/complex_to_equil.xyz_equil_sequence.out') # get time series, potential and kinetic energies
    timeSa, potSa, kinSa = getTinkerEnergy('outfiles/complex_to_sample.xyz_sampling.out')

    outputs = {}
    outputs['time equil'] = timeEq
    outputs['pot equil'] = potEq
    outputs['kin equil'] = kinEq
    outputs['time sampling'] = timeSa
    outputs['pot sampling'] = potSa
    outputs['kin sampling'] = kinSa
    outputs['binding score'] = bindingScore
    outputs['com profile'] = comProfile
    outputs['params'] = params
    np.save('bindingOutputs', outputs) # unpack with load then outputs.item()
