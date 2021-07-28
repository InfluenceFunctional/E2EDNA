import re
from utils import *
from seqfold import dg, fold
from shutil import copyfile


class e2edna():
    def __init__(self,sequence,params):
        self.params = params
        self.sequence = sequence

        if self.params['force field'] == 'AMOEBA':
            self.MM_params = '../params/amoebabio18.prm'

        self.setup()
        self.getCheckpoint(self.params['simulation type'])


    def setup(self):
        '''
        setup working directory
        copy in relevant xyz and keyfiles
        set simulation parameters in keyfiles
        move to relevant directory
        '''
        if not os.path.exists(self.params['workdir'] + '/params/amoebabio18.prm'): # all runs share a params directory - make sure it exists
            copyfile('lib/params/amoebabio18.prm',self.params['workdir'] + '/params/')

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
            copyfile(self.params['grablastframe'], self.workDir + '/grablastframe.in')
            copyfile(self.params['movesomething'], self.workDir + '/movesomething.in')
            copyfile(self.params['killWater'], self.workDir + '/killWater.in')
            copyfile(self.params['addIons'], self.workDir + '/addIons.in') # not necessary anymore
            copyfile(self.params['addSodium'], self.workDir + '/addSodium.in')
            copyfile(self.params['addChloride'], self.workDir + '/addChloride.in')
            copyfile(self.params['origin'], self.workDir + '/origin.in')
            copyfile(self.params['mmb template'], self.workDir + '/commands.template.dat')


            if self.params['outisde secondary structure']:
                copyfile('fullPairList.npy', self.workDir + '/fullPairList.npy')


        else:
            self.workDir = self.params['workdir'] + '/' + 'run%d' %self.params['run num']


        # move to working dir
        os.chdir(self.workDir)

        # save params to outputs
        outputs = {}
        outputs['params'] = self.params
        np.save('e2ednaOutputs', outputs)  # before running - save the inputs


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


    def getCheckpoint(self,mode):
        '''
        use the checkpoint file to identify if any work has been done in this directory, and if so, where to pick up
        :return:
        '''
        checkpointFile = mode + '_checkpoint.txt'
        try:
            f = open(checkpointFile,'r')#if it does, see how long it is
            text = f.read()
            f.close()
            self.checkpoints = len([m.start() for m in re.finditer('\n', text)])
            print("Resuming Run #%d"%int(self.params['run num']))

        except:
            f = open(checkpointFile,'w')
            f.write('Started!\n')
            f.close()
            self.checkpoints = 1


    def runFreeAptamer(self):
        '''
        run e2edna
        consult checkpoints and operating mode to determine correct actions
        :return:
        '''
        if self.checkpoints < 2: # if we have initialized
            # get secondary structure
            print("Get Secondary Structure")
            if os.path.exists("pre_fold.pdb"): # if we have a pre-folded structure, do nothing
                pass
            elif os.path.exists('fullPairList.npy'): # if we have a manually-input pair list, use that
                self.pairList = np.asarray(np.load('fullPairList.npy',allow_pickle=True).item()['list'])
            else:
                self.ssString, self.pairList = self.getSecondaryStructure(self.sequence) # otherwise predict secondary structure
            writeCheckpoint("Got Secondary Structure", self.params['simulation type'])


        if self.checkpoints < 3: # if we have the secondary structure
            # fold it!
            print("Folding Sequence")
            if os.path.exists("pre_fold.pdb"):
                copyfile('pre_fold.pdb','sequence.pdb')
            else:
                self.foldSequence(self.sequence, self.pairList)
                pass
            writeCheckpoint("Folded Sequence", self.params['simulation type'])


        if self.checkpoints < 4: # if we added the analyte
            print('Soak Sequence')
            copyfile('sequence.xyz','sequence_pre_soak.xyz')
            self.soak('sequence.xyz','params_combined.key')
            os.rename('sequence.xyz', 'sequence_soaked.xyz')
            writeCheckpoint('Sequence Soaked', self.params['simulation type'])


        if self.checkpoints < 5: # if we have soaked it
            # neutralize the sequence and add salt
            print('Neutralize Sequence')
            copyfile('sequence_soaked.xyz','sequence_to_neutralize.xyz')
            self.neutralize('sequence_soaked.xyz','params_combined.key')
            os.rename('sequence_soaked.xyz','sequence_neutralized.xyz')
            writeCheckpoint('Sequence Neutralized', self.params['simulation type'])


        if self.checkpoints < 6: # if we have neutralized
            print('Minimized Soaked Sequence')
            copyfile('sequence_neutralized.xyz','sequence_to_min.xyz')
            self.minimize('sequence_neutralized.xyz','minimize.key') # minimize the combined structure
            os.rename('sequence_neutralized.xyz','sequence_minimized.xyz')
            writeCheckpoint('Sequence Minimized', self.params['simulation type'])


        if self.checkpoints < 7: # if we've minimized the combined structure
            print('Equilibrate Sequence')
            copyfile('sequence_minimized.xyz','sequence_to_equil.xyz')
            self.equilibrate('sequence_to_equil.xyz','equilibrate.key')
            os.rename('sequence_to_equil.xyz','sequence_equil.xyz')
            os.rename('sequence_to_equil.arc','sequence_equil.arc')
            coarsenArcfile(self.MM_params, 'sequence_equil.arc', 0)
            writeCheckpoint('Sequence Equilibrated', self.params['simulation type'])


        if self.checkpoints <8: # always run sampling and binding analysis if everything else is done
            print('Run Dynamics')
            if os.path.exists('sequence_sampled.arc'): # if we already sampled, extend the trajectory
                pass
            else:
                getFinalFrame(self.params['archive path'], 'sequence_equil.xyz') # only works if there are less than 1k frames in equilibration trajectory
                copyfile('sequence_equil.xyz','sequence_to_sample.xyz')
            self.sampling('sequence_to_sample.xyz','dynamics.key') # sample binding
            copyfile('sequence_to_sample.xyz','sequence_sampled.xyz')
            copyfile('sequence_to_sample.arc','sequence_sampled.arc')
            killPeriodicInfo('sequence_sampled.arc')
            coarsenArcfile(self.MM_params, 'sequence_sampled.arc', 0)
            writeCheckpoint('Sequence Sampled', self.params['simulation type'])

            # in the same step compute the binding
            print('Do Analysis')
            analysisDict = self.trajectoryAnalysis('sequence','coarse_sequence_sampled.arc.xyz',self.params['reaction coordinates'], 0, 2) # output a binding score
            getAFrame(self.params['archive path'], 'sequence_sampled.arc', analysisDict['representative structure frames'][-1][-1]) # save the last representative frame
            os.rename('grabbedFrame.xyz','representativeSequence.xyz')
            killWater(self.params['xyzedit path'], 'representativeSequence.xyz')

            return analysisDict


    def runBinding(self):
        '''
        run e2edna
        consult checkpoints and operating mode to determine correct actions
        :return:
        '''

        if self.checkpoints < 2: # if we've folded the sequence
            # add the analyte!
            print('Add Analyte')
            self.addAnalyte('representativeSequence.xyz', 'analyte.xyz','params_combined.key')
            os.rename('representativeSequence.xyz','complex.xyz')
            writeCheckpoint('Added Analyte', self.params['simulation type'])


        if self.checkpoints < 3: # if we added the analyte
            print('Soak Complex')
            copyfile('complex.xyz','complex_pre_soak.xyz')
            self.soak('complex.xyz','params_combined.key')
            os.rename('complex.xyz', 'complex_soaked.xyz')
            writeCheckpoint('Complex Soaked', self.params['simulation type'])


        if self.checkpoints < 4: # if we have soaked it
            # neutralize the analyte
            print('Neutralize Complex')
            copyfile('complex_soaked.xyz','complex_to_neutralize.xyz')
            self.neutralize('complex_soaked.xyz','params_combined.key')
            os.rename('complex_soaked.xyz','complex_neutralized.xyz')
            writeCheckpoint('Complex Neutralized', self.params['simulation type'])


        if self.checkpoints < 5: # if we have neutralized
            print('Minimize Combined Complex')
            copyfile('complex_neutralized.xyz','complex_to_min.xyz')
            self.minimize('complex_neutralized.xyz','minimize.key') # minimize the combined structure
            os.rename('complex_neutralized.xyz','complex_minimized.xyz')
            writeCheckpoint('Complex Minimized', self.params['simulation type'])


        if self.checkpoints < 7: # if we've minimized the combined structure
            print('Equilibrate complex')
            copyfile('complex_minimized.xyz','complex_to_equil.xyz')
            self.equilibrate('complex_to_equil.xyz','equilibrate.key')
            os.rename('complex_to_equil.xyz','complex_equil.xyz')
            os.rename('complex_to_equil.arc','complex_equil.arc')
            coarsenArcfile(self.MM_params, 'complex_equil.arc', 0)
            writeCheckpoint('complex Equilibrated', self.params['simulation type'])


        if self.checkpoints <8: # always run sampling and binding analysis if everything else is done
            print('Run Dynamics')
            if os.path.exists('complex_sampled.arc'): # if we already sampled, resample
                pass
            else:
                getFinalFrame(self.params['archive path'], 'complex_equil.xyz') # only works if there are less than 1k frames in equilibration trajectory
                copyfile('complex_equil.xyz','complex_to_sample.xyz')

            self.sampling('complex_to_sample.xyz','dynamics.key') # sample binding
            copyfile('complex_to_sample.xyz','complex_sampled.xyz')
            copyfile('complex_to_sample.arc','complex_sampled.arc')
            killPeriodicInfo('complex_sampled.arc')
            coarsenArcfile(self.MM_params, 'complex_sampled.arc', 0)
            writeCheckpoint('complex Sampled', self.params['simulation type'])

            # in the same step compute the binding
            print('Compute Binding')
            analysisDict = self.trajectoryAnalysis('complex','coarse_complex_sampled.arc.xyz',self.params['reaction coordinates'], 0, 2) # output a binding score
            writeCheckpoint('Got Binding Score', self.params['simulation type'])

            return analysisDict


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

        replaceText('minimize.key', 'XX', '%.1f' % xrange)
        replaceText('minimize.key', 'YY', '%.1f' % yrange)
        replaceText('minimize.key', 'ZZ', '%.1f' % zrange)
        replaceText('equilibrate.key', 'XX', '%.1f' % xrange)
        replaceText('equilibrate.key', 'YY', '%.1f' % yrange)
        replaceText('equilibrate.key', 'ZZ', '%.1f' % zrange)
        replaceText('dynamics.key', 'XX', '%.1f' % xrange)
        replaceText('dynamics.key', 'YY', '%.1f' % yrange)
        replaceText('dynamics.key', 'ZZ', '%.1f' % zrange)


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
        if self.params['simulation type'] == 'free':
            self.params['num charges'] = 0
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
        # move the analyte to the side of the sequence
        if self.params['analyte position'] == 'random':
            np.random.seed()
            copyfile(sequence, './pre_comb_sequence.arc')
            u = mda.Universe('pre_comb_sequence.arc')
            coords = u.atoms.positions
            seqCog = u.atoms.center_of_geometry()
            offset = 20 # angstroms
            moveDist = np.ptp(coords)/2 + 20 # move it at least out of the sequence in longest dimension + 'offset' angstroms
            moveComp = np.sqrt(moveDist**2/3) # the components of the shift
            rands = np.random.randint(0,2,size=3)
            rands = (rands - 0.5) * 2
            moveComp = np.array((moveComp,moveComp,moveComp))
            moveComp = moveComp * rands # randomize the directional components of the move - the analyte can then approach from any of 8 directions
            xMove = moveComp[0] + seqCog[0]
            yMove = moveComp[1] + seqCog[1]
            zMove = moveComp[2] + seqCog[2]
        elif self.params['analyte position'] == 'manual':
            xMove = self.params['analyte coordinates'][0]
            yMove = self.params['analyte coordinates'][1]
            zMove = self.params['analyte coordinates'][2]


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
        os.system(self.params['dynamic path'] + ' ' + structure + ' -k ' + keyfile + ' {} {} {} 2 310 >> outfiles/'.format(timeSteps,dt,printAt) + structure + '_equil.out')
        finished = 0
        reruns = 0
        max_reruns = 10
        finished, missingFrames = checkIfDynamicsFinished(structure, (simTime * 1e3) // printAt)
        while (finished == 0) and (reruns < max_reruns):
            finished, missingFrames = checkIfDynamicsFinished(structure, (simTime*1e3)//printAt)
            missingTime = (missingFrames * printAt / 1e3)
            missingSteps = int(missingTime * 1e6 / dt)
            print('Rerunning Equilibration #%d' %reruns)
            os.system(self.params['dynamic path'] + ' ' + structure + ' -k ' + keyfile + ' {} {} {} 2 310 >> outfiles/'.format(missingSteps, dt, printAt) + structure + '_equil.out')
            reruns += 1


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
            missingTime = (missingFrames * printAt / 1e3)
            missingSteps = int(missingTime * 1e6 / dt)
            print('Rerunning Sampling #%d' %reruns)
            os.system(self.params['dynamic path'] + ' ' + structure + ' -k ' + keyfile + ' {} {} {} 2 310 >> outfiles/'.format(missingSteps, dt, printAt) + structure + '_sampling.out')
            reruns += 1


    def trajectoryAnalysis(self,type,trajectory,reactionCoordinates,equilibrationTime,sigma):
        '''
        follow a set of reaction coordinates, bin the trajectory and generate a free energy profile
        :param reactionCoordinates: list of pairs of integers representing the atom numbers of desired distances
        :param equilibrationTime: number of frames from the beginning of the trajectory to omit as 'equilibration'
        :return: rcTrajectories, free energy curves and bin centers, 'representative' structures
        '''
        boxSize = getPeriodicInfo(type + '_to_sample.xyz')
        # get the trajectories of the various reaction coordinates (atom-atom distances)
        rcTrajectories = extractTrajectory(trajectory,reactionCoordinates,boxSize,equilibrationTime)

        # convert trajectories into free energy curves
        freeEnergies = []
        freeEnergyAxes = []
        for i in range(len(reactionCoordinates)):
            axis, freeEnergy = getFreeEnergy(rcTrajectories[:,i],sigma)
            freeEnergies.append(freeEnergy) # in units of kT
            freeEnergyAxes.append(axis)

        representativeFrames = minimaAnalysis(rcTrajectories,freeEnergies,freeEnergyAxes)

        analysisDict = {}
        analysisDict['reaction coordinate trajectories'] = rcTrajectories
        analysisDict['free energies'] = freeEnergies
        analysisDict['free energy axes'] = freeEnergyAxes
        analysisDict['representative structure frames'] = representativeFrames

        return analysisDict