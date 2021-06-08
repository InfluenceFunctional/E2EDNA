'''
utilities
'''
import os
import numpy as np
import MDAnalysis as mda
import glob
import argparse
import MDAnalysis.analysis.distances as distances
import scipy.ndimage as ndimage



def customPDBConversion(pdb,xyz,atomNumSource):
    '''
    convert pdb's with nonstandard atom-types to tiner xyz format
    requires also a pdb with the atoms in the same order as the pdb to be converted
    ########### DOESN'T CORRECTLY CONNECT PHOSPHATES ##################
    :param pdb:
    :param xyz:
    :return:
    '''
    # convert pdb to xyz with correct bonds
    atom_lines = "atom-lines.txt"
    test_csv = "test.csv"

    read_prm('amoebabio18.prm', atom_lines)
    try:
        lines, AMOEBA = fix_params(atom_lines, test_csv)
    except:
        pass

    system = load_pdb(pdb, AMOEBA)
    system = convert_names_AMOEBA(system, lines)
    write_xyz(system, 'intermediate')

    # now we have to kill the error message
    replaceText('intermediate.xyz','ATOM TYPE NOT FOUND','')

    # and write the correct atom types
    # first get them from the original pdb
    f = open(atomNumSource, 'r')
    text = f.read()
    f.close()
    text = text.split('\n')

    atomNums = []
    for i in range(len(text)):
        line = text[i]
        if 'HETATM' in line:
            line_split = line.split(' ')
            ind = 0
            for j in range(len(line_split)):
                if line_split[j] != '':
                    ind += 1
                    if ind == 4:
                        atomNums.append(line_split[j])

    # now write them to the xyz
    f = open('intermediate.xyz', 'r')
    text = f.read()
    f.close()
    text = text.split('\n')

    for i in range(1, len(text) - 1):
        line = text[i]
        line = line.replace(' 0 ', ' ' + str(atomNums[i - 1]) + ' ')
        text[i] = line

    text = "\n".join(text)

    # print the edited xyz
    f = open(xyz + '.xyz', 'w')
    f.write(text)
    f.close()


def get_input():
    '''
    get the command line in put for the run-num. defaulting to a new run (0)
    :return:
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--run_num', type=int, default = 0)
    cmd_line_input = parser.parse_args()
    run = cmd_line_input.run_num

    return run


def findLine(file, string):
    # return the line number of a given string, indexing from 1
    f = open(file, 'r')
    text = f.read()
    f.close()
    text = text.split('\n')
    try:
        lineInd = 1
        for line in text:
            if line == string:
                lineNum = lineInd # index from 1
            lineInd += 1

        return lineNum
    except:
        raise ValueError("String not found in file!")


def replaceText(file, old_string, new_string):
    # search and replace text in a file, then save the new version
    f = open(file, 'r')
    text = f.read()
    f.close()
    text = text.replace(old_string, new_string)
    f = open(file, 'w')
    f.write(text)
    f.close()


def removetwo(target):
    if os.path.exists(target + "_2"):
        os.remove(target) # delete original
        os.rename(target + "_2",target) # append updated
    else:
        print("No _2 structure was made!")
        raise ValueError


def copyLine(file,line_number):
    # copy a line of text from a file and return it
    f = open(file, 'r')
    text = f.read()
    f.close()
    return text.split('\n')[line_number - 1] # copy a line from the file, indexing from 1


def addLine(file, string, line_number):
    # add a new line of text to a file
    f = open(file, 'r')
    text = f.read()
    f.close()
    lines = text.split('\n')
    lines.insert(line_number - 1, string) # replace line ### with the string, indexing from 1
    text = "\n".join(lines)
    f = open(file, 'w')
    f.write(text)
    f.close()


def writeCheckpoint(text):
    '''
    write some output to the checkpoint file
    :return:
    '''
    f = open('checkpoint.txt','a')
    f.write('\n' + text)
    f.close()


def killPeriodicInfo(arcfile):
    '''
    delete periodic information from arcfiles so they display properly in vmd
    this can also be done using tinker's 'archive'
    :param arcfile:
    :return:
    '''
    f = open(arcfile,'r')
    text = f.read()
    f.close()
    lines = text.split('\n')
    periodicInfo = lines[1]
    lines = list(filter((periodicInfo).__ne__,lines))
    text = "\n".join(lines)
    f = open(arcfile, 'w')
    f.write(text)
    f.close()


def getPeriodicInfo(xyzfile):
    '''
    extract the box dimensions from a Tinker xyz file or arcfile
    :param xyzfile:
    :return: [xdim, ydim, zdim]
    '''
    f = open(xyzfile,'r')
    text = f.read()
    f.close()
    lines = text.split('\n')
    periodicInfo = lines[1]
    periodicInfo = periodicInfo.split('.')
    xdim = float(periodicInfo[0][-3:])
    ydim = float(periodicInfo[1][-3:])
    zdim = float(periodicInfo[2][-3:])

    return [xdim, ydim, zdim]


def getFinalFrame(archivePath, structure):
    '''
    use MDA to find out how many frames there are in the arc file
    then use tinker to extract the final frame
    :param structure:
    :return:
    '''
    arc = structure.split(".xyz")[0]+'.arc'
    u = mda.Universe(arc)
    frames = str(u.trajectory.n_frames)
    replaceText('grablastframe.in', 'XXX', arc)
    replaceText('grablastframe.in', 'FRAME', frames) # customize .in file
    os.system(archivePath + ' < grablastframe.in > outfiles/lastFrameGrab.out')
    replaceText('grablastframe.in', frames, 'FRAME') # reset .in file
    replaceText('grablastframe.in', arc, 'XXX')
    #os.rename(arc,'equil_'+arc) # rename trajectory
    flen = len(frames)
    if flen == 1:
        framestr = str('00' + frames)
    elif flen == 2:
        framestr = str('0' + frames)
    elif flen == 3:
        framestr = frames
    if os.path.exists(structure):
        if os.path.exists(structure.split('xyz')[0] + framestr):
            os.remove(structure)  # kill pre-dynamic structure
            os.rename(structure.split('xyz')[0] + framestr, structure) # rename last frame # will work up to 999 steps
    else:
        raise ValueError('Missing a structure!')


def getAFrame(archivePath, structure, frame):
    '''
    then use tinker to extract a given frame (with frame number < 1000)
    :param structure:
    :return:
    '''
    arc = structure.split(".xyz")[0] + '.arc'
    replaceText('grablastframe.in', 'XXX', arc)
    replaceText('grablastframe.in', 'FRAME', frame)  # customize .in file
    os.system(archivePath + ' < grablastframe.in > outfiles/framegrabber.out')
    replaceText('grablastframe.in', frame, 'FRAME')  # reset .in file
    replaceText('grablastframe.in', arc, 'XXX')
    # os.rename(arc,'equil_'+arc) # rename trajectory
    flen = len(frame)
    if flen == 1:
        framestr = str('00' + frame)
    elif flen == 2:
        framestr = str('0' + frame)
    elif flen == 3:
        framestr = frame
    if os.path.exists(structure):
        if os.path.exists(structure.split('xyz')[0] + framestr):
            os.rename(structure.split('xyz')[0] + framestr, 'grabbedFrame.xyz')  # rename # will work up to 999 steps
    else:
        raise ValueError('Missing a structure!')


def coarsenArcfile(mm_file,arcfile,frameNum):
    '''
    generate an arcfile with a certain maximum number of frames to reduce overall file size
    a nice upgrade would be to exclude water
    :param arcfile:
    :return:
    '''
    if 'amber' in mm_file:
        waterO = ' 2001'
        waterH = ' 2002'
    elif 'amoeba' in mm_file:
        waterO = ' 349'
        waterH = ' 350'

    u = mda.Universe(arcfile)
    noSolvent = u.select_atoms("not type" + waterO)
    noSolvent = noSolvent.select_atoms("not type" +waterH)

    traj_length = u.trajectory.__len__()
    if frameNum == 0: # if the input is 0, just take all the frames
        deltaStep = 1
    else:
        deltaStep = traj_length // frameNum
    if deltaStep >= 1:
        pass
    elif deltaStep < 1:
        deltaStep = 1

    with mda.Writer('coarse_' + arcfile +'.xyz',noSolvent.atoms.n_atoms) as W:
        for ts in u.trajectory:
            if (ts.frame % deltaStep) == 0:
                W.write(noSolvent)


def combineTrajectories(file1,file2):
    '''
    combine two coarsened xyz trajectories end-to-end and save
    :param file1:
    :param file2:
    :return:
    '''
    f = open(file1,'r')
    text1 = f.read()
    f.close()
    f = open(file2,'r')
    text2 = f.read()
    f.close()
    text = text1[:-1] + text2

    f = open('combined_trajectory.xyz', 'w')
    f.write(text)
    f.close()


def getTinkerEnergy(outfile):
    '''
    record time and energy from a tinker dynamic outfile
    :param outfile:
    :return:
    '''
    f = open(outfile,'r')
    text = f.read()
    f.close()

    text = text.split('\n')

    times = []
    potential = []
    kinetic = []

    # line-by-line
    for i in range(len(text)):
        line = text[i]
        if "Picosecond" in line:
            splitLine = line.split(' ')
            for j in range(len(splitLine)):
                if splitLine[j] == 'Picosecond':
                    ind = j-1 # the time and word "picosecond" are always separated by one space
                    break

            try:
                times.append(float(splitLine[ind]))
            except:
                times.append(np.nan)

        if "Potential" in line:
            splitLine = line.split(' ')
            for j in range(len(splitLine)):
                if splitLine[j] == 'Kcal/mole':
                    ind = j-1 # the energy and unit are always separated by one space
                    break

            try:
                potential.append(float(splitLine[ind]))
            except:
                potential.append(np.nan)

        if "Kinetic" in line:
            splitLine = line.split(' ')
            for j in range(len(splitLine)):
                if splitLine[j] == 'Kcal/mole':
                    ind = j - 1  # the energy and unit are always separated by one space
                    break

            try:
                kinetic.append(float(splitLine[ind]))
            except:
                kinetic.append(np.nan)

    return np.asarray(times), np.asarray(potential), np.asarray(kinetic)


def checkIfDynamicsFinished(structure, expectedFrames):
    '''
    check if a dynamics run is finished
    :param structure:
    :return:
    '''
    filename = structure.split('.xyz')[0]
    arcfile = filename + '.arc'
    if glob.glob(filename + '.err*'): #os.path.exists(filename + '.err'): # if an error file exists, then it has failed at least once
        u = mda.Universe(arcfile)
        frames = u.trajectory.n_frames
        missingFrames = expectedFrames - frames
        if missingFrames == 0:
            finished = 1
        elif missingFrames > 0:
            finished = 0
        else:
            raise ValueError('Error in calculation of missing printout steps')
    else:
        finished = 1
        missingFrames = 0

    return finished, missingFrames


def evaluateBinding(trajectory):
    '''
    evaluate the docking trajectory to determine whether the DNA binds to the analyte
    :param trajectory:
    :return:
    '''
    u = mda.Universe(trajectory)  # load up the trajectory for analysis
    # do
    # some
    # analysis

    # analyze trajectory
    comDist = []  # center of mass distance
    for ts in u.trajectory:
        # select backbone phosphates
        bBoneCom = u.select_atoms('type 343').center_of_mass()
        # select UTP
        analyteCom = u.select_atoms('type 4??').center_of_mass()
        comDist.append(np.linalg.norm(analyteCom - bBoneCom))

    comDist = np.asarray(comDist)

    bindingScore = np.sum(comDist < 20) / len(comDist)  # time within 2 nm
    return comDist, bindingScore
    ## maybe we could also have something like distance do nearest base?
    '''
    # GC binding for the truncated aptamer - identify hydrogen bonding group and track center of geomtetry distance
    # standard lengths in a watson-crick helix for each bond are in the 1.7-1.8 A range
    cogDist = []
    for ts in u.trajectory:
        CH = u.atoms.select_atoms('index 89')
        CN = u.atoms.select_atoms('index 84')
        CO = u.atoms.select_atoms('index 83')
        GO = u.atoms.select_atoms('index 726')
        GH1 = u.atoms.select_atoms('index 733')
        GH2 = u.atoms.select_atoms('index 735')
        
        C_group = CH + CN + CO
        G_group = GO + GH1 + GH2
        
        cCog = C_group.center_of_geometry()
        gCog = G_group.center_of_geometry()
        cogDist.append(np.linalg.norm(gCog - cCog))

    '''


def trajAnalysis(trajectory,analyte):
    '''
    evaluate the docking trajectory to determine whether the DNA binds to the analyte
    also whether the DNA rearranges at critical base pairing sites
    :param trajectory:
    :return:
    '''
    u = mda.Universe(trajectory)  # load up the trajectory for analysis
    # do
    # some
    # analysis

    # customized binding for the truncated aptamer - identify hydrogen bonding group and track center of geomtetry distance, at several points
    # works for xyz or arcfiles
    # also track UTP distance from binding sites
    # standard lengths in a watson-crick helix for each bond are in the 1.7-1.8 A range
    dist1 = [] # mid-loop pair dist 28-33
    dist2 = [] # right loop pair dist 10-25
    dist3 = [] # left loop pair dist 9-35
    dist4 = [] # utp-top binding site ~13
    dist5 = [] # utp-left binding site ~36
    for ts in u.trajectory:
        d1_1 = u.atoms.select_atoms('index 889').center_of_geometry() # base 28 - guanine ring amino proton
        d1_2 = u.atoms.select_atoms('index 1049').center_of_geometry() # base 33 - thymine ring amino proton

        d2_1 = u.atoms.select_atoms('index 317').center_of_geometry() # base 10 - guanine ring amino proton
        d2_2 = u.atoms.select_atoms('index 785').center_of_geometry() # base 25 - cytosine ring nitrogen

        d3_1 = u.atoms.select_atoms('index 282').center_of_geometry() # base 9 - thymine ring amino proton
        d3_2 = u.atoms.select_atoms('index 1110').center_of_geometry() # base 35 - adenine ring amino



        dist1.append(np.linalg.norm(d1_1-d1_2))
        dist2.append(np.linalg.norm(d2_1-d2_2))
        dist3.append(np.linalg.norm(d3_1-d3_2))

        if analyte == True: # for false - just analyze the DNA
            d4_1 = u.atoms.select_atoms('index 1365').center_of_geometry()  # utp ribose ring oxygen
            d4_2 = u.atoms.select_atoms('index 407').center_of_geometry()  # base 13 - cytosine ring nitrogen

            d5_1 = u.atoms.select_atoms('index 1365').center_of_geometry()  # utp ribose ring oxygen
            d5_2 = u.atoms.select_atoms('index 1139').center_of_geometry()  # base 36 - cytosine ring nitrogen
            dist4.append(np.linalg.norm(d4_1-d4_2))
            dist5.append(np.linalg.norm(d5_1-d5_2))


    '''
    # nice figure
    plt.figure()
    plt.clf()
    plt.plot(dist1,label='inter-ring')
    plt.plot(dist2,label='ring top')
    plt.plot(dist3,label='ring left')
    plt.plot(dist4,label='UTP-site1')
    plt.plot(dist5,label='UTP-site2')
    plt.legend()
    '''

    return np.asarray(dist1), np.asarray(dist2), np.asarray(dist3), np.asarray(dist4), np.asarray(dist5)


def rolling_mean(input, run):
    output = np.zeros(len(input))
    for i in range(len(output)):
        if i < run:
            output[i] = np.average(input[0:i])
        else:
            output[i] = np.average(input[i - run:i])

    return output


def findTXYZEndSoluteEnd(structure,MM_params):
    '''
    find the final line of a tinker xyz file, to append ions
    :return: the line number
    '''
    # find the last solute atom in the xyz file (find the first water)
    # return the line number of a given string, indexing from 1

    f = open(structure, 'r')
    text = f.read()
    f.close()
    text = text.split('\n')
    if MM_params == '../params/amoebabio18.prm':
        waterO = ' 349 '
        waterH = ' 350 '
    elif MM_params == '../params/amber99.prm':
        waterO = ' 2001 '
        waterH = ' 2002 '

    # find the first line of the xyz file that is a water
    try:
        for i in range(len(text)):
            line = text[i]
            # if (' 349 ' in line) and (' 350 ' in text[i+1]) and (' 350 ' in text[i + 2]) and (' 349' in text[i + 3]): # look for waters, atom type 349 and 350 in amoebabio18
            if (waterO in line) and (waterH in text[i + 1]) and (waterH in text[i + 2]) and (waterO in text[i + 3]):  # look for waters, atom type 349 and 350 in amoebabio18
                max_solute_line = text[i - 1]  # index from 1
                break

    except:
        raise ValueError("No waters in structure!")

    solLine = max_solute_line.split(' ')
    for i in range(len(solLine)):  # extract the first integer
        try:
            max_solute_index = int(solLine[i])
            break
        except:
            aa = 1

    return max_solute_index


def saveOutputs(params,reactionCoordinates):
    '''
    save simulation outputs
    :return:
    '''
    outputs = {}
    outputs['reaction coordinates'] = reactionCoordinates
    outputs['params'] = params
    np.save('bindingOutputs', outputs) # do these first in case the analysis fails

    timeEq, potEq, kinEq = getTinkerEnergy('outfiles/complex_to_equil.xyz_equil.out') # get time series, potential and kinetic energies
    timeSa, potSa, kinSa = getTinkerEnergy('outfiles/complex_to_sample.xyz_sampling.out')

    outputs['time equil'] = timeEq
    outputs['pot equil'] = potEq
    outputs['kin equil'] = kinEq
    outputs['time sampling'] = timeSa
    outputs['pot sampling'] = potSa
    outputs['kin sampling'] = kinSa
    outputs['reaction coordinates'] = reactionCoordinates
    outputs['params'] = params
    np.save('bindingOutputs', outputs) # unpack with np.load('bindingOutputs.npy',allow_pickle=True) then outputs=outputs.item()


def minimaAnalysis(rcTrajectories,freeEnergies,freeEnergyAxes):
    '''
    analyze free energy profiles and isolate the separate and joint minima
    a more sophisticated approach would be to model the joint multidimensional free energy distribution
    rather than assuming separate linear contributions
    :param rcTrajectories:
    :param freeEnergies:
    :param freeEnergyAxes:
    :return:
    '''
    # identify free energy minima
    freeEnergyMinima = np.zeros(len(freeEnergies))
    for i in range(len(freeEnergies)):
        freeEnergyMinima[i] = freeEnergyAxes[i][np.argmin(freeEnergies[i])]

    # identify RC trajectory frames which come within xx% of the minima
    cutoff = 0.1  # criteria for being 'inside' the minima
    minimaFrames = np.zeros((len(rcTrajectories), len(freeEnergies)))
    for i in range(len(freeEnergies)):
        minimaFrames[:, i] = (np.abs(rcTrajectories[:, i] - freeEnergyMinima[i]) / freeEnergyMinima[0]) < cutoff

    jointMinima = np.nonzero((np.sum(minimaFrames, axis=1) == len(freeEnergies)).astype(int))  # find the frames which simultaneously satisfy all the criteria
    # with thorough sampling, would could more rigorously do this via multidimensional free energy analysis using histogramdd
    return jointMinima


def extractTrajectory(trajectory, reactionCoordinates, boxSize, equilibrationTime):
    '''
    ONLY FOR THE 'REPRESENTATIVE STRUCTURE' BINDING RUNS
    :param trajectory:
    :return:
    '''
    u = mda.Universe(trajectory)  # load up the trajectory for analysis
    # do
    # some
    # analysis

    rcTrajectories = np.zeros((u.trajectory.n_frames,len(reactionCoordinates))) # initialize RC trajectories

    tt = 0
    for ts in u.trajectory: # for each frame
        if tt > equilibrationTime:
            for i in range(len(reactionCoordinates)): # for each reaction coordinate
                d1 = u.atoms.select_atoms('index %d'%reactionCoordinates[i][0]) # atom 1
                d2 = u.atoms.select_atoms('index %d'%reactionCoordinates[i][1]) # atom 2
                rcTrajectories[tt, i] = distances.dist(d1,d2,box=(boxSize[0],boxSize[1],boxSize[2],90,90,90))[-1] # 'box' information accounts for periodicity - assuming cubic periodicity

        tt+= 1

    return rcTrajectories


def getFreeEnergy(trajectory,sigma):
    '''
    use histogramming and gaussian smoothing to generate free energy profile from a given trajectory
    :param trajectory:
    :return:
    '''
    pop, edges = np.histogram(trajectory, bins=100)
    bin_mids = edges[0:-1] + np.diff(edges)
    freeEnergy = -np.log(ndimage.gaussian_filter1d(pop/np.sum(pop), sigma))

    return bin_mids,freeEnergy


