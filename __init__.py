#### by Michael Szegedy, 2017 Jun
#### for personal use at Khare Lab at Rutgers University

'''"mszegedy's PyRosetta extension" is just some functions and classes I've
needed to use in multiple scripts. Do not take it too seriously.
'''

import os
import fcntl, itertools, json, re, time
import pyrosetta as pr

def add_constraints_from_file(pose, filename):
    '''Reads a new-style (non-enzdes) .cst file and adds the constraints in it
    to a pose.
    '''
    # For the curious: passing pose doesn't actually accomplish anything. It's
    # required because they were gonna have this function overwrite the
    # constraint set of the pose at the end, but as of the time of this
    # writing, that line is commented out.
    cst = pr.rosetta.core.scoring.constraints.ConstraintIO \
            .read_constraints_new(filename, pose.constraint_set(), pose)
    pose.constraint_set(cst)
    return cst

def atom_dist(pose, atom1, atom2):
    '''Takes in a pose and two AtomIDs. Returns the distance of the two atoms from
    each other.
    '''
    conformation = pose.conformation()
    return (conformation.residue(atom1.rsd()).xyz(atom1.atomno()) - \
            conformation.residue(atom2.rsd()).xyz(atom2.atomno())).norm()

def dump_numbered_pdb(pose, prefix, directory=None):
    '''Checks for files in the current directory with prefix 'prefix', and safely
    dumps pose to the next available numbered filename with that prefix.
    '''
    ## get directory
    if directory is None:
        directory = os.getcwd()
    else:
        directory = os.path.join(os.getcwd(), directory)
    ## get number of this particular PDB
    lockfile_path = os.path.join(directory, '.mpre-dump-lockfile-'+prefix+'~')
    while True:
        try:
            with open(lockfile_path, 'x') as lockfile:
                # We construct a list where each index, apart from 0, contains True
                # or False, corresponding to whether the file with the same number
                # as the index exists or not.
                #
                # NOTE FROM THE FUTURE: This was originally done because the list
                #   was cached in a file for reuse; then the file was found to
                #   become outdated too easily, so I turned it into a dummy
                #   lockfile and instead recalculated the table every time. For
                #   sure can this be rewritten to be more concise, but who really
                #   cares? Maybe the list'll become more useful in the future.
                numtable = [None]
                filenames_in_dir = os.listdir(directory)
                nummatcher = re.compile(re.escape(prefix)+r'(\d+)\.pdb')
                matchnums = []
                for filename in filenames_in_dir:
                    match = nummatcher.fullmatch(filename)
                    if match:
                        try:
                            matchnum = int(match.group(1))
                            if matchnum >= 1:
                                matchnums.append(matchnum)
                        except ValueError:
                            pass
                matchnums.sort()
                for matchnum in matchnums:
                    for x in range(matchnum-len(numtable)):
                        numtable.append(False)
                    numtable.append(True)
                ournum = None
                try:
                    ournum = numtable.index(False)
                    # numtable[ournum] = True # now redundant
                except ValueError:
                    ournum = len(numtable)
                    # numtable.append(True) # now redundant
                ## save PDB
                output_path = os.path.join(directory,
                                           prefix + str(ournum) + ".pdb")
                retval = pose.dump_pdb(output_path)
            os.remove(lockfile_path)
            return retval
        except FileExistsError:
            time.sleep(0.05)
        except:
            os.remove(lockfile_path)
            raise

def pose_from_file_with_params(filename, params):
    '''This will create a Pose() object from a PDB file that needs params. The
    params are to be given as a list of string literals of filenames.
    '''
    pose = pr.Pose()
    pr.generate_nonstandard_residue_set(pose, params)
    pr.pose_from_file(pose, filename)
    return pose

def pose_from_pdbstring_with_params(pdbstring, params):
    '''This will create a Pose() object from a PDB contents string that needs
    params. The params are to be given as a list of string literals of
    filenames.
    '''
    pose = pr.Pose()
    pr.generate_nonstandard_residue_set(pose, params)
    pr.rosetta.core.import_pose.pose_from_pdbstring(pose, pdbstring)
    return pose

def pose_from_file(filename, params=None):
    '''This will create a Pose() object from a PDB file that needs params. The
    params are to be given as a list of string literals of filenames.
    '''
    if params is not None:
        return pose_from_file_with_params(filename, params)
    else:
        return pr.pose_from_file(pr.Pose(), filename)

def pose_from_pdbstring(pdbstring, params=None):
    '''This will create a Pose() object from a PDB contents string that needs
    params. The params are to be given as a list of string literals of
    filenames.
    '''
    if params is not None:
        return pose_from_pdbstring_with_params(pdbstring, params)
    else:
        return pr.rosetta.core.import_pose.pose_from_pdbstring(pr.Pose(),
                                                               pdbstring)

def res_neighbors_p(pose, resnum1, resnum2, coarsep=False, bound=None):
    '''This compares the conformations of two residues in a pose, and returns a
    boolean of whether or not they are physically close to each other. There
    are two possible standards that the algorithm employs:
        coarsep=False: any two non-hydrogen atoms in the two residues'
          sections of protein are 4.2 Ao or less apart
        coarsep=True: CA atoms of residues are 8 Ao or less apart
    '''
    # get ready for some lisp bullshit
    bound = bound or \
            (coarsep and 8.) or \
            4.2
    conformation = pose.conformation()
    if coarsep:
        # dunno what kind of residue doesn't have CA as the second atom, but we
        # are now prepared even then
        return atom_dist(pose,
                         pr.AtomID(conformation.residue(resnum1).atom_index('CA'),
                                   resnum1),
                         pr.AtomID(conformation.residue(resnum2).atom_index('CA'),
                                   resnum2)) < bound
    else:
        comparisons = itertools.product(
                        (pr.AtomID(n, resnum1) \
                         for n \
                         in range(1,conformation.residue(resnum1).nheavyatoms()+1)),
                        (pr.AtomID(n, resnum2) \
                         for n \
                         in range(1,conformation.residue(resnum2).nheavyatoms()+1)))
        retval = False
        for comparison in comparisons:
            retval = retval or \
                     (atom_dist(pose, comparison[0], comparison[1]) < bound)
        return retval
