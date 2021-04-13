import os
import json
from klab.rosetta import input_files

def get_insertion(pdbpath, workspace):
    # hardcoding for now
    insertion_path = os.path.join(
            workspace.root_dir,
            'lucs', 'data', 'compatible'
            )
    model_number = os.path.basename(pdbpath).split('.')[0].split('_')[-1]
    insertion_json = os.path.join(
            insertion_path,
            'insertion_points_{}.json'.format(model_number)
            )
    with open(insertion_json, 'r') as f:
        insertion = json.load(f)
    # Return only first insertion for this script because there is only
    # one
    return insertion[0]

def relative_resfile(workspace, insertion):
    '''
    Convert a resfile to relative resfiles given a LUCS insertion.
    Note: Resfile/original PDB must have Rosetta numbering, i.e. start
    with the first residue having the number 1. Otherwise, the function
    will require the original pose and the LUCS pose.
    '''
    resfile_parser = input_files.Resfile(input_resfile=workspace.resfile_path)
    start = insertion['start']

    chain = 'A'
    designable = []
    for res in resfile_parser.design[chain]:
        if int(res) < insertion['start']:
            designable.append(int(res))
        else:
            designable.append(int(res) + (insertion['stop'] -
                insertion['start']))

    designable.extend([j for j in range(insertion['start'],
        insertion['stop'])])
    repackable = []
    if chain in resfile_parser.repack:
        for res in resfile_parser.repack[chain]:
            if int(res) < insertion['start']:
                repackable.append(int(res))
            else:
                repackable.append(int(res) + (insertion['stop'] -
                    insertion['start']))


    return designable, repackable


def relative_resfile_from_reference(workspace, insertion, pose,
        reference_pose):
    '''
    Convert a resfile to relative resfiles given a LUCS insertion.
    Note: Resfile/original PDB must have Rosetta numbering, i.e. start
    with the first residue having the number 1. Otherwise, the function
    will require the original pose and the LUCS pose.
    '''
    resfile_parser = input_files.Resfile(input_resfile=workspace.resfile_path)
    start = insertion['start']

    chain = 'A'
    designable = []
    for res in resfile_parser.design[chain]:
        if int(res) < insertion['start']:
            designable.append(int(res))
        else:
            designable.append(int(res) + (pose.size() -
                reference_pose.size()))

    designable.extend([j for j in range(insertion['start'],
        insertion['stop'])])
    repackable = []
    if chain in resfile_parser.repack:
        for res in resfile_parser.repack[chain]:
            if int(res) < insertion['start']:
                repackable.append(int(res))
            else:
                repackable.append(int(res) + (pose.size() -
                    reference_pose.size()))


    return designable, repackable
