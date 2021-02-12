import os
import json

def get_insertion(pdbpath):
    # hardcoding for now
    insertion_path = os.path.join(
            os.environ['HOME'],
            'cas', 'dels', 'rec2_helix',
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
