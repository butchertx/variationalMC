import json
import os


K = [-2.0, -1.6, -1.2, -1.0, -0.8, -0.4, 0.0, 0.4, 0.8, 1.0, 1.2, 1.6, 2.0]
for kidx, k in enumerate(K):
    with open('neel.json') as f:
        input_file = json.load(f)
    
    term_list = input_file['model']['bilinear']
    new_term_list = []
    for term in term_list:
        if term['type'] == "su2 biquadratic":
            term['coupling'] = k

        new_term_list.append(term)

    input_file['model']['bilinear'] = new_term_list
    
    rundir = "K" + str(kidx)
    try:
        os.mkdir(rundir)
    except OSError as error:
        print(f"Directory {rundir} already exists")

    with open(os.path.join(rundir, 'neel.json'), 'w') as outfile:
        json.dump(input_file, outfile)
