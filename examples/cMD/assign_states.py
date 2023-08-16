import numpy
from tqdm.auto import tqdm, trange

def assign_states(input_array):

    state_list = []
    
    for val in tqdm(input_array):
        # Assign source and target states
        if val[0] >= 25 and val[0] <= 90 and val[1] >= -55 and val[1] <=0: # C7ax
            state_list.append(0)
        elif val[0] >= -170 and val[0] <= -55 and val[1] >= 40 and val[1] <= 100: # C7eq
            state_list.append(1)
        else:
            state_list.append(2)

    return state_list
