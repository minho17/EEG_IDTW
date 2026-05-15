import numpy as np
from scipy.io import loadmat
from scipy.stats import ttest_ind,levene

def main():

    path_data = "./biomarker.mat" # File path
    mat_file = loadmat(path_data)

    biomarker = mat_file['EEG_biomarker'][0,:] # 'DiCohBI' in the paper
    ind0 = mat_file['index_HC'][0,:] # index for HC group after pre-processing
    ind1 = mat_file['index_MCI'][0,:] # index for MCI group after pre-processing

    [p_value,t,_] = t_test_ex(biomarker,ind0,ind1) 
    print(p_value, t)



def t_test_ex(data,ind0,ind1):

    x = data[ind0]
    y = data[ind1]

    ind = np.argwhere(x > -999)[:,0]
    x2 = x[ind] 
    ind = np.argwhere(y > -999)[:,0]
    y2 = y[ind]

    _, p2 = levene(x2,y2)

    if p2 > 0.05:
        t,p = ttest_ind(x2,y2)
    else:
        t,p = ttest_ind(x2,y2,equal_var= False)

    return p,t,p2


if __name__ == '__main__':
    main()
    