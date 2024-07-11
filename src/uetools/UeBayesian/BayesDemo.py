import numpy as np
from uetools.UeCase import Case
import os

# Synthetic data to calculate mean-square error between UEDGE and observation
psi_syn = np.array([0.880,0.891,0.901,0.910,0.919,0.930,0.939,0.947,0.955,0.964,0.971,0.979,0.989,0.996,1.006,1.014,1.022,1.033,1.048,1.072,1.119])
te_syn = np.array([627.,644.,660.,585.,529.,488.,549.,469.,482.,403.,418.,340.,181.,147.,53.,20.,18.,35.,6.,12.,2.])




def set_params_demo(params=np.array([1.,1.]), **kwargs):
    
    # set chi_e in the PF and SOL region based on input params
    from uedge import bbb, com
    kyepf, kyesol = params
    bbb.kye_use[0:com.ixpt1[0]+1,0:com.iysptrx+1] = kyepf
    bbb.kye_use[com.ixpt2[0]+1:com.nx+2,0:com.iysptrx+1] = kyepf
    bbb.kye_use[:,com.iysptrx+1:com.ny+2] = kyesol



def find_equilibrium_demo(uetools_case: Case, save_dir='.', **kwargs):
    
    # converge the UEDGE case and save the result
    uetools_case.converge(savefname='{}/initial'.format(save_dir), ii1max=200)
    uetools_case.save(savefname='{}/final.hdf5'.format(save_dir))



def loss_function_demo(uetools_case: Case, psi, te, **kwargs):
    
    from uedge import bbb, com
    
    # get psi_n and te from UEDGE
    psi_uedge = (com.psi[bbb.ixmp,:,0] - com.simagxs) / (com.sibdrys - com.simagxs)
    te_uedge = bbb.te[bbb.ixmp,:]/bbb.ev
    
    # interpolate UEDGE data to observed locations
    te_interpolated = np.interp(psi, psi_uedge, te_uedge, left=te_uedge[0], right=te_uedge[-1])
    
    # calculate mean-square-error
    return np.sqrt(np.sum( (te_interpolated-te)**2. / 2. ) / len(te))



def find_constraint_demo():
    
    return None
    
    

def objective_demo(
    uetools_case: Case, 
    params=np.array([1.,1.]), 
    save_dir='./bayes_opt',
    **kwargs
    ):
    
    # create new sub-folder in save_dir
    dir_list = []
    for item in os.listdir(save_dir):
        if os.path.isdir('{}/{}'.format(save_dir, item)): dir_list.append(item)
    sub_dir = '{}/{}'.format(save_dir, len(dir_list))
    os.makedirs(sub_dir)
    
    # set parameter
    set_params_demo(params, **kwargs)
    
    # calculate UEDGE equilibrium
    find_equilibrium_demo(uetools_case, save_dir=sub_dir, **kwargs)
    
    # calculate loss function by comparing UEDGE and observed profiles
    if uetools_case.get('iterm') == 1:
        target = - loss_function_demo(uetools_case, psi = psi_syn, te = te_syn, **kwargs)
    else:
        target = np.nan
    
    # optional: calculate the constraint if needed
    if uetools_case.optimizer.is_constrained:
        constraint = find_constraint_demo()
    else:
        constraint = None
        
    # return params, target, and constraint
    return params, target, constraint
    