import numpy as np
from uetools.UeCase import Case
import os

# Synthetic data to calculate mean-square error between UEDGE and observation
psi_syn = np.array([0.880,0.891,0.901,0.910,0.919,0.930,0.939,0.947,0.955,0.964,
                    0.971,0.979,0.989,0.996,1.006,1.014,1.022,1.033,1.048,1.072,1.119])
te_syn = np.array([627.,644.,660.,585.,529.,488.,549.,469.,482.,403.,418.,340.,181.,147.,53.,20.,18.,35.,6.,12.,2.])



class PhysicsOperationDemo:


    def set_params(self, params=np.array([1.,1.]), **kwargs):
        """ 
        A function acting as the interface between UEDGE and Bayesian optimization (BO). In each step, 
        BO algorithm will provide the next parameter(s) to evaluate. Then one can call this function
        to transform the next parameter(s) into UEDGE setting.
        
        In this demo, input parameter has two elements, representing cross field tarnsport coefficients
        chi_e in private flux and sol regions.
        
        Arguments:
        ----------
        params: np.array(N,) -- Required
            Each parameter BO provide is an np.array with N elements. Generally, they represents some
            transport coefficients that needs to be defined on UEDGE grids. 
        **kwargs:
            Other user defined parameters.
        """
        
        # set chi_e in the PF and SOL region based on input params
        from uedge import bbb, com
        
        kyepf, kyesol = params
        bbb.kye_use[0:com.ixpt1[0]+1,0:com.iysptrx+1] = kyepf
        bbb.kye_use[com.ixpt2[0]+1:com.nx+2,0:com.iysptrx+1] = kyepf
        bbb.kye_use[:,com.iysptrx+1:com.ny+2] = kyesol




    def find_equilibrium(self, uetools_case: Case, save_dir='.', **kwargs):
        """ 
        A function that calculates the UEDGE equilibrium. In the simplest case, one just need to
        call the Case.converge() method.
        
        Arguments:
        ----------
        uetools_case: uetools.Case -- Required
            A uetools Case container to be converged.
        save_dir: String -- Required
            The location where HDF5 files are saved. The default is the current location.
        **kwargs:
            Other user defined parameters.
            
        Return:
        -------
        convergence: boolean -- Required
            An indicator on whether such a case converges. 
        """
        
        # converge the UEDGE case and save the result
        uetools_case.converge(savefname='{}/initial'.format(save_dir), ii1max=200)
        
        # return the convergence
        if uetools_case.get('iterm') == 1: return True
        else: return False




    def loss_function(self, **kwargs):
        """
        A function quantifying the difference between UEDGE calculation and experimental profiles.
        
        In this example, the difference is defined as the mean-square error in the mid-plane Te profile 
        in UEDGE and experiment. The experiment data is given by psi_syn and te_syn.
        
        Argumetns:
        ----------
        **kwargs: 
            Other user defined parameters.
            
        Return:
        -------
        loss: float -- Required
        """
        
        from uedge import bbb, com
        
        # get psi_n and te from UEDGE
        psi_uedge = (com.psi[bbb.ixmp,:,0] - com.simagxs) / (com.sibdrys - com.simagxs)
        te_uedge = bbb.te[bbb.ixmp,:]/bbb.ev
        
        # interpolate UEDGE data to observed locations
        te_interpolated = np.interp(psi_syn, psi_uedge, te_uedge, left=te_uedge[0], right=te_uedge[-1])
        
        # calculate mean-square-error
        return np.sum( (te_interpolated-te_syn)**2. / 2. ) / len(te_syn)




    def find_constraint(self, **kwargs):
        """
        Optional function
        
        A function calculates the constrain of the system. In some situations, one may want to keep some quantities 
        within some range, but do not want to put them into the loss function. Those quantites can be the constraint of 
        the system.
        
        Argumetns:
        ----------
        **kwargs: 
            Other user defined parameters.
        """
        
        return None
        



    def probability_function(self, loss, **kwargs):
        """
        Optional function
        
        A function calculates probability for given loss function. When the loss is the mean-square error, the
        probability is calculated as P = exp(-loss).
        
        Argumetns:
        ----------
        loss: float or np.array
            The value of the loss function
        **kwargs: 
            Other user defined parameters.
            
        Return:
            -------
            probability: float -- Required
        """
        
        return np.exp( -loss )