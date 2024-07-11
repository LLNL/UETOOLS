"""
A module implements asynchronous Bayesian optimization. The following Python packages are needed:

multiprocessing: 
    Achieve asynchronous sampling. 

sklearn.gaussian_process: 
    Implementation of Gaussian process (GP) and kernels.
    https://github.com/scikit-learn/scikit-learn

bayes_opt: 
    A Bayesian optimization package capable to fitting GP, calculating acquisition function,
    and is compatible with inequality constrain. This package provide slightly more freedom for 
    user to define their GP, kernel, and acquisition functions. 
    https://github.com/bayesian-optimization/BayesianOptimization

    Notice: bayes_opt package tries to maximize target, but we need to minimize loss function. 
    So, in general, we need target = - loss, or - log(loss).
"""

import pickle, os, copy, time
import matplotlib.pyplot as plt
import numpy as np

from uedge import *
from scipy.stats import qmc
from IPython.utils import io

from bayes_opt import BayesianOptimization
from bayes_opt.util import NotUniqueError
from bayes_opt import UtilityFunction
from scipy.optimize import NonlinearConstraint

from sklearn.gaussian_process.kernels import Matern

from .BayesUtility import *
from .BayesDemo import objective_demo


class Bayesian():
    
    def bayes_opt(
        self,
        objective_function=objective_demo, # callable
        constraint=None, # NonlinearConstraint
        param_bounds=np.array([[0.,2.],[0.,2.]]), # np.array
        gp_kernel=None,
        initial_sample=5,
        N_processes=16,
        acq_functions={},
        async_bo_steps=10,
        constant_lier_steps=5,
        bo_plot_status=True,
        save_dir = './bayes_opt',
        **kwargs
        ):
        
        """ 
        An function to do batch BO automatically with default settings. 

        Inputs:
        -------
        objective_function:
            Function defined externally to calculate the UEDGE equilibrium and evaluate a loss function.
        gp_kernel:
            Kernels used in Gaussian regression model.
        initial_sample:
            Number of initial samples used in quasi-Monte Carlo method. The actual sample points is 2^(initial_sample).
        N_processes:
            Number of processes in multiprocess pool.
        acq_functions:
            Acquisition functions used in Bayesian optimization.
        async_bo_steps:
            The number of steps to be calculate in BO process.
        constant_lier_steps:
            The number of data points predicted in each acquisition function using the 'constant lier' approximation.
        bo_plot_status:
            Whether plot currect status in each step of BO. Only work for 2-D case. 
        **kwargs:
            Additional keywords.
        """
        
        self.param_bounds = param_bounds
        self.save_dir = save_dir
        
        # create folder if not exist
        if not os.path.exists(save_dir): os.makedirs(save_dir)
        
        # define Bayesian Optimization object
        self.optimizer = BayesianOptimization(f=None, 
                                              constraint=constraint,
                                              pbounds=array_to_param(param_bounds),
                                              verbose=0, 
                                              random_state=0)
        
        # Set kernels
        self.set_kernel(gp_kernel)
        
        # Calculate initial sampling through quasi-Monte Carlo method
        if self.verbose: print('\n====== Begin initial sampling through quasi-Monte Carlo method =====')
        self.calculate_initial_sampling(objective_function, m=initial_sample, N_process=N_processes, **kwargs)
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    def set_kernel(self, kernel=None):
        """
        Define kernel for Gaussian processes.
        
        Input:
        ------
        kernel: sklearn.gaussian_process.kernels
            The kernel to be used for GP.
        """

        # default kernel
        if kernel == None:
            kernel = Matern(length_scale=1., nu=2.5, length_scale_bounds=(0.2,10))

        # set kernel for GP
        self.optimizer.set_gp_params(kernel=kernel)

        # set kernel for GP to estimate constraint
        if self.optimizer.is_constrained:
            for i in range(len(self.optimizer.constraint._model)):
                self.optimizer.constraint._model[i].set_params(kernel=kernel)
        
        
        
        
        
        
        



    def calculate_initial_sampling(self, objective_function=objective_demo, m=5, N_process=16, **kwargs):
        """
        Search initial sampling asynchronously using quasi-Monte Carlo method (LPtau).
        
        Input:
        ------
        objective_function:
            A user-defined function to calculate equilibrium, target, and constraint (if applicable).
        m: 
            Number of initial sampling in LPtau method, N = 2^m.
        N_process:
            Number of multiprocess pool.
        **kwargs:
            Additional variables
        """

        # calculate initial sampling using LPtau
        sampler = qmc.Sobol(d=len(self.param_bounds), seed=0)
        sample = sampler.random_base2(m=m)
        param_grid = qmc.scale(sample, self.param_bounds[:,0], self.param_bounds[:,1])

        # define multiprocessing pool
        with io.capture_output() as captured:
            pool = QuietPool(processes=N_process)

        # start sampling
        job_status = {}
        for job_id in np.arange(len(param_grid)):
            time.sleep(0.1)
            job_status[job_id] = pool.apply_async(func=objective_function, 
                                                  kwds={'uetools_case': self, 
                                                        'params': param_grid[job_id], 
                                                        'save_dir': self.save_dir, 
                                                        **kwargs}, 
                                                  error_callback=print)

        # wait for all jobs and register data
        for job_id in np.arange(len(param_grid)):
            
            params, target, constraint = job_status[job_id].get()
            
            # TO DO: use a separate variable to restore result so that np.nan result is ignored.
            
            # try to register if not duplicated
            try: self.optimizer.register(params=params, target=target, constraint_value=constraint)                
            except NotUniqueError:
                if self.verbose: print('Unable to register duplicated data. job id = {}'.format(job_id))
                
        # print current best
        if self.verbose:
            print('Current best:')
            current_best = self.optimizer._space.max()
            for key in current_best.keys(): print('{}: {}'.format(key, current_best[key]))

        # close the pool after calculation
        pool.close()
    
    




