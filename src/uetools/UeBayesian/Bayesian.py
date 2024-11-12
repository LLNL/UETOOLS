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
from typing import Callable

from bayes_opt import BayesianOptimization
from bayes_opt.util import NotUniqueError
from bayes_opt import UtilityFunction
from scipy.optimize import NonlinearConstraint
from multiprocessing import Pool, TimeoutError


from sklearn.gaussian_process.kernels import Matern

from .BayesUtility import *
from .BayesDemo import *


class Bayesian():
    
    
    
    def __init__(self, case: Case, physics=PhysicsOperationDemo):
        """ 
        Initialize classes to define variables as place holder. 
        
        Arguments:
        ----------
        case: uetools.Case
            The interface between Bayesian object and Case object. It allows Bayesian to handle all UEDGE calculation.
        physics:
            A user-defined class setting the creteria of Bayesian optimization, which should includes the following
            required methods:
        
            set_params(params, **kwargs):
                Defines how given parameters are used in UEDGE calculation.
                
                Arguments:
                ----------
                params: np.array(N,) -- Required
                    Each parameter BO provide is an np.array with N elements. Generally, they represents some
                    transport coefficients that needs to be defined on UEDGE grids. 
                **kwargs:
                    Other user defined parameters.
                
                
            find_equilibrium(uetools_case, save_dir, **kwargs):
                Defines the method to calculate an equilibrium.
                
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


            loss_function(**kwargs):
                Defines how the loss function is defined, which will be minimized during 
                Bayesian optimization process. 
                
                Argumetns:
                ----------
                **kwargs: 
                    Other user defined parameters.
                    
                Return:
                -------
                loss: float -- Required
                
                
            find_constraint(**kwargs): -- Optional
                Define constraint of the system if needed. This will force the algorithm 
                to do a constrained optimization.
                
                Argumetns:
                ----------
                **kwargs: 
                    Other user defined parameters.
                    
                    
            probability_function(**kwargs): -- Optional
                Calculate the probability based for given loss. This is used to estimate the
                uncertainties from Bayesian optimization.
                
                Argumetns:
                ----------
                **kwargs: 
                    Other user defined parameters.
                    
                Return:
                -------
                probability: float -- Required
        """
        
        self.c = case
        self.physics = physics
        self.verbose = self.c.info['verbose']
        
        self.current_job_target = [] # number of calculation that converged
        self.total_job_target = [] # number of calculation that converged
        self.valid_jobs = 0 # number of calculations done during BO that converge and not repeated.
        self.bo_data = {} # dinctionary saving basic information for all finished jobs during BO.
        



    def bayes_opt(
        self,
        param_bounds=np.array([[0.,2.],[0.,2.]]),
        constraint=None,
        gp_kernel=None,
        initial_sample=5,
        N_initial_process=16,
        N_processes=16,
        acq_functions={},
        bo_steps=10,
        constant_lier_steps=5,
        bo_plot_status=True,
        save_dir='./bayes_opt',
        random_state=0,
        timeout=600.,
        **kwargs
        ):
        
        """ 
        An function to do batch BO automatically with default settings. 

        Arguments:
        ----------
        param_bounds: array (N,2)
            Lower and upper bounds for N value of parameters.
        constraint: scipy.optimize.NonlinearConstraint
            The constrains needed in the system. If no constraint is needed, constraint = None.
        gp_kernel:
            Kernels used in Gaussian regression model.
        initial_sample:
            Number of initial samples used in quasi-Monte Carlo method. The actual sample points is 2^(initial_sample).
            If initial_sample == 0, the algorithm will try to read from existing calculations.
        N_initial_process:
            Number of processes in multiprocess pool for initial sampling.
        N_processes:
            Number of processes in multiprocess pool.
        acq_functions:
            Acquisition functions used in Bayesian optimization.
        bo_steps:
            The number of steps to be calculate in BO process.
        constant_lier_steps:
            The number of data points predicted in each acquisition function using the 'constant lier' approximation.
        bo_plot_status:
            Whether plot currect status in each step of BO. Only work for 2-D case. 
        save_dir:
            The directory to save HDF5 for all calculations
        random_state:
            Seed of random number
        **kwargs:
            Additional keywords.
        """
        
        t1 = time.time()
        
        self.param_bounds = param_bounds
        self.save_dir = save_dir
        
        # create folder if not exist
        if not os.path.exists(save_dir): os.makedirs(save_dir)
        
        # define Bayesian Optimization object
        self.optimizer = BayesianOptimization(f=None, 
                                              constraint=constraint,
                                              pbounds=array_to_param(param_bounds),
                                              verbose=0, 
                                              random_state=random_state)
        
        # Set kernels
        self.set_gp_kernel(gp_kernel)
        
        # Calculate initial sampling through quasi-Monte Carlo method
        if initial_sample > 0:
            if self.verbose: print('\n====== Begin initial sampling through quasi-Monte Carlo method =====')
            self.calculate_initial_sampling(m=initial_sample, 
                                            N_process=N_initial_process, 
                                            random_state=random_state, 
                                            timeout=timeout, 
                                            **kwargs)
        else:
            if self.verbose: print('\n====== Begin reading calculated samples =====')
            self.read_existing_samples(save_dir=self.save_dir)
        
        # Define a list of acqusition functions based on initial results
        self.define_acq(acq_functions)
        
        # begin parallel Bayesian optimization
        if self.verbose: print('\n===== Begin Bayesian optimization =====')
        self.batch_BO_searching(N=bo_steps, 
                                step=constant_lier_steps, 
                                N_process=N_processes, 
                                plot_status=bo_plot_status,
                                timeout=timeout,
                                **kwargs)
        
        # print the final result of Bayesian optimization
        print('\n===== Bayesian optimization conclusion =====\n')
        total_time = (time.time() - t1) / 60.
        print('\nTotal BO time =  {:.2f} mins = {:.2f} hours'.format(total_time, total_time/60.))
        
        return self.BO_conclusion()
        
    
    
    
    def bo_objective(
        self,
        params=np.array([1.,1.]), 
        save_case=True,
        save_bo_data=True,
        **kwargs
        ):
        """ 
        The function to be called in parallel during Bayesian optimization. 
        
        Arguments:
        params:
            The parameters that Bayesian optimization works with.
        save_case:
            Whether the save the final HDF5 file.
        save_bo_data:
            Whether to save the bo data in pickle.
        """
        
        # create new sub-folder in save_dir
        dir_list = []
        for item in os.listdir(self.save_dir):
            if os.path.isdir('{}/{}'.format(self.save_dir, item)): dir_list.append(item)
        sub_dir = '{}/{}'.format(self.save_dir, len(dir_list))
        
        try: os.makedirs(sub_dir)
        except: pass
        
        # set parameter
        self.physics.set_params(params, **kwargs)
        
        # calculate UEDGE equilibrium and return the convergence
        try:
            convergence = self.physics.find_equilibrium(uetools_case=self.c, save_dir=sub_dir, **kwargs)
        except:
            convergence = False
        
        # calculate loss function by comparing UEDGE and observed profiles
        if convergence:
            target = - self.physics.loss_function(**kwargs)
        else:
            target = np.nan
        
        # optional: calculate the constraint if needed
        if self.optimizer.is_constrained:
            constraint = self.physics.find_constraint(**kwargs)
        else:
            constraint = None
            
        # save the current case
        if save_case: self.c.save('{}/final.hdf5'.format(sub_dir))
        
        # save the data for Bayesian optimzation 
        if save_bo_data:
            bo_data = {
                'params':params, 
                'target':target, 
                'constraint':constraint, 
                'sub_dir':'{}'.format(len(dir_list))
                }
            with open('{}/bo_data.pickle'.format(sub_dir), 'wb') as handle:
                pickle.dump(bo_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
            
        # return params, target, constraint, and the location of the save file
        return params, target, constraint, sub_dir
        
    
    
    
    def set_gp_kernel(self, kernel=None):
        """
        Define kernel for Gaussian processes.
        
        Arguments:
        ----------
        kernel: sklearn.gaussian_process.kernels
            The kernel to be used for GP.
        """

        # default kernel
        if kernel == None: kernel = Matern(length_scale=1., nu=2.5, length_scale_bounds=(0.2,10))

        # set kernel for GP
        self.optimizer.set_gp_params(kernel=kernel)

        # set kernel for GP to estimate constraint
        if self.optimizer.is_constrained:
            for i in range(len(self.optimizer.constraint._model)):
                self.optimizer.constraint._model[i].set_params(kernel=kernel)




    def calculate_initial_sampling(self, m=5, N_process=16, timeout=600., random_state=0, **kwargs):
        """
        Search initial sampling asynchronously using quasi-Monte Carlo method (LPtau).
        After calculation finished, convergent result will be stored into optimizer.
        
        Arguments:
        ----------
        m: 
            Number of initial sampling in LPtau method, N = 2^m.
        N_process:
            Number of multiprocess pool.
        timeout:
            Max wait time (s) for pool.
        random_state:
            Initial random seed. 
        **kwargs:
            Additional variables
        """

        # calculate initial sampling using LPtau
        sampler = qmc.Sobol(d=len(self.param_bounds), seed=random_state)
        sample = sampler.random_base2(m=m)
        param_grid = qmc.scale(sample, self.param_bounds[:,0], self.param_bounds[:,1])

        # define multiprocessing pool
        with io.capture_output() as captured:
            pool = QuietPool(processes=N_process)

        # start sampling
        job_status = {}
        for job_id in np.arange(len(param_grid)):
            time.sleep(0.5)
            job_status[job_id] = pool.apply_async(func=self.bo_objective, 
                                                  kwds={'params': param_grid[job_id], 
                                                        **kwargs}, 
                                                  error_callback=print)

        # wait for all jobs and register data
        for job_id in np.arange(len(param_grid)):
            
            try:
                params, target, constraint, sub_dir = job_status[job_id].get(timeout=timeout)
            except:
                print('Error getting result for job_id = {}'.format(job_id))
                target = np.nan
                
            self.current_job_target.append(target)
            
            # try to register if result is converge and not duplicated. 
            if not np.isnan(target):
                
                try:                     
                    # save data to optimizer
                    self.optimizer.register(params=params, target=target, constraint_value=constraint)                
                    
                    # save valid job locally
                    self.bo_data[self.valid_jobs] = {'params':params, 'target':target, 
                                                     'constraint':constraint, 'sub_dir':sub_dir}
                    self.valid_jobs += 1
                    
                except NotUniqueError: pass
                
        # print current best
        self.total_job_target.extend(self.current_job_target)
        if self.verbose: self.print_current_status()
        
        # close the pool after calculation
        pool.close()
        self.current_job_target = []




    def read_existing_samples(self, sample_list = None, save_dir = '.'):
        """ 
        Read existing samples that has been calculated save stored.
        
        Arguments:
        sample_list: list
            The list of initial sample to read.
        save_dir:
            The directory to save HDF5 for all calculations
        """
        
        # look for all folders in save_dir
        dir_list = []
        for item in os.listdir(save_dir):
            if os.path.isdir('{}/{}'.format(save_dir, item)): dir_list.append(item)
            
        if len(dir_list) == 0:
            print('Warning: no samples exist! Need to set initial_sample > 0!')
            return None
            
        # read data from calculated cases
        if sample_list is None: sample_list = list(range(len(dir_list)))
        
        for i in sample_list:
            # if folder i is calculated and in the directory
            if '{}'.format(i) in dir_list:
                sub_dir = '{}/{}'.format(save_dir, i)
                
                try:
                    with open('{}/bo_data.pickle'.format(sub_dir), 'rb') as handle:
                        bo_data = pickle.load(handle)
                        
                    params = bo_data['params']
                    target = bo_data['target']
                    constraint = bo_data['constraint']
                    sub_dir = bo_data['sub_dir']
                    
                    self.current_job_target.append(target)
                    
                    # try to register if result is converge and not duplicated. 
                    if not np.isnan(target):
                        
                        try: 
                            # save data to optimizer
                            self.optimizer.register(params=params, target=target, constraint_value=constraint)                
                            
                            # save valid job locally
                            self.bo_data[self.valid_jobs] = bo_data
                            self.valid_jobs += 1
                            
                        except NotUniqueError: pass
                except FileNotFoundError: pass
                
        # print current best
        self.total_job_target.extend(self.current_job_target)
        if self.verbose: self.print_current_status()
        self.current_job_target = []
        



    def define_acq(self, acq_functions: dict = {}):
        """
        Initialize the acquisition functions used in async searching process.
        The default version can only be defined after initial sampling.

        Arguments:
        ----------
        acq_functions: dict
            Dictionary whose elements are bayes_opt.UtilityFunction object.
        """

        # Define acqusition functions if no input
        if len(acq_functions.keys()) == 0:

            E_max = np.abs(self.optimizer._space.target).max()

            # Expected improvement (EI)
            acq_functions[0] = UtilityFunction(kind="ei", xi=0.)
            acq_functions[1] = UtilityFunction(kind="ei", xi=-0.2*E_max)
            acq_functions[2] = UtilityFunction(kind="ei", xi=0.1*E_max)
            acq_functions[3] = UtilityFunction(kind="ei", xi=0.2*E_max)

            # Upper Confidence Bound (UCB)
            acq_functions[4] = UtilityFunction(kind="ucb", kappa=1.)
            acq_functions[5] = UtilityFunction(kind="ucb", kappa=4.)
            # acq_functions[6] = UtilityFunction(kind="ucb", kappa=4.)
            # acq_functions[7] = UtilityFunction(kind="ucb", kappa=16.)
        
        self.acq_functions = acq_functions
        self.next_point = self.optimizer.suggest(self.acq_functions[0])    
    
    
    
    
    def multi_suggest(self, acq_function, step=4, return_optimizer=False):
        """
        For a given obervation data and acquisition function, predict #step points to be evaluate 
        based on the "constant liar" approximation (https://doi.org/10.1007/978-3-642-10701-6_6).

        Traditionaly, Bayesian optimization fit N observed point with GP. With given acquisition function,
        it predict the next point to evaluate. The constant liar simply assume the N+1 point will return 
        the value exactly the same as GP predict, then use the acquisition function to predict the N+2 point 
        to evaluate. Same procedure can be repeated multiple times.

        Since the "constant liar" approximation is a sequential job, it can not be parallized.

        Arguments:
        ----------
        acq_function: bayes_opt.UtilityFunction
            Acquisition functions used in this calculation. 
        step:
            The number of points predicted from the "constant liar" approximation.
        return_optimizer:
            Whether return the optimizer (bayes_opt.BayesianOptimization) after finish.
        """

        suggest_points = {}
        tmp_optimizer = copy.deepcopy(self.optimizer)

        # calculate #asyc_step steps for each acquisition functions
        for i in range(step):

            # get prediction for next step
            suggest_points[i] = tmp_optimizer.suggest(acq_function)
            x = param_to_array(suggest_points[i]).reshape(1,-1)

            # assume the result is exactly what GP could predict
            approx_target = tmp_optimizer._gp.predict(x)[0]

            # if we have constrain, also predict the value of constrains
            if tmp_optimizer.is_constrained:
                approx_constraint = tmp_optimizer.constraint.approx(x)[0]
            else:
                approx_constraint = None

            # register approximated values
            try:
                tmp_optimizer.register(params = suggest_points[i], 
                            target = approx_target, 
                            constraint_value = approx_constraint)
            
            # if multiple points are same, they will not be able to be registered. 
            except NotUniqueError:
                if self.verbose: print('Unable to register due to duplicated points during multi_suggest.')

        if return_optimizer:
            return suggest_points, tmp_optimizer
        else:
            return suggest_points




    def async_suggest(self, step=4):
        """
        Calculate the suggested points for multiple acquisition functions asynchronously.

        Arguments:
        ----------
        step:
            The number of points predicted from the "constant liar" approximation. 
        """

        # start a multiprocessing pool
        with io.capture_output() as captured:
            pool = QuietPool(processes = len(self.acq_functions))

        # async suggesting
        acq_suggest_jobs = {}
        acq_suggest = {}

        for i in range(len(self.acq_functions)):
            acq_suggest_jobs[i] = pool.apply_async(self.multi_suggest, 
                                                   args=(self.acq_functions[i], step), 
                                                   error_callback=print)
            
        # collect results
        for i in range(len(self.acq_functions)):
            tmp_suggests = acq_suggest_jobs[i].get()

            for j in range(step):
                acq_suggest[i*step + j] = tmp_suggests[j]

        # close the pool
        pool.close()

        return acq_suggest




    def find_non_repeat(self, data, N=16, d_min=0.05):
        """
        Give a collection of data, find #N non-repeated data points where the distance of each points 
        are larger than d_min.

        Suggestions from different acquisition functions may be very similar. Therefore, only the 
        non-repeated data points need to be evaluated. If the number of non-repeated points is less 
        than N, add some random points to evaluate.

        Arguments:
        ----------
        data: dict
            Data points given by async_suggest().
        N:
            Number of non-repeated point to be collected.
        d_min:
            Minimum (Euclidean) distance between non-repeated points.
        random_replenish:
            When the number of non-repeated points is less than N, whether to use random sample points
            to fill in the slots.
        """

        non_repeat_data = {}
        non_repeat_data[0] = data[0]

        for i in range(len(data.keys())-1):

            # collect non-repeat data until enough data is collected
            if len(non_repeat_data.keys()) < N:
                repeated = False

                for j in range(len(non_repeat_data)):
                    distance = 0.
                    
                    # calculate the distance between points
                    for param in non_repeat_data[j]:
                        distance += (non_repeat_data[j][param] - data[i][param])**2.
                    if distance < d_min**2.: repeated = True

                if not repeated:
                    non_repeat_data[len(non_repeat_data.keys())] = data[i]

        if len(non_repeat_data.keys()) < N:

            if self.verbose:
                print('Warning: not enough non-repeated data, missing {} points. Add random points instead.'.format(N-len(non_repeat_data.keys())))

            # add random data 
            dim = len(self.param_bounds)

            while len(non_repeat_data.keys()) < N:
                random_data = [np.random.uniform(self.param_bounds[i,0], self.param_bounds[i,1]) for i in range(dim)]
                non_repeat_data[len(non_repeat_data.keys())] = array_to_param(random_data)

        return non_repeat_data




    def batch_BO_searching(self, N=20, step=4, N_process=16, plot_status=True, timeout=600., **kwargs):
        """
        Proceed to Bayesian Optimization for #N steps. In every step, run #step number of points for each 
        acquisition function and collected #N_process of non-repeated points. All #N_process points are
        run in parallel using multiprocessing. The data are saved in self.save_dir.

        Strictly speaking, the current version only uses pool.apply_async, but it is not exactly asynchronous
        searching because it waited all job to finish before calculating the next step. Future version need to 
        modify it and make it truly asynchronous.

        Arguments:
        ----------
        N:
            Total number of steps in Bayesian optimization.
        step:
            Number of points generated from each acquisition function using the the "constant liar" 
            approximation.
        N_process:
            Number of processes in multiprocessing pool. Same number of non-repeated points should be
            generated at each step.
        plot_status:
            Whether to plot the currect status of BO. Only work for 2-D searching.
        timeout:
            Abort if UEDGE does not converge within time.
        **kwargs:
            Additional arguments
        """

        # async searching
        for i in range(N):

            if self.verbose: print('\nBO step {} ====='.format(i+1))

            # record the time in each step
            t1 = time.time()

            # dict to save the job data and record the job id
            job_list = {}
            total_suggeested = self.async_suggest(step=step)
            next_points = self.find_non_repeat(total_suggeested, N=N_process, d_min=0.05)
            
            t2 = time.time()
            # print the time used
            if self.verbose:
                print('\nTime for suggestion: {:.2f} mins'.format((t2-t1)/60.))
            

            # plot status only for 2-D inference
            if plot_status and (self.optimizer.space.bounds.shape[0]==2):
                plot_async_2D(param_bounds=self.param_bounds, 
                              optimizer=self.optimizer, 
                              acq_function = self.acq_functions[0],
                              next_points=next_points, 
                              title='Step = {}'.format(i+1),
                              save_folder=self.save_dir,
                              **kwargs)
                
            # start a multiprocessing pool
            with io.capture_output() as captured:
                pool = QuietPool(processes = N_process)

            # submit all jobs
            for asyc_job_id in range(len(next_points)):
                time.sleep(0.5)
                job_list[asyc_job_id] = pool.apply_async(self.bo_objective, 
                                                         kwds={'params': param_to_array(next_points[asyc_job_id]), 
                                                               **kwargs},
                                                         error_callback=print)

            # wait for all jobs and register data
            for asyc_job_id in range(len(next_points)):
                
                try:
                    params, target, constraint, sub_dir = job_list[asyc_job_id].get(timeout=timeout)
                except:
                    print('Error getting result for job_id = {}'.format(asyc_job_id))
                    target = np.nan
                    
                self.current_job_target.append(target)
                
                # try to register if result is converge and not duplicated. 
                if not np.isnan(target):
                    
                    try: 
                        # save data to optimizer
                        self.optimizer.register(params=params, target=target, constraint_value=constraint)                
                        
                        # save valid job locally
                        self.bo_data[self.valid_jobs] = {'params':params, 'target':target, 
                                                        'constraint':constraint, 'sub_dir':sub_dir}
                        self.valid_jobs += 1
                        
                    except NotUniqueError: pass
                    
            # print current best
            self.total_job_target.extend(self.current_job_target)
            if self.verbose: self.print_current_status()

            # just to let optimizer fit data
            tmp = self.optimizer.suggest(self.acq_functions[0])

            # print the time used
            if self.verbose:
                print('\nTime for UEDGE: {:.2f} mins'.format((time.time()-t2)/60.))

            # close the pool after calculation
            pool.close()
            self.current_job_target = []




    def print_current_status(self):
        """ 
        Print current best fit that has the smallest loss
        """
        
        print('\nCurrent best:')
        current_best = self.optimizer._space.max()
        for key in current_best.keys(): 
            print('    {}: {}'.format(key, current_best[key]))
        
        # Convergence analysis
        print('\nConvergence analysis:')
        
        total_job = len(self.current_job_target)
        unconverged_job = np.isnan(np.array(self.current_job_target)).sum()
        converged_job = total_job-unconverged_job
        print('    Current: {:.1f} % jobs ({}/{}) converged.'.format(100.*converged_job/total_job, converged_job, total_job))
        
        total_job = len(self.total_job_target)
        unconverged_job = np.isnan(np.array(self.total_job_target)).sum()
        converged_job = total_job-unconverged_job
        print('    Total: {:.1f} % jobs ({}/{}) converged.'.format(100.*converged_job/total_job, converged_job, total_job))
        
        # GP kernel
        print('\nCurrent GP kernel:')
        try:
            print('    {}'.format(self.optimizer._gp.kernel_))
        except:
            self.optimizer.suggest(UtilityFunction(kind="ei", xi=0.))
            print('    {}'.format(self.optimizer._gp.kernel_))
            




    def BO_conclusion(self, verbose = True):
        """
        Return the conclusion of the Bayesian optimization, which includes:
        1. The minimum loss founded and the paramters associated;
        2. The directory of the HDF5 file asscociated with the best parameter.
        """

        # 1. Find max target
        target_max = self.optimizer._space.max()['target']
        ind = 0
        
        while self.bo_data[ind]['target'] != target_max:
            ind += 1
            
        # the container of conclusion
        conclusion = self.bo_data[ind]
        
        if getattr(self.physics, 'probability_function', None) is not None:
            conclusion['uncertainty'] = self.get_uncertainty()
        
        if verbose:
            for key in conclusion.keys(): 
                print('{}: {}'.format(key, conclusion[key]))

        return conclusion
    
    
    
    
    def get_uncertainty(self):
        """ 
        Return the uncertainty during the BO processes when self.probability_function is defined.
        """
        
        if getattr(self.physics, 'probability_function', None) is not None:
            
            if self.optimizer.is_constrained:
                
                # get data that within constraint
                ind = np.where(np.array([r['allowed'] for r in self.optimizer.res])==True)[0]
                valid_params = self.optimizer._space.params[ind,:]
                valid_loss = - self.optimizer._space.target[ind]
                
            else:
                
                # when no constraint, use all data
                valid_params = self.optimizer._space.params
                valid_loss = - self.optimizer._space.target
                
            # calculate probability
            P = self.physics.probability_function(valid_loss).reshape(-1,1)
            
            # calculate uncertainties
            params_max = param_to_array(self.optimizer._space.max()['params']).reshape(1,-1)
            params_uncertainty = np.sqrt( np.sum( (valid_params-params_max)**2. * P, axis=0) / np.sum(P) )
            
            return params_uncertainty
