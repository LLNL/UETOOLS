"""
A module implements asynchronous Bayesian optimization. The following Python packages are needed:

multiprocessing: 
    Achieve asynchronous sampling. 

sklearn.gaussian_process: 
    Implementation of Gaussian process (GP) and kernels.
    https://github.com/scikit-learn/scikit-learn

bayes_opt: 
    A Bayesian optimization package capable to fitting GP, calculating acquisition function. 
    This package provide slightly more freedom for user to define their GP, kernel, and 
    acquisition functions. 
    https://github.com/bayesian-optimization/BayesianOptimization

    Notice: bayes_opt package tries to maximize target, but we need to minimize loss function. 
    So, in general, we need target = - loss, or - log(loss).
"""

import pickle, os, copy, time
import numpy as np

from uedge import *
from scipy.stats import qmc
from IPython.utils import io

from bayes_opt import BayesianOptimization
from bayes_opt.acquisition import ExpectedImprovement, UpperConfidenceBound, ProbabilityOfImprovement, ConstantLiar
from bayes_opt.util import NotUniqueError

from sklearn.gaussian_process.kernels import Matern
from sklearn.gaussian_process import GaussianProcessRegressor
from bayes_opt import TargetSpace

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
                params: dict -- Required
                    Each parameter BO provide is a dict with N elements. Generally, they represents some
                    transport coefficients that needs to be defined on UEDGE grids. The key words should be
                    {'x0':1., 'x1':1., ...}
                **kwargs:
                    Other user defined parameters.
                
                
            find_equilibrium(case, save_dir, **kwargs):
                Defines the method to calculate an equilibrium.
                
                Arguments:
                ----------
                case: uetools.Case -- Required
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
        
        self.get = case.get # Links the subclass to the case
        self.c = case
        self.physics = physics
        self.verbose = self.c.info['verbose']
        
        self.current_job_target = [] # number of calculation that converged
        self.total_job_target = [] # number of calculation that converged
        self.valid_jobs = 0 # number of calculations done during BO that converge and not repeated.
        self.bo_data = {} # dinctionary saving basic information for all finished jobs during BO.
        
        np.set_printoptions(precision=3)
        



    def bayes_opt(
        self,
        param_bounds={'x0': (0., 2.), 'x1': (0., 2.)},
        gp_kernel=None,
        initial_sample=5,
        N_initial_process=16,
        N_processes=16,
        acq_functions={},
        bo_steps=10,
        constant_lier_steps=4,
        constant_lier_strategy='max',
        bo_plot_status=False,
        save_dir='./bayes_opt',
        random_state=0,
        timeout=600.,
        d_min_start=1.e-3,
        d_min_decay=3, 
        d_min_lim=1.e-3,
        **kwargs
        ):
        
        """ 
        An function to do batch BO automatically with default settings. 

        Arguments:
        ----------
        param_bounds: dict
            Lower and upper bounds for N value of parameters.
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
        constant_lier_strategy:
            The strategy in the constant lier approximation
        bo_plot_status:
            Whether plot currect status in each step of BO. Only work for 2-D case. 
        save_dir:
            The directory to save HDF5 for all calculations
        random_state:
            Seed of random number
        timeout:
            Max time to run for each async job
        d_min_start:
            d_min controls the minimum distance between batch samples. d_min_start is the initial distance
        d_min_decay:
            After a certain steps, decrease d_min by 2.
        d_min_lim:
            The minimum d_min allowed.
        **kwargs:
            Additional keywords.
        """
        
        t1 = time.time()
        
        self.param_bounds = param_to_array(param_bounds)
        self.save_dir = save_dir
        
        # create folder if not exist
        if not os.path.exists(save_dir): os.makedirs(save_dir)
        
        # define Bayesian Optimization object
        self.optimizer = BayesianOptimization(f=None, 
                                              pbounds=param_bounds,
                                              verbose=0, 
                                              random_state=random_state,
                                              allow_duplicate_points=True)
        
        # Set kernels
        self.set_gp_kernel(gp_kernel)
        
        # Calculate initial sampling through quasi-Monte Carlo method
        if initial_sample > 0:
            if self.verbose: print('\n============ Begin initial sampling through quasi-Monte Carlo method ===========\n')
            self.calculate_initial_sampling(m=initial_sample, 
                                            N_process=N_initial_process, 
                                            random_state=random_state, 
                                            timeout=timeout+600, 
                                            **kwargs)
        else:
            if self.verbose: print('\n============ Begin reading calculated samples ===========')
            self.read_existing_samples(save_dir=self.save_dir)
        
        # Define a list of acqusition functions based on initial results
        self.define_acq(acq_functions)
        
        # begin parallel Bayesian optimization
        if self.verbose: print('\n=========== Begin Bayesian optimization ===========')
        self.batch_BO_searching(N=bo_steps, 
                                constant_lier_steps=constant_lier_steps,
                                constant_lier_strategy=constant_lier_strategy,
                                N_process=N_processes, 
                                plot_status=bo_plot_status,
                                timeout=timeout,
                                d_min_start=d_min_start,
                                d_min_decay=d_min_decay, 
                                d_min_lim=d_min_lim,
                                **kwargs)
        
        # print the final result of Bayesian optimization
        print('\n=========== Bayesian optimization conclusion ===========\n')
        total_time = (time.time() - t1) / 60.
        print('\nTotal BO time =  {:.2f} mins = {:.2f} hours'.format(total_time, total_time/60.))
        
        return self.BO_conclusion()
        
    
    
    
    def bo_objective(
        self,
        params={'x0':1., 'x1':1.}, 
        save_case=True,
        save_bo_data=True,
        **kwargs
        ):
        """ 
        The function to be called in parallel during batched Bayesian optimization. 
        
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
        
        # initialize save file
        self.c.restore_save(self.c.info['savefile'])
        
        # set parameter
        self.physics.set_params(params, **kwargs)
        
        # calculate UEDGE equilibrium and return the convergence
        try:
            convergence = self.physics.find_equilibrium(case=self.c, save_dir=sub_dir, **kwargs)
        except:
            convergence = False
        
        # calculate loss function by comparing UEDGE and observed profiles
        if convergence:
            target = - self.physics.loss_function(**kwargs)
        else:
            target = np.nan
            
        # save the current case
        if save_case: self.c.save('{}/final.hdf5'.format(sub_dir))
        
        # save the data for Bayesian optimzation 
        if save_bo_data:
            bo_data = {
                'params':params, 
                'target':target, 
                'sub_dir':'{}'.format(len(dir_list))
                }
            with open('{}/bo_data.pickle'.format(sub_dir), 'wb') as handle:
                pickle.dump(bo_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
            
        # return params, target, and the location of the save file
        return params, target, sub_dir
        
    
    
    
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
        
        # record the time in each step
        t1 = time.time()

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
            time.sleep(1.)
            job_status[job_id] = pool.apply_async(func=self.bo_objective, 
                                                  kwds={'params': array_to_param(param_grid[job_id]), 
                                                        **kwargs}, 
                                                  error_callback=print)

        # wait for all jobs and register data
        time_left = timeout
        for job_id in np.arange(len(param_grid)):
            
            t_async = time.time()
            
            try:
                params, target, sub_dir = job_status[job_id].get(timeout=time_left)
            except:
                print('Error (time out) getting result for job_id = {}, params = {}'.format(job_id, param_grid[job_id]))
                target = np.nan
                
            time_left -= time.time() - t_async
            if time_left <= 0.: time_left = 5.
                
            self.current_job_target.append(target)
            
            # try to register if result is converge and not duplicated. 
            if not np.isnan(target):
                
                try:                     
                    # save data to optimizer
                    self.optimizer.register(params=params, target=target)                
                    
                    # save valid job locally
                    self.bo_data[self.valid_jobs] = {'params':params, 
                                                     'target':target, 
                                                     'sub_dir':sub_dir}
                    self.valid_jobs += 1
                    
                except NotUniqueError: pass
                
        # print current best
        self.total_job_target.extend(self.current_job_target)
        if self.verbose: self.print_current_status()
        
        t2 = time.time()
        if self.verbose: print('\nTime for initial sampling = {:.2f} mins.'.format((t2-t1)/60.))
        
        # close the pool after calculation
        pool.close()
        pool.terminate()
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
                    sub_dir = bo_data['sub_dir']
                    
                    self.current_job_target.append(target)
                    
                    # try to register if result is converge and not duplicated. 
                    if not np.isnan(target):
                        
                        try: 
                            # save data to optimizer
                            self.optimizer.register(params=params, target=target)                
                            
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
            Dictionary whose elements belongs to bayes_opt.acquisition
        """

        # Define acqusition functions if no input
        if len(acq_functions.keys()) == 0:

            E_mean = np.abs(self.optimizer._space.target).mean()
            
            acq_functions[0] = ExpectedImprovement(xi=0.)
            acq_functions[1] = UpperConfidenceBound()
            acq_functions[2] = ProbabilityOfImprovement(xi=0.)
            acq_functions[3] = ExpectedImprovement(xi=0.02*E_mean)
            acq_functions[4] = ExpectedImprovement(xi=1.*E_mean)
        
        self.acq_functions = acq_functions
        self.next_point = self.optimizer.suggest()    
    



    def async_suggest(self, step=4, strategy='max'):
        """
        Calculate the suggested points for multiple acquisition functions asynchronously.

        Arguments:
        ----------
        step:
            The number of points predicted from the "constant liar" approximation. 
        strategy: 
            The strategy of constant lier approximation. See documents in bayes_opt.acquisition.ConstantLiar.
            Possible input:
            'min', 'max', 'mean': set one single strategy
            ['min', 'max', ...]: set multiple strategies
            
        """
        
        if isinstance(strategy, str):
            strategy = [strategy]
            
        N_acq = len(self.acq_functions)
        N_strat = len(strategy)

        # start a multiprocessing pool
        with io.capture_output() as captured:
            pool = QuietPool(processes = N_acq * N_strat)

        # async suggesting
        acq_suggest_jobs = {}
        acq_suggest = {}

        for i in range(N_acq):
            for j in range(N_strat):
                                
                acq_suggest_jobs[i*N_strat + j] = pool.apply_async(multi_suggest, 
                                                                 args=(self.optimizer, 
                                                                       self.acq_functions[i],
                                                                       step,
                                                                       strategy[j]), 
                                                                 error_callback=print)
            
        # collect results
        for i in range(N_acq):
            for j in range(N_strat):
                
                try:
                    tmp_suggests = acq_suggest_jobs[i*N_strat+j].get(timeout=300)
                    
                    for k in range(step):
                        acq_suggest[(i*N_strat+j)*step + k] = tmp_suggests[k]
                
                except:
                    print('Error getting async suggestion #{}'.format(i*N_strat+j))

        # close the pool
        pool.close()
        pool.terminate()

        return acq_suggest




    def find_non_repeat(self, data, N=16, d_min=0.05, strategy='max', print_sample_detail=True):
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
        strategy: 
            The strategy of constant lier approximation
        print_sample_detail:
            Show the sample belongs to which acq as (#acq, #aync_step).
        """
        
        if isinstance(strategy, str):
            strategy = [strategy]
            
        N_acq = len(self.acq_functions)
        N_strat = len(strategy)

        non_repeat_data = {}
        non_repeat_data[0] = data[0]
        
        data_used = [0]
        async_step = int(len(data)/(N_acq*N_strat))

        for i in range(1, len(data.keys())):

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
                    data_used.append(i)
            
        # show which sample from whcih acq        
        if print_sample_detail:
            print('\nSampling with minimum distance = {:.2f}:'.format(d_min))
            for i in range(len(data_used)): 
                
                step_number = int(data_used[i]%async_step)
                strategy_number = int(((data_used[i] - step_number)/async_step)%N_strat)
                acq_number = int(((data_used[i] - step_number)/async_step - strategy_number)/N_strat)
                
                print('Acq #{}, CL_step #{}, CL_strategy = {}, val = ('.format(acq_number, step_number, strategy[strategy_number]), end='')
                
                for key in non_repeat_data[i].keys(): print('{: 7.3f},'.format(non_repeat_data[i][key]), end='')
                print(');')
            

        if len(non_repeat_data.keys()) < N:

            if self.verbose:
                print('Warning: not enough non-repeated data, missing {} points. Add random points instead.'.format(N-len(non_repeat_data.keys())))

            # add random data 
            dim = len(self.param_bounds)

            while len(non_repeat_data.keys()) < N:
                random_data = [np.random.uniform(self.param_bounds[i,0], self.param_bounds[i,1]) for i in range(dim)]
                non_repeat_data[len(non_repeat_data.keys())] = array_to_param(random_data)

        return non_repeat_data




    def batch_BO_searching(self, 
                           N=20, 
                           constant_lier_steps=4,
                           constant_lier_strategy='max',
                           N_process=16, 
                           plot_status=True, 
                           timeout=600.,
                           d_min_start=1.e-3,
                           d_min_decay=3, 
                           d_min_lim=1.e-3,
                           non_converge_replacement='mean',
                           **kwargs):
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
        constant_lier_steps:
            Number of points generated from each acquisition function using the the "constant liar" 
            approximation.
        constant_lier_strategy:
            Strategy of constant lier approximation. Can use a float number, or 'min', 'max', 'mean'.
        N_process:
            Number of processes in multiprocessing pool. Same number of non-repeated points should be
            generated at each step.
        plot_status:
            Whether to plot the currect status of BO. Only work for 2-D searching.
        timeout:
            Abort if UEDGE does not converge within time.
        d_min_start:
            d_min controls the minimum distance between batch samples. d_min_start is the initial distance
        d_min_decay:
            After a certain steps, decrease d_min by 2.
        d_min_lim:
            The minimum d_min allowed.
        non_converge_replacement:
            If the code does not converge, replace the target with some value.
            'mean' / 'min' / 'max': mean / min / max of previous targets
        **kwargs:
            Additional arguments
        """
        
        # distence between difference samples decreases as BO precess
        d_min = d_min_start

        # async searching
        for i in range(N):

            if self.verbose: print('\n\n========== BO step {} =========='.format(i+1))

            # record the time in each step
            t1 = time.time()
            
            # reduce d_min if set
            if (i>0) and (i%d_min_decay==0): d_min /= 2.
            if d_min < d_min_lim: d_min = d_min_lim

            # dict to save the job data and record the job id
            job_list = {}
            total_suggeested = self.async_suggest(step=constant_lier_steps, strategy=constant_lier_strategy)
            next_points = self.find_non_repeat(total_suggeested, N=N_process, strategy=constant_lier_strategy, d_min=d_min)
            
            t2 = time.time()
            # print the time used
            if self.verbose:
                print('\nTime for suggestion: {:.1f} seconds\n'.format(t2-t1))
            

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
                time.sleep(1.)
                job_list[asyc_job_id] = pool.apply_async(self.bo_objective, 
                                                         kwds={'params': next_points[asyc_job_id],
                                                               **kwargs},
                                                         error_callback=print)

            # wait for all jobs and register data
            time_left = timeout
            for asyc_job_id in range(len(next_points)):
                
                t1 = time.time()
                
                try:
                    params, target, sub_dir = job_list[asyc_job_id].get(timeout=time_left)
                except:
                    print('Error (time out) getting result for job_id = {}, params = {}'.format(asyc_job_id, param_to_array(next_points[asyc_job_id])))
                    target = np.nan
                    params = next_points[asyc_job_id]
                    
                time_left -= time.time() - t1
                if time_left <= 0.: time_left = 5.
                    
                self.current_job_target.append(target)
                
                # try to register if result is converge and not duplicated. 
                if not np.isnan(target):
                    _target = target
                # If not converged, replance the target with something else
                elif non_converge_replacement == 'mean':
                    _target = self.optimizer.space.target.mean()
                    
                try: 
                    # save data to optimizer
                    self.optimizer.register(params=params, target=_target)                
                    
                    # save valid job locally
                    self.bo_data[self.valid_jobs] = {'params':params, 
                                                     'target':target, 
                                                     'sub_dir':sub_dir}
                    self.valid_jobs += 1
                    
                except NotUniqueError: pass
                    
                
                    
            # print current best
            self.total_job_target.extend(self.current_job_target)
            if self.verbose: self.print_current_status()

            # just to let optimizer fit data
            tmp = self.optimizer.suggest()

            # print the time used
            if self.verbose:
                print('\nUEDGE time: {:.2f} mins'.format((time.time()-t2)/60.))

            # close the pool after calculation
            pool.close()
            pool.terminate()
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
            self.optimizer.suggest()
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
            
            valid_params = self.optimizer._space.params
            valid_loss = - self.optimizer._space.target
                
            # calculate probability
            P = self.physics.probability_function(valid_loss).reshape(-1,1)
            
            # calculate uncertainties
            params_max = param_to_array(self.optimizer._space.max()['params']).reshape(1,-1)
            params_uncertainty = np.sqrt( np.sum( (valid_params-params_max)**2. * P, axis=0) / np.sum(P) )
            
            return params_uncertainty


    
    
    
    
def multi_suggest(optimizer: BayesianOptimization, acq_function, step=4, strategy='max'):
    """ 
    Generate multiple suggestion points in one step for paralleled evaluation. Multiple points are generated
    using the bayes_opt.acquisition.ConstantLiar function based based on the "constant liar" approximation 
    (https://doi.org/10.1007/978-3-642-10701-6_6).
    
    Arguments:
    ----------
    optimizer: 
        The BayesianOptimization objective that contains all information
    acq_function:
        Acquisition function used in multi suggestion
    step:
    step:
        The number of points predicted from the "constant liar" approximation.
    """
    
    suggest_points = {}
    gp = copy.deepcopy(optimizer._gp)
    space = copy.deepcopy(optimizer.space)
    CL_acq = ConstantLiar(base_acquisition=acq_function, strategy=strategy)
    
    total_suggestion_time = 0.

    # calculate #asyc_step steps for each acquisition functions
    for i in range(step):
        
        t1 = time.time()

        # get prediction for next step
        suggest_points[i] = array_to_param(CL_acq.suggest(gp=gp, target_space=space, n_random=100000))
            
        one_step_time = time.time() - t1
        total_suggestion_time += one_step_time
        print('step = {}, time = {:.2f} s'.format(i, one_step_time))
        
    print('Total suggestion time = {:.1f} s.'.format(total_suggestion_time))

    return suggest_points