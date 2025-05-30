import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import os


def param_to_array(param):
    """
    Auxiliary function. Convert dict to numpy array.
    parameter: dict
        Example: 
        Parameter values: {'x0':0, 'x1':1} >> np.array([0,1])
        Paramter range: {'x0':[0,1], 'x1':[2,3]} >> np.array([[0,1],[2,3]])
    """
    
    keys = list(param.keys())
    x0 = np.reshape(np.array(param[keys[0]]), -1)
    x = np.zeros((len(keys), len(x0)))

    for i in range(len(keys)):
        x[i] = param['x{}'.format(i)]
    
    return np.squeeze(x)



def array_to_param(x):
    """
    Auxiliary function. Convert np.array to dict. 
    x: np.array
        Example: 
        Parameter values: np.array([0,1]) >> {'x0':0, 'x1':1}
        Paramter range: np.array([[0,1],[2,3]]) >> {'x0':[0,1], 'x1':[2,3]}
    """

    param = {}
    for i in range(len(x)):
        param['x{}'.format(i)] = x[i]

    return param



def plot_async_2D(param_bounds, optimizer, acq_function, next_points=None, title='', save_folder=None, show_fig=True, **kwargs):
    """
    Plot the GP mean estimation and one of the acquisition functions.
    White dots: observed data
    Pink dots: data to be observe in this step
    Red dots: data that does not satisfy the constraint
    """
    
    x0_lim = param_bounds[0]
    x1_lim = param_bounds[1]

    x = np.linspace(x0_lim[0], x0_lim[1], 200)
    y = np.linspace(x1_lim[0], x1_lim[1], 200)
    xy = np.array([[x_i, y_j] for y_j in y for x_i in x])
    X, Y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(13,5))
    grid = fig.add_gridspec(ncols=2, nrows=1)

    fig.subplots_adjust(wspace=0.25, hspace=0.3)
    fig.suptitle(title, fontsize=20)

    # Estimation of target
    ax0 = plt.subplot(grid[0,0])
    Z_est = optimizer._gp.predict(xy).reshape(X.shape)
    plt.contourf(X, Y, Z_est, levels=20)
    plt.title('Target');
    plt.colorbar();
    plt.xlabel('x0', fontsize=12)
    plt.ylabel('x1', fontsize=12)

    # Estimate acq
    ax1 = plt.subplot(grid[0,1])
    gp_mean, gp_std = optimizer._gp.predict(xy, return_std=True)
    acq_function.y_max = optimizer.space.target.max()
    acq_est = acq_function.base_acq(gp_mean, gp_std).reshape(X.shape)
    plt.contourf(X, Y, acq_est, levels=20)
    plt.title('Acquisition #0')
    plt.colorbar()
    plt.xlabel('x0', fontsize=12)
    plt.ylabel('x1', fontsize=12)

    # show all data points
    res = optimizer.res
    keys = list(res[0]['params'].keys())
    x_ = np.array([r["params"][keys[0]] for r in res])
    y_ = np.array([r["params"][keys[1]] for r in res])
    
    if optimizer.is_constrained:
        a_ = np.array([r["allowed"] for r in res])
    else:
        a_ = np.array([True for r in res])

    for ax in [ax0, ax1]:

        ax.scatter(x_[a_], y_[a_], c='white', s=20, edgecolors='black')
        ax.scatter(x_[~a_], y_[~a_], c='red', s=20, edgecolors='black')
        
        xx, yy = param_to_array(optimizer._space.max()['params'])
        ax.scatter([xx], [yy], marker='*', c='orange', s=200, edgecolors='black', zorder=10)

        if next_points is not None:
            for i in next_points.keys():
                next_point = param_to_array(next_points[i])
                ax.scatter(next_point[0], next_point[1], marker='D', c='magenta', s=20, edgecolors='black')
            
    plt.scatter([],[], c='white', s=20, edgecolors='black', label='Data calculated')
    plt.scatter([], [], marker='*', c='orange', s=100, edgecolors='black', label='Current best')
    if next_points is not None:
        plt.scatter([],[], marker='D', c='magenta', s=20, edgecolors='black', label='Next evaluations')
    if optimizer.is_constrained:
        plt.scatter([],[], c='red', s=20, edgecolors='black')
    ax1.legend(bbox_to_anchor=(1.85, 0.6))
    
    if save_folder is not None: plt.savefig('{}/{}.pdf'.format(save_folder, title), bbox_inches='tight')
        
    if not show_fig: plt.close()



def mute():
    """
    Redirect output to /dev/null by modifying the C library
    records. This should also redirect output from Fortran code.

    Taken from https://stackoverflow.com/a/978264

    To use this, create a pool and pass this function as the initializer:

        import multiprocessing as mp
        pool = mp.Pool(processes=2, initializer=mute)

    """
    # open 2 fds
    null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
    # save the current file descriptors to a tuple
    save = os.dup(1), os.dup(2)
    # put /dev/null fds on 1 and 2
    os.dup2(null_fds[0], 1)
    os.dup2(null_fds[1], 2)



def QuietPool(processes: int = 2, *args, **kwargs):
    """
    Create a multiprocessing pool that suppresses output.
    Adds an initializer argument to redirect stdout to /dev/null.

    processes: int
        Maximum number of processes in the pool, has to be less than max CPU number
    """

    processes = np.min([processes,os.cpu_count()])

    return mp.Pool(processes=processes, initializer=mute, *args, **kwargs)