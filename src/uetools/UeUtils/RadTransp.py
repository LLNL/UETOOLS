class RadTransp():
    def __init__(self, case):
        self.setuserdiff = case.input.userdiff
        pass


    def profile_function(self, *args, plot=False):
        from numpy import zeros, arange, interp, linspace
        from scipy.interpolate import splrep, BSpline
        
        self.set_psinormc()
        res = len(args)
        x = linspace(self.psinormc[0], self.psinormc[-1], res)
        y = args

        if self.bayopt['interpolator'] == 'spline':
#            print(' USING SPLINE INTERPOLATOR')
            tck = splrep(x, y)
            xnew = arange(psi0, psif, (psif-psi0)/250) 
            
            if plot:
                f=self.plotprofile(points[:,0], points[:,1], marker='o', linestyle='')
                f.get_axes()[0].plot(xnew, BSpline(*tck)(xnew), '-')
                f.get_axes()[0].plot(self.psinormc, BSpline(*tck)(self.psinormc), 'k-')
            return BSpline(*tck)(self.psinormc)
        elif self.bayopt['interpolator'] == 'linear':
#            print(" USING LINEAR INTERPOLATOR")
            return interp(self.psinormc, x, y)

    def set_bayesian_profiles(self, x):
        res = int(len(x)/2)
        
        dif_use = self.getue('dif_use', cpy=False) 
        kye_use = self.getue('kye_use', cpy=False) 
        kyi_use = self.getue('kyi_use', cpy=False) 
    
#        print('x', x)
#        print('kx', x[:8])
#        print('dx', x[8:])
        ky = self.profile_function(*x[:res])
        dif = self.profile_function(*x[res:])
        ky[ky<=0] = 1e-2
        dif[dif<=0] = 1e-2
#        print('ky', ky)
#        print('dif', dif)

        for i in range(self.ixpt1+1,self.ixpt2+1):
            for j in range(dif_use.shape[-1]):
                dif_use[i,:,j] =  dif
                kye_use[i] = ky
        if self.bayopt['setkyi']:
            kyi_use = kye_use
        return dif_use, kye_use, kyi_use
#        self.setue('dif_use', dif_use)
#        self.setue('kye_use', kye_use)
#        self.setue('kyi_use', kyi_use)


    def blackbox(self, *args, savedir='test', savename='test{}', initres=1e4):
        from numpy import interp, sum, mean, array, exp, nan_to_num
        from glob import glob
        from os import path

        try:
            print('Accessing Bayopt')
            self.__getattribute__('bayopt')
            print('Bayopt present - set values')
            savedir = self.bayopt['savedir']
            savename = self.bayopt['savename']
            initres = self.bayopt['initres']
        except:
            pass


        [dif, kye, kyi] = self.set_bayesian_profiles(*args)
        self.iter += 1
        var = {
            'dif_use': {'target': dif},
            'kye_use': {'target': kye},
        }
        if self.bayopt['setkyi']:
            var['kyi_use'] = {'target': kyi},
        print(" STARTING NEW BLACKBOX ITERATION = {}".format(self.iter))
        # TODO: identify the best match and always restore that
        print(" RESTORING SAVE FILE")
        if len(self.x) == 0:
            self.restore_save(self.savefile)
        else:
            i = array(self.y).argmin() + 1
            bestdir = '/'.join([savedir, savename.format(i)])
            # TODO: search for best match to current params instead?
            lastfile = max(glob('/'.join([bestdir,'*.hdf5'])), key=path.getctime)
            self.restore_save(lastfile)
        print(" SETTING DIFFUSIVITIES")
        if len(self.x) == 0:
            self.setuserdiff(self.diff_file)
        else:
            self.setuserdiff(lastfile)
        print(" STARTING CONTINUATION SOLVER")
        self.setue('iterm', 0)
        print(" CHECKING FNRM")
        self.populate()
        try:
            success = self.continuation_solve(var, savedir = '/'.join(
                    [savedir, savename.format(self.iter)]), 
                    initres=initres,# itermx=7
            )
        except Exception as error:
            print("Continuation solve failed:")
            print(f"    {error}")
            return 10
            
        
        
        ss_res_t = sum((interp(self.psinormc, self.target_tex, self.target_te) \
                - self.get('te')[self.get('ixmp')]/self.get('ev'))**2)
        ss_tot_t = sum(self.target_te - mean(self.target_te)**2)
        ss_t = ss_res_t/ss_tot_t
        rr_t = 1-ss_t

        ss_res_n = sum((interp(self.psinormc, self.target_nex, self.target_ne) \
                - self.get('ne')[self.get('ixmp')])**2)
        ss_tot_n = sum(self.target_ne - mean(self.target_ne)**2)
        ss_n = ss_res_n/ss_tot_n
        rr_n = 1-ss_n
        
        # TODO: Store x0, y0, x, y to save
        # TODO: Punish exponentially for >10% radial power transport
        fey = (self.get('feey')+self.get('feiy'))[self.get('ixpt1')[0]+1:self.get('ixpt2')[0]+1]
        radfraction = sum(fey[:,-2])/sum(fey[:,0])
        radval = 1000*exp(radfraction*20 - 13)*self.bayopt['optimize']
        radval = nan_to_num(radval)
        print(" BLACKBOX EVALUATION COMPLETED")
        print(f'    Radiated fraction: {radfraction}')
        print(f'    Radiated fraction impact: {radval}')
        print('    Least squares, n: {}'.format(abs(rr_n)))
        print('    Least squares, T: {}'.format(abs(rr_t)))
        print('    Punishment for incompleteness: {}'.format( 1000*(1-self.lastsuccess)*(success is not True)))
        val =  10*(abs(rr_n) + abs(rr_t) - 2) + 1000*(1-self.lastsuccess)*(success is not True) + radval

        val = min(val, 1e3)
        print('==========================')
        print(f'Total function value: {val}')
        self.x.append(args)
        self.y.append(val)
        return val
        # TODO: Converge using continuation solve
        # Restore successful save and diffusion coeffs
        # Penalize if continuation solve does not succeed

#        self.converge(savename='/'.join([savedir, savename]))



        # TODO: define optimization function


    def optimize(self, ne=None, te=None, psine=None, psite=None, random_state=None,
        n_calls=10, n_initial_points=5, initres=1e3, savedir='test', savename='test{}', 
        oldstate = None,  resolution = 8, optimize_radtransp=True, interpolator='linear',
        setkyi=False):
        from skopt.space import Real
        from skopt import gp_minimize, load
        from skopt.callbacks import CheckpointSaver
        from numpy import interp
        from scipy.interpolate import splrep, BSpline
        from numpy import arange, ones, linspace
        self.iter = 0
        if te is None:
            x = [0.95, 0.98, 1, 1.01, 1.03, 1.06, 1.2]
            xfine = arange(x[0], x[-1], (x[-1]-x[0])/250)
            y = [250, 230, 100, 30, 10, 3, 1]
#            te = BSpline(*splrep(x, y))(xfine)
            te = interp(xfine, x, y)
            psite = xfine
        if ne is None:
            x = [0.95, 0.98, 1, 1.01, 1.03, 1.06, 1.2]
            xfine = arange(x[0], x[-1], (x[-1]-x[0])/250)
            y = [2.2e19, 1.9e19, 1.5e19, 1e19, .8e19, .5e19, .3e18]
#            ne = BSpline(*splrep(x, y))(xfine)
            ne = interp(xfine, x, y)
            psine = xfine
            
        self.set_psinormc()
        self.x = []
        self.y = []
        self.target_ne = ne
        self.target_nex = psine
        self.target_te = te
        self.target_tex = psite

        self.bayopt = { 'savedir': savedir,
                        'savename': savename,
                        'initres': initres,
                        'setkyi': setkyi,
                        'interpolator': interpolator,
                        'optimize_radtransp': optimize_radtransp,
        }
    
        bounds = [
            Real(0.02, .98, name='kx1'),
            Real(0.02, .98, name='kx2'),
            Real(0.02, .98, name='kx3'),
            Real(1e-2, 10, name='ky1'),
            Real(1e-2, 10, name='ky2'),
            Real(1e-2, 10, name='ky3'),
            Real(1e-2, 10, name='ky4'),
            Real(1e-2, 10, name='ky5'),
            Real(0.02, .98, name='dx1'),
            Real(0.02, .98, name='dx2'),
            Real(0.02, .98, name='dx3'),
            Real(1e-2, 10, name='dy1'),
            Real(1e-2, 10, name='dy2'),
            Real(1e-2, 10, name='dy3'),
            Real(1e-2, 10, name='dy4'),
            Real(1e-2, 10, name='dy5'),
        ]

        kfit = interp(  linspace(self.psinormc[0], self.psinormc[-1], resolution), 
                        self.psinormc, self.get('kye_use')[self.get('ixmp')]
        )
        dfit = interp(  linspace(self.psinormc[0], self.psinormc[-1], resolution), 
                        self.psinormc, self.get('dif_use')[self.get('ixmp'),:,0]
        )


        x0 = list(kfit) + list(dfit)
        y0 = None

        if oldstate is not None:
            res = load(oldstate)
            x0 = res.x_iters
            y0 = res.func_vals

        # TODO: 
        # Fit parametrized functions to current values
        # Call black box with best fits to get initial guess
        # Seed x0 y0 using best fits

        # TODO: Add check for manual abort
        
        bounds = []
        for val in kfit:
#            bounds.append(Real(val*0.1, val*3))
            bounds.append(Real(1e-2, 10))
        for val in dfit:
#            bounds.append(Real(val*0.1, val*3))
            bounds.append(Real(1e-2, 10))


        res = gp_minimize(  self.blackbox,
                            bounds,
                            n_calls = n_calls,
                            n_initial_points = n_initial_points,
                            noise = 1e-10,
                            callback = [CheckpointSaver(f"{savedir}/checkpoint.pkl", compress=9)],
                            x0 = x0,
                            y0 = y0,
                            random_state = random_state,
                            initial_point_generator='lhs')


# self.blackbox([0.4,0.2,0.5,1,.1,.3,2,2.5,0.5,0.2,0.5,1,.1,.3,2,2.5])
        '''
        self.set_psinormc()
        self.plot(xfine, te, True, color='r')
        self.plot(self.psinormc, self.get('te')[self.get('ixmp')]/1.602e-19, color='k')
        self.plot(xfine, ne, True, color='r')
        self.plot(self.psinormc, self.get('ne')[self.get('ixmp')], color='k')
        '''


    
    def iterate(self, savefname, psite, teexp, psine, neexp, 
        sibdrys=None, simagxs=None, iters=10, steps=4, **kwargs):
        from numpy import linspace
        i = 0
        if steps != 0:
            for frac in linspace(1/steps, 1, steps):
                self.step(frac, psite, teexp, psine, neexp, 
                    sibdrys=sibdrys, simagxs=simagxs, **kwargs)
                self.converge(savefname='{}_radtransp_frac{}_iter{}'.format(\
                    savefname, str(frac).replace('.', ''), i+1), dtreal=1e-10)
                if self.get('iterm') != 1:
                    print('Convergence failed!')
                    return
                i +=1
        for itr in range(i, iters):
            self.step(1, psite, teexp, psine, neexp, 
                sibdrys=sibdrys, simagxs=simagxs, **kwargs)
            self.converge(savefname='{}_radtransp_iter{}'.format(savefname, i+1), 
                dtreal=1e-10)
            if self.get('iterm') != 1:
                print('Convergence failed!')
                return
            i +=1


    def step(self, frac, psite, teexp, psine, neexp, sibdrys=None, 
        simagxs=None, **kwargs):
        from uedge import bbb
        diff = self.diffrad(psine, neexp, sibdrys=sibdrys, simagxs=simagxs, **kwargs)
        kye = self.kyerad(psite, teexp, psine, neexp, sibdrys=sibdrys, 
            simagxs=simagxs, **kwargs)
        #kyi = self.kyirad(psite, teexp, psine, neexp, sibdrys=sibdrys, 
        #    simagxs=simagxs, **kwargs)

        for ix in range(self.get('ixpt1')[0]+1, self.get('ixpt2')[0]+1):
            for isp in range(bbb.dif_use.shape[-1]):
                bbb.dif_use[ix,:,isp] = frac*diff + (1-frac)*bbb.dif_use[ix,:,isp]
            bbb.kye_use[ix] = frac*kye + (1-frac)*bbb.kye_use[ix]
#            bbb.kyi_use[ix] = frac*kyi + (1-frac)*bbb.kyi_use[ix]
        bbb.dif_use[:,:,1] = 0
        bbb.kyi_use = bbb.kye_use

    def diffrad(self, psiexp, neexp, sibdrys=None, simagxs=None,
        lbound=5e-2, ubound=200, **kwargs):
        from scipy.interpolate import interp1d
        from numpy import gradient
        expfunc = interp1d(psiexp, neexp)
        # Assert required psi values are calculated
        try:
            self.psinormc
        except:
            self.set_psinormc(sibdrys, simagxs)
        try:
            self.psinormf
        except:
            self.set_psinormf(sibdrys, simagxs)
        ixmp = self.get('ixmp')
        gamma = (self.get('fniy')[:,:,0]/self.get('sy'))[ixmp,:]
        uefunc = interp1d(self.psinormf, gamma, kind='linear',
            fill_value='extrapolate')
        diffnew = -uefunc(self.psinormc)/gradient(\
#            expfunc(self.psinormc), self.psinormc)
            expfunc(self.psinormc), self.get('yyc'))
        # Fix boundaries to not use GC values
        diffnew[0] = diffnew[1]
        diffnew[-1] = diffnew[-2]
        diffnew[diffnew<lbound] = lbound
        diffnew[diffnew>ubound] = ubound

        from matplotlib.pyplot import subplots
        f, ax = subplots()
        ax.plot(self.psinormc, self.get('dif_use')[ixmp,:,0], 'k-')
        ax.plot(self.psinormc, diffnew, 'r-')
        return diffnew

    def kyerad(self, psite, teexp, psine, neexp, sibdrys=None, simagxs=None,
        lbound=5e-2, ubound=200, **kwargs):
        from scipy.interpolate import interp1d
        from numpy import gradient
        expte = interp1d(psite, teexp*1.602e-19)
        expne = interp1d(psine, neexp)
        # Assert required psi values are calculated
        try:
            self.psinormc
        except:
            self.set_psinormc(sibdrys, simagxs)
        try:
            self.psinormf
        except:
            self.set_psinormf(sibdrys, simagxs)
        ixmp = self.get('ixmp')
        # Turn off convective gas flow?
#        cfloyeold = self.get('cfloye')
#        self.setue('cfloye', 0)
#        self.populate(verbose=False)
        # Store conductive flux only
        gamma = (self.get('fniy')[:,:,0]/self.get('sy'))[ixmp,:]
        # Restore defaults
#        self.setue('cfloye', cfloyeold)
#        self.populate(verbose=False)
        qerad = (self.get('feey')/self.get('sy'))[ixmp,:]
        gammafunc = interp1d(self.psinormf, gamma, kind='linear',
            fill_value='extrapolate')
        qefunc = interp1d(self.psinormf, qerad, kind='linear',
            fill_value='extrapolate')
        kyenew = -(qefunc(self.psinormc) - \
            (5/2)*gammafunc(self.psinormc)*self.get('te')[ixmp]) /\
            (expne(self.psinormc)*gradient(expte(self.psinormc),
            self.get('yyc')))
#            (expne(self.psinormc)*gradient(expte(self.psinormc),self.psinormc))
        kyenew[0] = kyenew[1]
        kyenew[-1] = kyenew[-2]
        kyenew[kyenew<lbound] = lbound
        kyenew[kyenew>ubound] = ubound
        
        from matplotlib.pyplot import subplots
        f, ax = subplots()
        ax.plot(self.psinormc, self.get('kye_use')[ixmp,:], 'k-')
        ax.plot(self.psinormc, kyenew, 'r-')
        return kyenew


    def kyirad(self, psiti, tiexp, psine, neexp, sibdrys=None, simagxs=None,
        lbound=5e-2, ubound=200, **kwargs):
        from scipy.interpolate import interp1d
        from numpy import gradient
        expti = interp1d(psiti, tiexp*1.602e-19)
        expne = interp1d(psine, neexp)
        # Assert required psi values are calculated
        try:
            self.psinormc
        except:
            self.set_psinormc(sibdrys, simagxs)
        try:
            self.psinormf
        except:
            self.set_psinormf(sibdrys, simagxs)
        ixmp = self.get('ixmp')
        # Turn off convective gas flow?
#        cfloyiold = self.get('cfloyi')
#        self.setue('cfloyi', 0)
#        self.populate(verbose=False)
        # Store conductive flux only
        gamma = (self.get('fniy')[:,:,0]/self.get('sy'))[ixmp,:]
        # Restore defaults
#        self.setue('cfloyi', cfloyiold)
        self.populate(verbose=False)
        qirad = (self.get('feiy')/self.get('sy'))[ixmp,:]
        gammafunc = interp1d(self.psinormf, gamma, kind='linear',
            fill_value='extrapolate')
        qifunc = interp1d(self.psinormf, qirad, kind='linear',
            fill_value='extrapolate')
        kyinew = -(qifunc(self.psinormc) - \
            (5/2)*gammafunc(self.psinormc)*self.get('ti')[ixmp]) /\
            (expne(self.psinormc)*gradient(expti(self.psinormc),
            self.get('yyc')))
#            (expne(self.psinormc)*gradient(expti(self.psinormc),self.psinormc))
        kyinew[0] = kyinew[1]
        kyinew[-1] = kyinew[-2]
        kyinew[kyinew<lbound] = lbound
        kyinew[kyinew>ubound] = ubound
        
        from matplotlib.pyplot import subplots
        f, ax = subplots()
        ax.plot(self.psinormc, self.get('kyi_use')[ixmp,:], 'k-')
        ax.plot(self.psinormc, kyinew, 'r-')
        return kyinew

