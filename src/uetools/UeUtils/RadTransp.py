class RadTransp():


    def profile_function(self, x1, x2, x3, y1, y2, y3, y4, y5, plot=False):
        from numpy import zeros, arange, interp
        from scipy.interpolate import splrep, BSpline
        
        self.set_psinormc()
        psi0 = self.psinormc[0]
        psif = self.psinormc[-1]
        points = zeros((6,2))
        points[0] = [psi0, y1]
        points[1] = [psi0 + x1*(1-psi0), y2]
        points[2] = [points[1,0] + x2*(1-points[1,0]), y2]
        points[3] = [1, y3]
        points[4] = [1 + x3*(psif - 1), y4]
        points[5] = [psif, y5]
        if self.bayopt['interpolator'] == 'spline':
#            print(' USING SPLINE INTERPOLATOR')
            tck = splrep(points[:,0], points[:,1])
            xnew = arange(psi0, psif, (psif-psi0)/250) 
            
            if plot:
                f=self.plotprofile(points[:,0], points[:,1], marker='o', linestyle='')
                f.get_axes()[0].plot(xnew, BSpline(*tck)(xnew), '-')
                f.get_axes()[0].plot(self.psinormc, BSpline(*tck)(self.psinormc), 'k-')
            return BSpline(*tck)(self.psinormc)
        elif self.bayopt['interpolator'] == 'linear':
#            print(" USING LINEAR INTERPOLATOR")
            return interp(self.psinormc, points[:,0], points[:,1])

    def set_bayesian_profiles(self, x):
        dif_use = self.getue('dif_use', cpy=False) 
        kye_use = self.getue('kye_use', cpy=False) 
        kyi_use = self.getue('kyi_use', cpy=False) 
    
#        print('x', x)
#        print('kx', x[:8])
#        print('dx', x[8:])
        ky = self.profile_function(*x[:8])
        dif = self.profile_function(*x[8:])
        ky[ky<=0] = 1e-2
        dif[dif<=0] = 1e-2
#        print('ky', ky)
#        print('dif', dif)

        for i in range(self.ixpt1+1,self.ixpt2+1):
            for j in range(dif_use.shape[-1]):
                dif_use[i,:,j] =  dif
                kye_use[i] = ky
        kyi_use = kye_use
        return dif_use, kye_use, kyi_use
#        self.setue('dif_use', dif_use)
#        self.setue('kye_use', kye_use)
#        self.setue('kyi_use', kyi_use)


    def eval_fit(self, x, vals):
        from numpy import sum, mean
        try:
            fit = self.profile_function(*x)
        except:
            return 1e5
        ss_res = sum((fit - vals)**2)
        ss_tot = sum((vals - mean(vals))**2)
        return abs(ss_res/ss_tot)

    def fit_coefs(self, x0, vals):
        from scipy.optimize import minimize

        return minimize(lambda x : self.eval_fit(x, vals), x0,
                method='Nelder-Mead', options={'disp': True, 
                'fatol': 1e-4}).x

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
            'kyi_use': {'target': kyi},
        }
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
        radval = 1000*exp(radfraction*20 - 13)*self.bayopt['optimize_radtransp']
        radval = nan_to_num(radval)
        print(" BLACKBOX EVALUATION COMPLETED")
        print(f'    Radiated fraction: {radfraction}')
        print(f'    Radiated fraction impact: {radval}')
        print('    Least squares, n: {}'.format(abs(rr_n)))
        print('    Least squares, T: {}'.format(abs(rr_t)))
        print('    Punishment for incompleteness: {}'.format( 1000*(1-self.lastsuccess)*(success is not True)))
        val =  abs(rr_n) + abs(rr_t) - 2 + 1000*(1-self.lastsuccess)*(success is not True) + radval

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


    def optimize_transport(self, ne=None, te=None, random_state=None,
        n_calls=10, n_initial_points=5, initres=1e3, savedir='test', savename='test{}', oldstate = None, interpolator='spline',
        optimize_radtransp=True):
        from skopt.space import Real
        from skopt import gp_minimize, load
        from skopt.callbacks import CheckpointSaver
        from numpy import interp
        from scipy.interpolate import splrep, BSpline
        from numpy import arange
        self.iter = 0
        if te is None:
            x = [0.95, 0.98, 1, 1.01, 1.03, 1.06, 1.2]
            xfine = arange(x[0], x[-1], (x[-1]-x[0])/250)
            y = [250, 230, 100, 30, 10, 3, 1]
#            te = BSpline(*splrep(x, y))(xfine)
            te = interp(xfine, x, y)
        if ne is None:
            x = [0.95, 0.98, 1, 1.01, 1.03, 1.06, 1.2]
            xfine = arange(x[0], x[-1], (x[-1]-x[0])/250)
            y = [2.2e19, 1.9e19, 1.5e19, 1e19, .8e19, .5e19, .3e18]
#            ne = BSpline(*splrep(x, y))(xfine)
            ne = interp(xfine, x, y)            
            
        self.x = []
        self.y = []
        self.target_ne = ne
        self.target_nex = xfine
        self.target_te = te
        self.target_tex = xfine

        self.bayopt = { 'savedir': savedir,
                        'savename': savename,
                        'initres': initres,
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

        kfit = self.fit_coefs([0.5, 0.1, 0.8, 2, 0.5, 1, 2, 2], self.get('kye_use')[self.get('ixmp')])
        dfit = self.fit_coefs([0.5, 0.1, 0.8, 2, 0.5, 1, 2, 2], self.get('dif_use')[self.get('ixmp'),:,0])
        dfit[3:] *=1.005
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
        
 
        bounds = [
            Real(kfit[0]*0.5, min(0.99, kfit[0]*1.5), name='kx1'),
            Real(kfit[1]*0.5, min(0.99, kfit[1]*1.5), name='kx2'),
            Real(kfit[2]*0.5, min(0.99, kfit[2]*1.5), name='kx3'),
            Real(kfit[3]*0.1, kfit[3]*3, name='kx4'),
            Real(kfit[4]*0.1, kfit[4]*3, name='kx5'),
            Real(kfit[5]*0.1, kfit[5]*3, name='kx6'),
            Real(kfit[6]*0.1, kfit[6]*3, name='kx7'),
            Real(kfit[7]*0.1, kfit[7]*3, name='kx8'),
            Real(dfit[0]*0.5, min(0.99, dfit[0]*1.5), name='dx1'),
            Real(dfit[1]*0.5, min(0.99, dfit[1]*1.5), name='dx2'),
            Real(dfit[2]*0.5, min(0.99, dfit[2]*1.5), name='dx3'),
            Real(dfit[3]*0.1, dfit[3]*3, name='dx4'),
            Real(dfit[4]*0.1, dfit[4]*3, name='dx5'),
            Real(dfit[5]*0.1, dfit[5]*3, name='dx6'),
            Real(dfit[6]*0.1, dfit[6]*3, name='dx7'),
            Real(dfit[7]*0.1, dfit[7]*3, name='dx8'),
        ]

  


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


    
    def iterate_radtransp(self, savefname, psite, teexp, psine, neexp, 
        sibdrys=None, simagxs=None, iters=10, steps=4, **kwargs):
        from numpy import linspace
        i = 0
        if steps != 0:
            for frac in linspace(1/steps, 1, steps):
                self.step_radtransp(frac, psite, teexp, psine, neexp, 
                    sibdrys=sibdrys, simagxs=simagxs, **kwargs)
                self.converge(savefname='{}_radtransp_frac{}_iter{}'.format(\
                    savefname, str(frac).replace('.', ''), i+1), dtreal=1e-10)
                if self.get('iterm') != 1:
                    print('Convergence failed!')
                    return
                i +=1
        for itr in range(i, iters):
            self.step_radtransp(1, psite, teexp, psine, neexp, 
                sibdrys=sibdrys, simagxs=simagxs, **kwargs)
            self.converge(savefname='{}_radtransp_iter{}'.format(savefname, i+1), 
                dtreal=1e-10)
            if self.get('iterm') != 1:
                print('Convergence failed!')
                return
            i +=1


    def step_radtransp(self, frac, psite, teexp, psine, neexp, sibdrys=None, 
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

