

class RadTransp():
    
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

