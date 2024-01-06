class PostProcessors:
    


    def xvert(self, arr):
        from numpy import zeros_like
        """" Returns values of arr on x-vertices (staggered grid)

        index      x-1      x       x+1  
            _________________________
        iy  |  ix-1 |   ix  | ix+1  |
            |_______|_______|_______|
        """
        try:
            self.ixp1
        except:
            self.ixp1 = self.get('ixp1')
        retvar = zeros_like(arr)
        retvar[:-1] = 0.5*(arr[:-1] + arr[1:])
        if self.snull:
            # Account for the X-point cuts
            for ix in [self.ixpt1, self.ixpt2]:
                for iy in range(self.iysptrx+1):
                    retvar[ix, iy] = 0.5*(
                        arr[ix, iy] +  arr[self.ixp1[ix, iy], iy])
        else:
            raise ValueError('dnull xvert not implemented')
        return retvar

    def dxv(self, arr):
        from numpy import zeros
        """ Returns differentiaties arr at the vertices
        """
        return
        try:
            self.ixm1
        except:
            self.ixm1 = self.get('ixm1')
        retvar = zeros_like(arr)
        retvar[:-1] = (arr[1:] + arr[:-1])
        if self.snull:
            # Account for the X-point cuts
            for ix in [self.ixpt1, self.ixpt2]:
                for iy in range(self.iysptrx+1):
                    retvar[ix, iy] = 0.5*(
                        arr[ix, iy] +  arr[self.ixp1[ix, iy], iy])
        else:
            raise ValueError('dnull xvert not implemented')
        return retvar


        return

    def dxc(self, arr):
        from numpy import zeros
        """ Returns differentiaties arr at cell centers
        """
        return


    def yvert(self, arr):
        if self.snull:
            1
        else:
            raise ValueError('dnull yvert not implemented')
        return
    

    def set_psinormc(self, simagxs=None, sibdrys=None):
        if simagxs is None:
            simagxs = self.get('simagxs')
        if sibdrys is None:
            sibdrys = self.get('sibdrys')
        self.psinormc =  (self.get('psi')[self.get('ixmp'),:,0]-simagxs)/ \
            (sibdrys-simagxs)

    def set_psinormf(self, simagxs=None, sibdrys=None):
        if simagxs is None:
            simagxs = self.get('simagxs')
        if sibdrys is None:
            sibdrys = self.get('sibdrys')
        psi = self.get('psi')[self.get('ixmp')]
        self.psinormf = ((0.5*(psi[:,3] + psi[:,4])) -simagxs)/ \
            (sibdrys-simagxs)


    def calc_forcebalance(self, cutlo=1e-300):
        """ Calculates the force-balance terms
        """
        from numpy import minimum, zeros, zeros_like, sum
        from sys import modules
        self.forcebalance = {}

        # Load all variables necessary to the local namespace
        for var in ['ev', 'misotope', 'natomic', 'mi', 'zi', 'rrv', 'gpex', 
            'gpix', 'gtex', 'gtix', 'pri', 'netap', 'qe', 'fqp', 'sx',
            'alfe', 'loglambda', 'betai', 'ex', 'vol', 'up', 'volmsor',
            'pondomfpari_use',
        ]:
            setattr(modules[__name__], var, self.get(var))       

        # Set up necessary arrays
        den = zeros((misotope, self.get('nchstate')) + self.get('ne').shape )
        gradt = zeros((misotope, self.get('nchstate')) + self.get('ne').shape )
        gradp = zeros((misotope, self.get('nchstate')) + self.get('ne').shape )
        upi = zeros(self.get('ne').shape + (sum(natomic),))
        upi_gradp = zeros(self.get('ne').shape + (sum(natomic),))
        upi_alfe = zeros(self.get('ne').shape + (sum(natomic),))
        upi_betai = zeros(self.get('ne').shape + (sum(natomic),))
        upi_ex = zeros(self.get('ne').shape + (sum(natomic),))
        upi_volmsor = zeros(self.get('ne').shape + (sum(natomic),))
        taudeff = zeros(self.get('ne').shape + (sum(natomic),))

        # Get necessary data
        ni = self.get('ni')
        den[0, 0] = self.xvert(self.get('ne'))
        den[1, 0] = self.xvert(ni[:, :, 0])
        tempa = self.xvert(self.get('te'))
        tif = self.xvert(self.get('ti'))

        # Get flux-limit factor
        ltmax = minimum(    abs(tempa/(rrv*gtex + cutlo)),
                            abs(tif/(rrv*gtix + cutlo)),
                            abs(den[0, 0]*tempa/(rrv*gpex + cutlo))
        )
        # Approx Coulomb mfps
        lmfpe = 1e16*( (tempa/ev)**2/(den[0, 0] + cutlo))
        lmfpi = 1e16*( (tif/ev)**2/(den[0, 0] + cutlo))
        # Get ion-pressure gradient scale-lengths
        ltmax = minimum(
                    ltmax, 
                    abs(self.xvert(pri[:,:,0])/(rrv*gpix[:,:,0] + cutlo))
        )
        # Get flux-limit factor
        flxlimf = 1/ (1 + self.get('fricflf') * ((lmfpe + lmfpi)/\
                (ltmax + cutlo))**2)
        gradp[0, 0] = rrv * gpex
        gradt[0, 0] = rrv * gtex
        # Get total impurity density
        # TODO: What happens to misotope and natomic when ishymol>0??
        ionspecies = [natomic[x] for x in range(2, misotope)]
#        dztot = sum(
#            self.xvert(self.get('ni'))[2:2+sum(ionspecies)], 
#            axis=2
#        )
        ifld = 0
        # Loop over the different impurity species
        for misa in range(2, misotope):
            # Loop over each impurity species charge state
            for nz in range(natomic[misa]):
                ifld += 1
                # Skip neutral atoms
                ifld += (self.get('ziin')[ifld] < 1e-10)
                # Calculate the hydrogen-impurity scattering rate
                den[misa, nz] = self.xvert(ni[:, :, ifld])
                zeffv = self.xvert(self.get('zeff'))
                gradt[misa, nz] = rrv*gtix
                gradp[misa, nz] = rrv*gpix[:,:,ifld] \
                    - pondomfpari_use[:, :, ifld]
                if self.get('is_z0_imp_const') == 0:
                    z0 = den[0, 0] * zeffv / (den[1, 0] +cutlo) - 1
                else:
                    z0 = z0_imp_const
                taud = self.get('cftaud')*5.624e54*mi[0]**0.5*mi[ifld]*\
                    tif**1.5 / (loglambda*den[misa, nz]*zi[ifld]**2*
                    (mi[0]+mi[ifld]) + cutlo)
                taudeff[:,:,ifld] = flxlimf*taud*den[misa, nz]*(1+2.65*z0)*\
                    (1+0.285*z0)/(den[0, 0] * (1+0.24*z0) * \
                    (1+0.93*z0) + cutlo)

                # Store the components of the force-balance equation
                upi_gradp[:,:,ifld] = -gradp[misa, nz]/ (den[misa,nz] + cutlo)
                upi_gradp[:,:,ifld] *= (taudeff[:,:,ifld]/mi[0])
                upi_alfe[:,:,ifld] = alfe[ifld]*gradt[0, 0] 
                upi_alfe[:,:,ifld] *= (taudeff[:,:,ifld]/mi[0])
                upi_betai[:,:,ifld] = betai[ifld]*gradt[misa, nz]
                upi_betai[:,:,ifld] *= (taudeff[:,:,ifld]/mi[0])
                upi_ex[:,:,ifld] = qe*zi[ifld]*rrv*ex
                upi_ex[:,:,ifld] *= (taudeff[:,:,ifld]/mi[0])
                upi_volmsor[:,:,ifld] = volmsor[:, :, ifld] / \
                    (den[misa, nz] * vol + cutlo)
                upi_volmsor[:,:,ifld] *= (taudeff[:,:,ifld]/mi[0])
                # Solve force-balance equation for impurity velocity
                upi[:, :, ifld] = up[:, :, 0] + upi_gradp[:,:,ifld] + \
                    upi_alfe[:,:,ifld] + upi_betai[:,:,ifld] + \
                    upi_ex[:,:,ifld] + upi_volmsor[:,:,ifld]

        self.forcebalance['upi'] = upi
        self.forcebalance['up'] = up[:,:,0]
        self.forcebalance['upi_gradp'] = upi_gradp
        self.forcebalance['upi_alfe'] = upi_alfe
        self.forcebalance['upi_betai'] = upi_betai
        self.forcebalance['upi_ex'] = upi_ex
        self.forcebalance['taudeff'] = taudeff
        self.forcebalance['upi_volmsor'] = upi_volmsor
        


        return upi, self.forcebalance
                 
    def calc_momentum_balance(self):
        """ Calculates the momentum and force-balance terms based on rdmom_res_up_v8.py
        originally written by Tom Rognlin, Menglong Zhao and Filippo Scotti

        A script that computes components of partial_up/partial_t,
        Units of terms = [m/s**2] = (vpar/time) i.e. acceleration.
        
        8 components correspond to array names:
        (fmoxv,floyv,pridxv,efieldfv,thermfv,momfricv,momizv,momcxrecv)
        Index associated with the component names above: (1,2,3,4,5,6,7,8).

        Right-hand side of total momentum eqn is up_net(ix,iy,ifld).
        Amplitude order of up_net components are in ordcmpname

        """
        from numpy import zeros,abs,empty,int
        from sys import modules
        self.momentum ={"momentum":{},"force":{},"amplitude":{},"rank":{}}

        # Load all variables necessary to the local namespace
        for var in ['ev', 'misotope', 'natomic', 'mi', 'zi', 'rrv', 'gpex', 
            'gpix', 'gtex', 'gtix', 'pri', 'netap', 'qe', 'fqp', 'sx',
            'alfe', 'loglambda', 'betai', 'ex', 'vol', 'up', 'volmsor',
            'pondomfpari_use','nx','ny','nusp','ixm1','ixp1','fmix','fmixy',
            'gxf','volv','ni','uu','fmiy','cpiup','qe','frici','fricnrl',
            'cfmsor','msor','msorxr','psor','psorxr','iysptrx','dx','rr',
            'nuiz','iigsp','cfneut','nm','cfupcx','nurc','nucx','fnix',
            'fniy',
        ]:
            setattr(modules[__name__], var, self.get(var))       

        # Set up necessary arrays
        # Components of momentum eqn
        fmox = zeros((nx+2,ny+2,nusp))
        fmoxconv = zeros((nx+2,ny+2,nusp))
        fmoxvis = zeros((nx+2,ny+2,nusp))
        fmoy = zeros((nx+2,ny+2,nusp))
        pridx = zeros((nx+2,ny+2,nusp))
        efieldf = zeros((nx+2,ny+2,nusp))
        thermf = zeros((nx+2,ny+2,nusp))
        momfric = zeros((nx+2,ny+2,nusp))
        momiz = zeros((nx+2,ny+2,nusp))
        momcxrec = zeros((nx+2,ny+2,nusp))
        sxv = zeros((nx+2,ny+2))

        # Components of upar eqn
        fmoxv = zeros((nx+2,ny+2,nusp))
        fmoyv = zeros((nx+2,ny+2,nusp))
        pridxv = zeros((nx+2,ny+2,nusp))
        efieldfv = zeros((nx+2,ny+2,nusp))
        thermfv = zeros((nx+2,ny+2,nusp))
        momfricv = zeros((nx+2,ny+2,nusp))
        momizv = zeros((nx+2,ny+2,nusp))
        momcxrecv = zeros((nx+2,ny+2,nusp))

       

        # Sum of components
        up_net = zeros((nx+2,ny+2,nusp))
        mom_net = zeros((nx+2,ny+2,nusp))

        xparrb = zeros((nx+2,ny+2))    # parallel distance from right bdry

        # Compute components of momentum equation; these are resmo/volv
        for iis in range(nusp):  # loop over all momentum eqns
            for iy in range(1,ny+1):
                for ix in range(1,nx):
                    ix1 = ixm1[ix,iy]
                    ix2 = ixp1[ix,iy]
                    vratio = sx[ix,iy]/(gxf[ix,iy]*volv[ix,iy])
            
                    fmox[ix,iy,iis] = ( -(fmix[ix2,iy,iis] - fmixy[ix2,iy,iis]) + \
                                (fmix[ix, iy,iis] - fmixy[ix, iy,iis]) )/volv[ix,iy]
                    sxv[ix,iy] = 0.5*(sx[ix1,iy]+sx[ix,iy])
                    sxv[ix2,iy] = 0.5*(sx[ix,iy]+sx[ix2,iy])
                    fmoxconv[ix,iy,iis] =mi[iis]*(-ni[ix2,iy,iis]*uu[ix2,iy,iis]*up[ix2,iy,iis]* \
                                                                        sxv[ix2,iy] \
                                +ni[ix,iy,iis]*uu[ix,iy,iis]*up[ix,iy,iis]*sxv[ix,iy])/ \
                                                                            volv[ix,iy]
                    fmoxvis[ix,iy,iis] = fmox[ix,iy,iis] - fmoxconv[ix,iy,iis]

                    fmoy[ix,iy,iis] = ( -fmiy[ix,iy,iis] + fmiy[ix,iy-1,iis])/volv[ix,iy]

                    pridx[ix,iy,iis] = -cpiup[iis]*gpix[ix,iy,iis]*rrv[ix,iy]*vratio

                    efieldf[ix,iy,iis] =  qe*zi[iis]*0.5*(ni[ix2,iy,iis]+ni[ix,iy,iis])* \
                                                ex[ix,iy]*rrv[ix,iy]*vratio

                    thermf[ix,iy,iis] = frici[ix,iy,iis]*vratio

                    momfric[ix,iy,iis] = fricnrl[ix,iy,iis]

                    momiz[ix,iy,iis] = cfmsor*msor[ix,iy,iis]/volv[ix,iy]

                    momcxrec[ix,iy,iis] = cfmsor*msorxr[ix,iy,iis]/volv[ix,iy]

        # For hydrogen ion (iis=0) and atoms (iis=1), charge-exchange terms given below:
        #iigsp - 1  ... python indexing
        for iy in range(1,ny+1):
            for ix in range(1,nx):
                ix2 = ixp1[ix,iy]
                momiz[ix,iy,0] = cfneut*0.25*(nuiz[ix,iy,0]+nuiz[ix2,iy,0])* \
                                        (nm[ix,iy,iigsp-1]+nm[ix2,iy,iigsp-1])* \
                                                            up[ix,iy,iigsp-1]
                momcxrec[ix,iy,0] = cfneut*cfupcx*0.25* \
                                        (nucx[ix,iy,0]+nucx[ix2,iy,0])* \
                                        (nm[ix,iy,iigsp-1]+nm[ix2,iy,iigsp-1])* \
                                            (up[ix,iy,iigsp-1]-up[ix,iy,0]) - \
                                cfneut*0.25*(nurc[ix,iy,0]+nurc[ix2,iy,0])* \
                                            (nm[ix,iy,0]+nm[ix2,iy,0])* \
                                                up[ix,iy,0]
                momiz[ix,iy,1] = -momiz[ix,iy,0]
                momcxrec[ix,iy,1] = - momcxrec[ix,iy,0]

        # Add all terms from the momentum density equation
        mom_net = fmox+fmoy+pridx+efieldf+thermf+momfric+momiz+momcxrec

        self.momentum['momentum']['fmox'] = fmox
        self.momentum['momentum']['fmoy'] = fmoy
        self.momentum['momentum']['pridx'] = pridx
        self.momentum['momentum']['efieldf'] = efieldf
        self.momentum['momentum']['thermf'] = thermf
        self.momentum['momentum']['momfric'] = momfric
        self.momentum['momentum']['momiz'] = momiz
        self.momentum['momentum']['momcxrec'] = momcxrec
        self.momentum['momentum']['mom_net'] = mom_net

        ###############################################################################
        ################ next compute  components for partial(v_par)/partial_t ########
        # Compute components of vpar equation                                  ########
        # Note: because of staggered up mesh wrt ni mesh, when adding the contribution
        # from mi*up*dn/dt, must average over 2 poloidal cells; should still conserve
        # up because each 2 cells in simple 50/50 average is in steady-state.

        for iis in range(nusp):  # loop over all momentum eqns
            for iy in range(1,ny+1):
                for ix in range(1,nx):
                    ix1 = ixm1[ix,iy]
                    ix2 = ixp1[ix,iy]
                    nbar = 0.5*(ni[ix,iy,iis]+ni[ix2,iy,iis])
            
                    # Explicitly showning fnix(ix,iy,ifld) even though it cancels
                    fmoxv[ix,iy,iis] = ( fmox[ix,iy,iis]/mi[iis] \
                                - up[ix,iy,iis]*0.5*(fnix[ix1,iy,iis]-fnix[ix ,iy,iis] \
                                                    +fnix[ix ,iy,iis]-fnix[ix2,iy,iis]) )\
                                                                            /nbar

                    fmoyv[ix,iy,iis] = ( fmoy[ix,iy,iis]/mi[iis] \
                                -up[ix,iy,iis]*0.5*(-fniy[ix ,iy,iis]+fniy[ix ,iy-1,iis] \
                                                    -fniy[ix2,iy,iis]+fniy[ix2,iy-1,iis]) )\
                                                                            /nbar

                    pridxv[ix,iy,iis] = pridx[ix,iy,iis]/(mi[iis]*nbar)

                    efieldfv[ix,iy,iis] =  efieldf[ix,iy,iis]/(mi[iis]*nbar)

                    thermfv[ix,iy,iis] = thermf[ix,iy,iis]/(mi[iis]*nbar)

                    momfricv[ix,iy,iis] = momfric[ix,iy,iis]/(mi[iis]*nbar)

                    momizv[ix,iy,iis] = momiz[ix,iy,iis]/(mi[iis]*nbar) \
                                -up[ix,iy,iis]*0.5*(psor[ix,iy,iis]+psor[ix2,iy,iis])/nbar

                    momcxrecv[ix,iy,iis] = momcxrec[ix,iy,iis]/(mi[iis]*nbar) \
                                    -up[ix,iy,iis]*0.5*(psorxr[ix,iy,iis]+psorxr[ix2,iy,iis])\
                                                                            /nbar

        # For hydrogen ion (is=1) and atoms (is=2), charge-exchange terms given below:
        for iis in range(2):
            for iy in range(1,ny+1):
                for ix in range(1,nx):
                    ix2 = ixp1[ix,iy]
                    nbar = 0.5*(ni[ix,iy,iis]+ni[ix2,iy,iis])
                    momizv[ix,iy,iis] = momiz[ix,iy,iis]/(mi[iis]*nbar) \
                                -up[ix,iy,iis]*0.5*(psor[ix,iy,iis]+psor[ix2,iy,iis])/nbar

                    momcxrecv[ix,iy,iis] = momcxrec[ix,iy,iis]/(mi[iis]*nbar) \
                            -up[ix,iy,iis]*0.5*(psorxr[ix,iy,iis]+psorxr[ix2,iy,iis])/nbar

        # Add all terms from the momentum energy density equation
        up_net = fmoxv+fmoyv+pridxv+efieldfv+thermfv+momfricv+momizv+momcxrecv

        self.momentum['force']['fmoxv'] = fmoxv
        self.momentum['force']['fmoyv'] = fmoyv
        self.momentum['force']['pridxv'] = pridxv
        self.momentum['force']['efieldfv'] = efieldfv
        self.momentum['force']['thermfv'] = thermfv
        self.momentum['force']['momfricv'] = momfricv
        self.momentum['force']['momizv'] = momizv
        self.momentum['force']['momcxrecv'] = momcxrecv
        self.momentum['force']['up_net'] = up_net

        # Compute absolute values
        afmoxv = abs(fmoxv)
        afmoyv = abs(fmoyv)
        apridxv = abs(pridxv)
        aefieldfv = abs(efieldfv)
        athermfv = abs(thermfv)
        amomfricv = abs(momfricv)
        amomizv = abs(momizv)
        amomcxrecv = abs(momcxrecv)

        self.momentum['amplitude']['afmoxv'] = afmoxv
        self.momentum['amplitude']['afmoyv'] = afmoyv
        self.momentum['amplitude']['apridxv'] = apridxv
        self.momentum['amplitude']['aefieldfv'] = aefieldfv
        self.momentum['amplitude']['athermfv'] = athermfv
        self.momentum['amplitude']['amomfricv'] = amomfricv
        self.momentum['amplitude']['amomizv'] = amomizv
        self.momentum['amplitude']['amomcxrecv'] = amomcxrecv


        cmpname = empty(8,dtype=('U',9)) #component name
        maxcomp = zeros((nx+2,ny+2,nusp))  # Maxim component amplitude code
        ordcmpname = empty((nx+2,ny+2,nusp,8),dtype=('U',9)) #cmp order by name
        ordcmpidx = zeros((nx+2,ny+2,nusp,8),dtype=int)      #cmp order by index

        #this can be simplified...
        for iis in range(nusp):   #identify max component
            for iy in range(1,ny+1):
                for ix in range(1,nx):
                    ordcmpidx[ix,iy,iis,0] = 1
                    ordcmpname[ix,iy,iis,0] = "fmoxv"
                    mxc = afmoxv[ix,iy,iis]
                    if (afmoyv[ix,iy,iis] > mxc):
                        ordcmpidx[ix,iy,iis,0] = 2
                        ordcmpname[ix,iy,iis,0] = "fmoyv"
                        mxc = afmoyv[ix,iy,iis]
                    if (apridxv[ix,iy,iis] > mxc): 
                        ordcmpidx[ix,iy,iis,0] = 3
                        ordcmpname[ix,iy,iis,0] = "pridxv"
                        mxc = apridxv[ix,iy,iis]
                    if (aefieldfv[ix,iy,iis] > mxc):
                        ordcmpidx[ix,iy,iis,0] = 4
                        ordcmpname[ix,iy,iis,0] = "efieldfv"
                        mxc = aefieldfv[ix,iy,iis]
                    if (athermfv[ix,iy,iis] > mxc):
                        ordcmpidx[ix,iy,iis,0] = 5
                        ordcmpname[ix,iy,iis,0] = "thermfv"
                        mxc = athermfv[ix,iy,iis]
                    if (amomfricv[ix,iy,iis] > mxc):
                        ordcmpidx[ix,iy,iis,0] = 6
                        ordcmpname[ix,iy,iis,0] = "momfricv"
                        mxc = amomfricv[ix,iy,iis]
                    if (amomizv[ix,iy,iis] > mxc):
                        ordcmpidx[ix,iy,iis,0] = 7
                        ordcmpname[ix,iy,iis,0] = "momizv"
                        mxc = amomizv[ix,iy,iis]
                    if (amomcxrecv[ix,iy,iis] > mxc):
                        ordcmpidx[ix,iy,iis,0] = 8
                        ordcmpname[ix,iy,iis,0] = "momcxrecv"
                        mxc = amomcxrecv[ix,iy,iis]
                    maxcomp[ix,iy,iis] = mxc

        for iis in range(nusp):   #identify second largest component
            for iy in range(1,ny+1):
                for ix in range(1,nx):
                    maxval = 0.
                    if (ordcmpidx[ix,iy,iis,0] != 1):
                        maxval = afmoxv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,1] = 1
                        ordcmpname[ix,iy,iis,1] = "fmoxv"
                    if (ordcmpidx[ix,iy,iis,0] != 2 and maxval<afmoyv[ix,iy,iis]):
                        maxval = afmoyv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,1] = 2
                        ordcmpname[ix,iy,iis,1] = "fmoyv"
                    if (ordcmpidx[ix,iy,iis,0] != 3 and maxval<apridxv[ix,iy,iis]):
                        maxval = apridxv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,1] = 3
                        ordcmpname[ix,iy,iis,1] = "pridxv"
                    if (ordcmpidx[ix,iy,iis,0] != 4 and maxval<aefieldfv[ix,iy,iis]):
                        maxval = aefieldfv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,1] = 4
                        ordcmpname[ix,iy,iis,1] = "efieldfv"
                    if (ordcmpidx[ix,iy,iis,0] != 5 and maxval<athermfv[ix,iy,iis]):
                        maxval = athermfv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,1] = 5
                        ordcmpname[ix,iy,iis,1] = "thermfv"
                    if (ordcmpidx[ix,iy,iis,0] != 6 and maxval<amomfricv[ix,iy,iis]):
                        maxval = amomfricv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,1] = 6
                        ordcmpname[ix,iy,iis,1] = "momfricv"
                    if (ordcmpidx[ix,iy,iis,0] != 7 and maxval<amomizv[ix,iy,iis]):
                        maxval = amomizv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,1] = 7
                        ordcmpname[ix,iy,iis,1] = "momizv"
                    if (ordcmpidx[ix,iy,iis,0] != 8 and maxval<amomcxrecv[ix,iy,iis]):
                        maxval = amomcxrecv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,1] = 8
                        ordcmpname[ix,iy,iis,1] = "momcxrecv"

        for iis in range(nusp):   #identify third largest component
            for iy in range(1,ny+1):
                for ix in range(1,nx):
                    maxval = 0.
                    if (ordcmpidx[ix,iy,iis,0] != 1 and ordcmpidx[ix,iy,iis,1] != 1):
                        maxval = afmoxv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,2] = 1
                        ordcmpname[ix,iy,iis,2] = "fmoxv"
                    if (ordcmpidx[ix,iy,iis,0] != 2 and ordcmpidx[ix,iy,iis,1] != 2 and maxval<afmoyv[ix,iy,iis]):
                        maxval = afmoyv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,2] = 2
                        ordcmpname[ix,iy,iis,2] = "fmoyv"
                    if (ordcmpidx[ix,iy,iis,0] != 3 and ordcmpidx[ix,iy,iis,1] != 3 and maxval<apridxv[ix,iy,iis]):
                        maxval = apridxv[ix,iy,iis] 
                        ordcmpidx[ix,iy,iis,2] = 3
                        ordcmpname[ix,iy,iis,2] = "pridxv"
                    if (ordcmpidx[ix,iy,iis,0] != 4 and ordcmpidx[ix,iy,iis,1] != 4 and maxval<aefieldfv[ix,iy,iis]):
                        maxval = aefieldfv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,2] = 4
                        ordcmpname[ix,iy,iis,2] = "efieldfv"
                    if (ordcmpidx[ix,iy,iis,0] != 5 and ordcmpidx[ix,iy,iis,1] != 5 and maxval<athermfv[ix,iy,iis]):
                        maxval = athermfv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,2] = 5
                        ordcmpname[ix,iy,iis,2] = "thermfv"
                    if (ordcmpidx[ix,iy,iis,0] != 6 and ordcmpidx[ix,iy,iis,1] != 6 and maxval<amomfricv[ix,iy,iis]):
                        maxval = amomfricv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,2] = 6
                        ordcmpname[ix,iy,iis,2] = "momfricv"
                    if (ordcmpidx[ix,iy,iis,0] != 7 and ordcmpidx[ix,iy,iis,1] != 7 and maxval<amomizv[ix,iy,iis]):
                        maxval = amomizv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,2] = 7
                        ordcmpname[ix,iy,iis,2] = "momizv"
                    if (ordcmpidx[ix,iy,iis,1] != 8 and ordcmpidx[ix,iy,iis,2] != 8 and maxval<amomcxrecv[ix,iy,iis]):
                        maxval = amomcxrecv[ix,iy,iis]
                        ordcmpidx[ix,iy,iis,2] = 8
                        ordcmpname[ix,iy,iis,2] = "momcxrecv"

        self.momentum['rank']['ordcmpidx'] = ordcmpidx
        self.momentum['rank']['ordcmpname'] = ordcmpname
        self.momentum['rank']['maxcomp'] = maxcomp

    def pradpltwl(self):
        from numpy import zeros, pi, cos
        from math import atan2

        ny = self.get('ny')
        nx = self.get('nx')
        nxpt = self.get('nxpt')
        ixlb = self.get('ixlb')
        ixrb = self.get('ixrb')
        rm = self.get('rm')
        zm = self.get('zm')
        sx = self.get('sx')
        angfx = self.get('angfx')
        vol = self.get('vol')
        eeli = self.get('eeli')
        ebind = self.get('ebind')
        ev = self.get('ev')
        psor = self.get('psor')
        erlrc = self.get('erlrc')
        isimpon = self.get('isimpon')
        # Initialize arrays
        self.pwr_pltz = zeros((ny+2, 2*nxpt)) 
        self.pwr_plth = zeros((ny+2, 2*nxpt)) 
        self.pwr_wallz = zeros((nx+2))
        self.pwr_wallh = zeros((nx+2))
        self.pwr_pfwallz = zeros(((nx+2)*nxpt)) 
        self.pwr_pfwallh = zeros(((nx+2)*nxpt)) 
        if (isimpon > 0):
            prdu = self.get('prad')
        else:
            prdu = zeros((nx+2, ny+2))
        nj = self.get('nxomit')
        for ip in range(0, 2*nxpt):
            ixv = (1-(ip % 2))*ixlb[0] + (ip % 2)*(ixrb[0]+1)
            print(ixv)
            for iyv in range(1, ny+1):
                for iy in range(1, ny+1):
                    for ix in range(1, nx+1):
                        theta_ray1 = atan2( 
                                zm[ixv+nj, iyv, 1]-zm[ix+nj, iy, 0],
                                rm[ixv+nj, iyv, 1]-rm[ix+nj, iy, 0]
                            )
                        theta_ray2 = atan2(
                                zm[ixv+nj, iyv, 3]-zm[ix+nj, iy, 0],
                                rm[ixv+nj, iyv, 3]-rm[ix+nj, iy, 0]
                            )
                        if ((ix==10) and (iy==10) and (iyv==10)):
                            print(theta_ray1, theta_ray2)
                        dthgy = abs(theta_ray1-theta_ray2)
                        frth = min(dthgy, 2*pi-dthgy)/(2*pi)
                        sxo = sx[ixv, iyv]/cos(angfx[ixv, iyv])
                        if ((ix==10) and (iy==10) and (iyv==10)):
                            print(dthgy, frth, sxo)
                        self.pwr_pltz[iyv, ip] += prdu[ix,iy]*vol[ix,iy]*frth/sxo
                        self.pwr_plth[iyv, ip] += ((eeli[ix,iy] - ebind*ev)*psor[ix,iy,0]+erlrc[ix, iy]) * frth/sxo
        self.pwr_pltz[0, ip] = self.pwr_pltz[1, ip]
        self.pwr_pltz[-1, ip] = self.pwr_pltz[-2, ip]
        self.pwr_plth[0, ip] = self.pwr_plth[1, ip]
        self.pwr_plth[-1, ip] = self.pwr_plth[-2, ip]
