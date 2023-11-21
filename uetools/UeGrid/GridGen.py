from matplotlib.pyplot import ion
ion()
''' PLOT GRID VISUALIZATION '''


class GridGen():
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def plot_flx(self, first=None, last=None, ax=None):
        ''' Plots flux surfaces as defined in grid setup
            Based on plotflx.bas by Rensink/Rognlien/Porter 
        '''
        from uedge import com, grd, flx
        from matplotlib.pyplot import subplots
        flx.flxrun()
        if ax is None:
            f, ax = subplots(figsize=(7,9))
        if first is None:
            first = 0
        if last is None:
            last = 2*(com.nycore[0]+com.nysol[0]+com.nyout[0])+2
        first = max(0, first)
        # Plot vessel if present in EQDSK files
        if com.nlim > 0:
            ax.plot(com.xlim, com.ylim, 'k-', linewidth=2)
        # Plot target plates, if specified
        try:
            ax.plot(grd.rplate1, grd.zplate1, 'ro-')
        except:
            pass
        try:
            ax.plot(grd.rplate2, grd.zplate2, 'b.-')
        except:
            pass
        # Plot the flux surfaces within the sepcified range
        for i in range(first,last+1):
            # Plot SOL flux surfaces for each half-mesh
            if ((i >= com.jmin[0]-1) and (i < com.jsptrx[0])) or \
                ((i >= com.jsptrx[1]) and (i <= com.jmax[1])):
                ax.plot(com.xcurve[:,i][abs(com.xcurve[:,i])>0], 
                        com.ycurve[:,i][abs(com.ycurve[:,i])>0], 'k-', linewidth=0.3)
            # Plot CORE/PFR flux surfaces for each half-mesh
            elif ((i >= com.jsptrx[0]) and (i <= com.jmax[0])) or \
                    ((i >= com.jmin[1]-1) and (i <= com.jsptrx[1]+1)):
                ax.plot(com.xcurve[:,i][abs(com.xcurve[:,i])>0][:flx.ijumpf[i]], 
                        com.ycurve[:,i][abs(com.ycurve[:,i])>0][:flx.ijumpf[i]], 
                        'k-', linewidth=0.3)
                ax.plot(com.xcurve[:,i][abs(com.xcurve[:,i])>0][flx.ijumpf[i]:], 
                        com.ycurve[:,i][abs(com.ycurve[:,i])>0][flx.ijumpf[i]:], 
                        'k-', linewidth=0.3)
        ax.set_aspect('equal')
        ax.set_xlabel('Horizontal position [m]')
        ax.set_ylabel('Vertical position [m]')
        return ax.get_figure()

    def plot_poloidal_distribution(self, yaxis='log', ylim = [1e-5, 5]):
        ''' Plots poloidal settings of model 
            Plot before executing grdrun?
        '''
        from uedge import com, flx, grd
        from matplotlib.pyplot import figure
        from numpy import linspace
        
        f = figure(figsize=(14,9))
        gs = f.add_gridspec(2,2, height_ratios=[6,1], wspace=0, hspace=0.3, top=0.98, bottom=0.05)
        axul = f.add_subplot(gs[0,0])
        axur_buf = f.add_subplot(gs[0,1])
        axur = axur_buf.twinx()
        axb = f.add_subplot(gs[1,:])
        x = linspace(0,1,200)
        if yaxis == 'lin':
            pli = axul.plot
            plo = axur.plot
        elif yaxis == 'log':
            pli = axul.semilogy
            plo = axur.semilogy
        else:
            print('yaxis option not recognized, using linear')
            pli = axur.plot
            plo = axur.plot
        flx.flxrun()
        grd.grdrun()
        nxtot = sum(com.nxleg[0])+sum(grd.nxuse)
        if grd.kxmesh == 0:
            print('Manual mesh seeding interface TBI')
            return None
        elif grd.kxmesh == 1:
            func = grd.xfcn
        elif grd.kxmesh == 2:
            func = grd.xfcn2
        elif grd.kxmesh == 3:
            func = grd.xfcn3
        elif grd.kxmesh == 4:
            func = grd.xfcn4
        else:
            print('grd.kxmesh option not recognized, aborting')
            return None
        if grd.kxmesh != 4:
            pli(x, [func(a) for a in x], 'k-')
            pli([a/nxtot for a in range(nxtot)], [func(a/nxtot) for a in range(nxtot)], 'r.')
            plo(x, [func(1)-func(a) for a in x], 'k-')
            plo([a/nxtot for a in range(nxtot)], [func(1)-func(a/nxtot) for a in range(nxtot)], 'r.')
            plot_distro(func, ax=axb)
        else:
            pli(x, [func(a, nxtot) for a in x], 'k-')
            pli([a/nxtot for a in range(nxtot)], [func(a/nxtot, nxtot) for a in range(nxtot)], 'r.')
            plo(x, [func(1, nxtot)-func(a, nxtot) for a in x], 'k-')
            plo([a/nxtot for a in range(nxtot)], [func(1, nxtot)-func(a/nxtot, nxtot) for a in range(nxtot)], 'r.')
            plot_distro(func, nxtot, ax=axb)
       
        for ax in [axul, axur]:
            ax.axvline(com.nxleg[0,0]/nxtot, color='r', linewidth=3)
            ax.axvline(1-com.nxleg[0,1]/nxtot, color='r', linewidth=3)
            ax.set_ylim(ylim)
        
        axul.set_xlabel('Normalized distance between targets, index space')
        axur_buf.set_xlabel('Normalized distance between targets, index space')
        axur_buf.set_yticks([])
        axul.set_ylabel('Distance from inner target [m]')
        axur.set_ylabel('Distance from outer target [m]')
        axul.set_xlim(0,0.5)
        axur.set_xlim(0.5-1e-6,1)
        return ax.get_figure()

    def plot_distro(self, func, *args, ax=None, fname='',
                outpath='/Users/holm10/Documents/fusion/analysis/gridcomp_22/figs'):
        from uedge import com, flx, grd
        from matplotlib.pyplot import subplots

        if ax is None:
            f, ax = subplots(1, 1, figsize=(10,3))
        ax2 = ax.twinx()    
        nxtot = sum(com.nxleg[0])+sum(grd.nxuse)
        for ix in range(nxtot):
            x = func(ix/nxtot,*args)
            ax.axvline(func(ix/nxtot,*args), color='k', linewidth=1)
        ax.axvline(func(com.nxleg[0,0]/nxtot, *args), color='r', linewidth=3)
        ax.axvline(func(1-com.nxleg[0,1]/nxtot, *args), color='r', linewidth=3)
        ax.set_xlabel('Poloidal distance along separatrix [m]')
        ax.set_ylabel('Inner target')
        ax2.set_ylabel('Outer target')    
        ax.set_yticks([])
        ax2.set_yticks([])
        ax.set_title('Poloidal cell distribution along flux-tube')
        ax.set_ylim(0,1)
        ax.set_xlim(0,func(1,*args))
        return ax.get_figure()

    def plot_efit(aeqdskfname='aeqdsk', geqdskfname='neqdsk', ax=None, ncontour=80, color='grey', linestyle='solid', sepcolor='k', linewidth=0.5):
        from matplotlib.pyplot import subplots
        from copy import deepcopy
        from numpy import linspace
        from uedge import com, flx
        from os.path import exists
        from scipy.interpolate import interp2d

        # Backup original pointers
        oldaeqdskfname = deepcopy(com.aeqdskfname)
        oldgeqdskfname = deepcopy(com.geqdskfname)
        # Set new file paths
        com.aeqdskfname = aeqdskfname
        com.geqdskfname = geqdskfname


        # Check whether the aeqdsk file can be located: if not, do not execute aeqdsk()
        if exists(aeqdskfname):
            flx.aeqdsk()

        if exists(geqdskfname):
            flx.neqdsk()
        else:
            print('EFIT geqdsk file "{}" not found.\nAborting...'.format(\
                geqdskfname))
            return

        if ax is None:
            f, ax = subplots(figsize=(7,9))

        # Reconstruct EFIT grid
        x = linspace(0, com.xdim, com.nxefit)+com.rgrid1
        y = linspace(0, com.zdim, com.nyefit)-(com.zdim*0.5-com.zmid)

        # Interpolate on EFIT grid to find X-points
        interp = interp2d(x, y, com.fold.transpose())

        ax.contour(x, y, com.fold.transpose(), ncontour, colors=color, 
            linewidths=linewidth, linestyles=linestyle) 

        # Check whether the upper X-point exists
        if (com.rseps2 >= x.min()) and (com.rseps2 <= x.max()):
            if (com.zseps2 >= y.min()) and (com.zseps2 <= y.max()):
                upperxpoint = interp(com.rseps2, com.zseps2)
                ax.contour(x, y, com.fold.transpose(), [upperxpoint], 
                    colors=sepcolor, linewidths=1, linestyles='solid') 
        
        # Check whether the lower X-point exists
        if (com.rseps >= x.min()) and (com.rseps <= x.max()):
            if (com.zseps >= y.min()) and (com.zseps <= y.max()):
                lowerxpoint = interp(com.rseps, com.zseps)
                ax.contour(x, y, com.fold.transpose(), [lowerxpoint], 
                    colors=sepcolor, linewidths=1, linestyles='solid') 
        try:
            ax.plot(com.xlim, com.ylim, 'k-', linewidth=2)
        except:
            pass
        ax.set_aspect('equal')
        ax.set_xlabel('Horizontal position [m]')
        ax.set_ylabel('Vertical position [m]')
        ax.set_title(com.runid[0].decode('UTF-8').strip())

        # Restore original pointers
        com.aeqdskfname = oldaeqdskfname
        com.geqdskfname = oldgeqdskfname

        return ax.get_figure()


    def find_efitdata(geqdsk, ncontour=250, **kwargs):
        from uedge import com
        from scipy.interpolate import interp2d
        from scipy.optimize import fmin
        from numpy import linspace, array, gradient, sum
        from matplotlib.pyplot import subplots, ginput

        f, ax = subplots(figsize=(7,9))
        ax.plot([],[], 'ko')
        ax.plot([],[], 'ro')
        ax.plot([],[], 'bo')
        plot_efit('', geqdsk, ax, ncontour=ncontour, **kwargs)

        print('Manually identify the following points in order')
        print(' - Magnetic axis')
        print(' - Lower X-point')
        print(' - Upper X-point')
        print('Choose the points by clicking on them in the figure')
        print('Undo by right-clicking')

        textbox = f.get_axes()[0].text(.87, 1.5, 'Manually identify the '
            'following points in order\n'
            ' - Magnetic axis\n - Lower X-point\n - Upper X-point\nChoose'
            ' the points by clicking on them in the figure\nUndo by '
            'right-clicking',backgroundcolor='w', zorder=10)
        magx, lxpt, uxpt = ginput(3,0)
        textbox.set_visible(False)
        
        # Reconstruct EFIT grid
        x = linspace(0, com.xdim, com.nxefit)+com.rgrid1
        y = linspace(0, com.zdim, com.nyefit)-(com.zdim*0.5-com.zmid)
        
        gradient = sum(abs(array(gradient(com.fold))), axis=0)
        flux = interp2d(x, y, com.fold.transpose())
        gradflux = interp2d(x, y, gradient.transpose())
        
        def evalfunc(coords, func):
            return func(*coords)

    #    ax1.contour(x, y, gradient.transpose(), 200, colors='r', 
    #        linewidths=.5, linestyles='solid') 

        magx = fmin(evalfunc, array(magx), args=(flux,), disp=False)
        lxpt = fmin(evalfunc, array(lxpt), args=(gradflux,), disp=False)
        uxpt = fmin(evalfunc, array(uxpt), args=(gradflux,), disp=False)
        ax.plot(*magx, 'ko')
        ax.plot(*uxpt, 'ro')
        ax.plot(*lxpt, 'bo')
        com.simagx = flux(*magx)[0]
        com.sibdry1 = flux(*lxpt)[0]
        com.sibdry2 = flux(*uxpt)[0]
        ax.legend(['Magnetic axis = {:.4}'.format(com.simagx),
            'Lower X-point = {:.4}'.format(com.sibdry1),
            'Upper X-point = {:.4}'.format(com.sibdry2)])

    #    print(*flux(*magx))
    #    print(*flux(*uxpt))
    #    print(*flux(*lxpt))
    #    print(*flux(1.173,.9465))


    def gendnull(aeqdskfile, neqdskfile,psi0min1=0.970, psi0max_outer=1.04, 
        psi0max_inner=1.03, psi0min2_lower=1.01, psi0min2_upper=0.993, 
        slp2fac=0.8, slp3fac=0.8, nycore=16, nysol=14, nyout=14,
        nxlegu=[9,11], nxcoreu=[8,10], slpxtu=1.2, nxxptu=1, nxmodu=4, alfxptuu=0.5, 
        alfxptlu=1, isrefxptnu=1, nsmoothxu=3, ismmonu=3,nsmoothu=3, wtmesh1u=0, 
        dmix0u=0.5,
        nxlegl=[9,11], nxcorel=[8,10], slpxtl=1.2, nxxptl=1, nxmodl=4, alfxptul=0.5, 
        alfxptll=1, isrefxptnl=1, nsmoothxl=3, ismmonl=3,nsmoothl=3, wtmesh1l=0, 
        dmix0l=0.5
        ):
        from gridgen import plot_efit, plot_flx
        from Forthon import gallot, gchange
        from numpy import array, pi
        from copy import deepcopy
    #    from uedge.gridue import write_gridue
        from uedge import com, bbb, flx, grd
        # This script implements the procedure for generating a full up/down
        # asymmetric double-null mesh.

        # The following files should be present when executing this script:
        #	rdmesh.dnull.top
        #	rdmesh.dnull.bot
        #	plate.d3d_top
        #	plate.d3d_bot

        # The following steps generate a mesh that includes both x-points:
        #	1) read the eqdsk to identify x-points and separatrices
        #	2) specify the radial distribution of flux surfaces
        #	3) construct the mesh for the upper half of the configuration
        #	4) construct the mesh for the lower half of the configuration
        #	5) combine the two halves of the mesh
        #	6) define guard cells at the divertor target plates
        #	7) compute magnetic field values for each cell
        #	8) write out the gridue file

        # NOTE: the guard cells at the radial boundaries are not defined by
        # the above procedure; that is done by the bbb package as part of the
        # exmain command.

        # NOTE: must generate bottom-half mesh last in order to get 
        #	normal-orientation eqdsk data in call to subroutine
        #	magnetics before writing the final 'gridue' file.

        ####################################################################
        ### STEP 1: read the eqdsk to identify x-points and separatrices ###
        ####################################################################


        # Set paths to EFIT files
        com.aeqdskfname = aeqdskfile
        com.geqdskfname = neqdskfile
        
        
        # mesh construction --
        bbb.mhdgeo=1
        # Read eqdsks to identify x-points
        com.geometry="dnull"
        flx.readefit()

        # Input factor to ensure that psi0sep is slightly outside of LCFS
        flx.psifac=1.0001

        # Set normalized separatrix flux values
        flx.psi0sep1 = flx.psifac*(com.sibdry1-com.simagx)/(com.sibdry-com.simagx)
        flx.psi0sep2 = flx.psifac*(com.sibdry2-com.simagx)/(com.sibdry-com.simagx)

        if flx.psi0sep1 < flx.psi0sep2:
            print('\n*** Lower divertor primary, dRsep < 0 ***')
        else:
            print('\n *** Upper divertor primary, dRsep > 0 ***')


        # NOTE: sibdry is the un-normalized flux at the last closed flux surface (LCFS)

        ## Diagnostic plots of the two separatrices and limiter surface:
    #    plot_efit(aeqdskfile, neqdskfile)

        ################################################################
        ### STEP 2: specify the radial distribution of flux surfaces ###
        ################################################################

        # Input flux boundary values
        flx.psi0min1 = psi0min1		# must be less than 1.00
        flx.psi0max_outer = psi0max_outer
            # must be greater than max(psi0sep1,psi0sep2)
        flx.psi0max_inner = psi0max_inner	# must be greater than max(psi0sep1,psi0sep2)
        flx.psi0min2_lower = psi0min2_lower	# must be less than psi0sep1
        flx.psi0min2_upper = psi0min2_upper	# must be less than psi0sep2
        flx.psi0max = flx.psi0max_outer	# don't change this line

        if psi0max_outer < max(flx.psi0sep1,flx.psi0sep2):
            print("\n*** INPUT ERROR: psi0max_outer less than '\
                'max(psi0sep1,psi0sep2)\n")
            return
        elif psi0max_inner < max(flx.psi0sep1,flx.psi0sep2):
            print("\n*** INPUT ERROR: psi0max_inner less than '\
                'max(psi0sep1,psi0sep2)\n")
            return
        elif psi0min2_lower >= flx.psi0sep1:
            print("\n*** INPUT ERROR: psi0min2_lower greater than psi0sep1\n")
            return
        if psi0min2_upper >= flx.psi0sep2:
            remark("\n*** INPUT ERROR: psi0min2_upper greater than psi0sep2\n")
            return

        # Input number of cells in each radial region
        com.nycore[0] = nycore	# closed flux surfaces
        com.nysol[0] = nysol	# open to one divertor only; =0 for balanced dnull
        com.nyout[0] = nyout	# open to both divertors

        # For balanced (but up/down asymmetric) double null configuration
        # (OK if psi0sep1,2 are not too different and radial mesh is not too fine):
        if com.nysol[0] == 0:
            flx.psi0sep1 = max(flx.psi0sep1, flx.psi0sep2)
            flx.psi0sep2 = flx.psi0sep1

        # Shape factors for radial mesh
        flx.slp2fac = slp2fac	# =1 for uniform mesh between separatrices
                    # <1 concentrates mesh near primary separatrix
        flx.slp3fac = slp3fac# =1 for uniform mesh between separatrices
                    # <1 concentrates mesh near secondary separatrix

        ### NOTE: for shot 125852 at 3325 ms, treated as balanced dnull (nysol=0)
        ### a) separatrix looks a little "kinked" near upper x-point
        ### b) radial resolution is limited to nycore <= 10 for slp2fac=slp3fac=1.0
        ### c) radial resolution is limited to slp2fac,slp3fac > 0.75 for nycore=10

        # Option for radial mesh construction
        flx.kymesh=0	# user-specified flux contour values for radial mesh

        # Set dimensions for flx package arrays
        com.nym = com.nycore[0] + com.nysol[0] + com.nyout[0]
        com.jdim = 2*com.nym + 3
        gallot('Flxin', 0)
        gallot('Workdn', 0)

    #        call allot("flx.psitop",0)
    #        call allot("flx.psibot",0)
    #        call gallot("flx.Workdn", 0)

        # Compute flux contour values:
        t = array(range(0, com.nym+1), dtype='float64')
        r2p = array(0.)
        r3p = array(0.)
        jj1 = 0
        t1 = jj1
        jj2 = jj1 + com.nycore[0]
        t2 = jj2
        jj3 = jj2 + com.nysol[0]
        t3 = jj3
        jj4 = jj3 + com.nyout[0]
        t4 = jj4


        # First, psi0_mp_outer:
        r1 = flx.psi0min1
        r2 = min(flx.psi0sep1, flx.psi0sep2)
        r3 = max(flx.psi0sep1, flx.psi0sep2)
        r4 = flx.psi0max_outer

        if t3 > t2:	# use 3-region model for unbalanced dnull
            flx.rho3dn(t, flx.psi0_mp_outer, com.nym+1, t1, t2, t3, t4, r1, r2, r3, 
                r4, flx.slp2fac, flx.slp3fac, r2p, r3p)
        else:			# use 2-region model for balanced dnull
           r2p=slp2fac*(r2-r1)/(t2-t1)
           flx.rho5(t, flx.psi0_mp_outer, com.nym+1, t1, t2, t4, r1, r2, r4, r2p)
           r3p=r2p


        # Next, psi0_mp_inner:
        flx.psi0_mp_inner[:jj3+1] = flx.psi0_mp_outer[:jj3+1]
        r4 = flx.psi0max_inner
        flx.rho1l(t[jj3:], flx.psi0_mp_inner[jj3:], jj4-jj3, t3, t4, 
            r3, r4, r3p)

        # Next, psi0_dp_upper_outer:
        if flx.psi0sep1 < flx.psi0sep2:	# lower divertor is primary
           jjs=jj3
           rs=r3
           rsp=deepcopy(r3p)
           ts=t3
        else:					# upper divertor is primary
           jjs=jj2
           rs=r2
           rsp=deepcopy(r2p)
           ts=t2

        flx.psi0_dp_upper_outer[jjs:jj4+1] = flx.psi0_mp_outer[jjs:jj4+1]
        r1= flx.psi0min2_upper
        flx.rho1r(t, flx.psi0_dp_upper_outer, jjs-jj1, t1, ts, r1, rs, rsp)

        # Next, psi0_dp_upper_inner:
        flx.psi0_dp_upper_inner[jj1:jjs+1] = flx.psi0_dp_upper_outer[jj1:jjs+1]
        flx.psi0_dp_upper_inner[jjs:jj4+1] = flx.psi0_mp_inner[jjs:jj4+1]

        # Next, psi0_dp_lower_outer:
        if flx.psi0sep1 < flx.psi0sep2:	# lower divertor is primary
           jjs=jj2
           rs=r2
           rsp=deepcopy(r2p)
           ts=t2
        else:					# upper divertor is primary
           jjs=jj3
           rs=r3
           rsp=deepcopy(r3p)
           ts=t3

        flx.psi0_dp_lower_outer[jjs:jj4+1] = flx.psi0_mp_outer[jjs:jj4+1]
        r1 = flx.psi0min2_lower
        flx.rho1r(t, flx.psi0_dp_lower_outer, jjs-jj1, t1, ts, r1, rs, rsp)

        # Next, psi0_dp_lower_inner:
        flx.psi0_dp_lower_inner[jj1:jjs+1] = flx.psi0_dp_lower_outer[jj1:jjs+1]
        flx.psi0_dp_lower_inner[jjs:jj4+1] = flx.psi0_mp_inner[jjs:jj4+1]

        ## Diagnostic plot of radial distribution of flux surfaces at outboard midplane
        ## plot ylim xlim style=dashed scale=equal labels=blank
        ## plotz (fold-simagx)/(sibdry-simagx) xold yold lev=psi0_mp_outer
        ## frame 2.26 2.36 zmagx-.05 zmagx+.05

        # flux contour search parameters:
        flx.istchkon=1
        flx.dtheta_exclude    = array([.75,.50])*pi
        flx.dtheta_overlap_pf = array([.05,.05])*pi
        flx.altsearch=1

        #echo=oldecho

        
        ##########################################################################
        ### STEP 3: construct the mesh for the upper half of the configuration ###
        ##########################################################################

        halfmesh(upper=True, nxleg=nxlegu, nxcore=nxcoreu, slpxt=slpxtu, nxxpt=nxxptu, 
            nxmod=nxmodu, alfxptu=alfxptuu, alfxptl=alfxptlu, isrefxptn=isrefxptnu, 
            nsmoothx=nsmoothxu, ismmon=ismmonu, nsmooth=nsmoothu, wtmesh1=wtmesh1u, 
            dmix0=dmix0u)
        # save mesh in temporary arrays:
        grd.nxmu=deepcopy(com.nxm)
        grd.nymu=deepcopy(com.nym)
        gchange("Dnull_temp",0)
        grd.rmu=deepcopy(com.rm)
        grd.zmu=deepcopy(com.zm)
        grd.ixpt1u=deepcopy(com.ixpt1)
        grd.ixtopu=deepcopy(grd.ixtop)
        grd.ixpt2u=deepcopy(com.ixpt2)
        grd.iysptrxu=deepcopy(com.iysptrx1)


        ##########################################################################
        ### STEP 4: construct the mesh for the lower half of the configuration ###
        ##########################################################################


        halfmesh(upper=False, nxleg=nxlegl, nxcore=nxcorel, slpxt=slpxtl, nxxpt=nxxptl, 
            nxmod=nxmodl, alfxptu=alfxptul, alfxptl=alfxptll, isrefxptn=isrefxptnl, 
            nsmoothx=nsmoothxl, ismmon=ismmonl, nsmooth=nsmoothl, wtmesh1=wtmesh1l, 
            dmix0=dmix0l)
        # save mesh in temporary arrays:
        grd.nxmb=deepcopy(com.nxm)
        grd.nymb=deepcopy(com.nym)
        gchange("Dnull_temp",0)
        grd.rmb=deepcopy(com.rm)
        grd.zmb=deepcopy(com.zm)
        grd.ixpt1b=deepcopy(com.ixpt1)
        grd.ixtopb=deepcopy(grd.ixtop)
        grd.ixpt2b=deepcopy(com.ixpt2)
        grd.iysptrxb=deepcopy(com.iysptrx1)

        ##################################################
        ### STEP 5: combine the two halves of the mesh ###
        ##################################################

        com.geometry="dnull"

        # allocate space for full double null mesh data
        com.nxpt=2
        gchange("Xpoint_indices",0)
        com.nxm = grd.nxmb + grd.nxmu - 2
        com.nym = grd.nymb
        gchange("RZ_grid_info",0)

        # set separatrix indices (poloidal indices are set below in map subroutines)
        if flx.psi0sep1 < flx.psi0sep2:
           com.iysptrx1[0] = com.nycore[com.igrid-1]
           com.iysptrx2[0] = com.nycore[com.igrid-1] + com.nysol[com.igrid-1]
        else:
           com.iysptrx1[0] = com.nycore[com.igrid-1] + com.nysol[com.igrid-1]
           com.iysptrx2[0] = com.nycore[com.igrid-1]

        com.iysptrx1[1] = com.iysptrx2[0]
        com.iysptrx2[1] = com.iysptrx1[0]

        # NOTE: must map bottom half first, then top half, to get indexing correct:

        # map dnbot data onto full mesh
        grd.mapdnbot()
        # map dntop data onto full mesh
        grd.mapdntop()

        ################################################################
        ### STEP 6: define guard cells at the divertor target plates ###
        ################################################################

        # construct guard cells at target plates
        grd.add_guardc_tp()

        ###########################################################
        ### STEP 7: compute magnetic field values for each cell ###
        ###########################################################

        # get magnetics data
        grd.magnetics(0,com.nxm+1,1,com.nym)

        #########################################
        ### STEP 8: write out the gridue file ###
        #########################################

        # write out gridue file
        try:
            if (com.isgriduehdf5 == 1):
                write_gridue()
            else:
                grd.writednf("gridue", '')
        except:
            grd.writednf("gridue", '')


    def halfmesh(upper=True, nxleg=[9,11], nxcore=[8,10], slpxt=1.2, nxxpt=1, 
        nxmod=4, alfxptu=0.5, alfxptl=1, isrefxptn=1, nsmoothx=3, ismmon=3,
        nsmooth=3, wtmesh1=0, dmix0=0.5):
        ''' Constructs upper dnull mesh '''
        from copy import deepcopy
        from numpy import array
        from Forthon import gchange
        from gridgen import plot_efit, plot_flx
        from uedge import com, bbb, flx, grd
        # This is top half of double-null
        com.geometry="dnbot"
        flx.neqdsk()	# read simagx and sibdry from EFIT
        # Use un-normalized psi values in psitop and psibot
        for n in range(1,com.nym+2):
            flx.psitop[n-1] = com.simagx + (com.sibdry-com.simagx)*\
                flx.psi0_mp_inner[com.nym-n+1]
            flx.psitop[com.jdim-n] = com.simagx + (com.sibdry-com.simagx)*\
                flx.psi0_mp_outer[com.nym-n+1]
    #	   psitop(n) = simagx+(sibdry-simagx)*psi0_mp_inner(nym-n+1)
    #	   psitop(jdim-n+1) = simagx+(sibdry-simagx)*psi0_mp_outer(nym-n+1)
        flx.psitop[com.nym+1] = com.simagx

        if upper is True:
            for n in range(1, com.nym+2):
                flx.psibot[n-1] = com.simagx + (com.sibdry-com.simagx)*flx.psi0_dp_upper_inner[com.nym-n+1]
                flx.psibot[com.jdim-n] = com.simagx + (com.sibdry-com.simagx)*\
                    flx.psi0_dp_upper_outer[com.nym-n+1]
        else:
            for n in range(1, com.nym+2):
                flx.psibot[n-1] = com.simagx + (com.sibdry-com.simagx)*flx.psi0_dp_lower_inner[com.nym-n+1]
                flx.psibot[com.jdim-n] = com.simagx + (com.sibdry-com.simagx)*\
                    flx.psi0_dp_lower_outer[com.nym-n+1]
        flx.psibot[com.nym+1] = com.simagx

        # For upper half only:
        if upper is True:
            flx.iseqdskr=1
            flx.psitop = -flx.psitop
            flx.psibot = -flx.psibot
        else:
            flx.iseqdskr=0


        # inputs for x-mesh:
        com.nxleg[0] = nxleg
        com.nxcore[0] = nxcore
        grd.kxmesh=1
        # mesh refinement near x-point
        grd.slpxt=slpxt
        com.nxxpt=nxxpt
        grd.nxmod=nxmod
        grd.alfxptu=alfxptu
        grd.alfxptl=alfxptl
        grd.isrefxptn = isrefxptn
        grd.nsmoothx = nsmoothx

        # make mesh non-orthogonal
        com.ismmon=ismmon;
        grd.istream=0;
        grd.iplate=1;
        grd.nsmooth=nsmooth;
        grd.wtmesh1=wtmesh1;
        grd.dmix0=dmix0

        # construct the mesh:
        flx.flxrun()

        # Find the separatrix strike points on the divertor plates
        grd.isspnew=1		#=1: calculate new strike points
                #=0: use rvsin,zvsin,rvsout,zvsout from aeqdsk
        # To calculate strike points, replace (xlim,ylim) arrays from EFIT
        # with simplified plate definitions:
        _xlim = deepcopy(com.xlim)
        _ylim = deepcopy(com.ylim)
        _nlim = deepcopy(com.nlim)

        if upper is True:
            plate = 'upper'
        else:
            plate = 'lower'
        grd.nplate1, rplate1, zplate1 = readplatefile('plate1_{}.dat'.format(plate))
        grd.nplate2, rplate2, zplate2 = readplatefile('plate2_{}.dat'.format(plate))
        gchange("Mmod",0)

        if upper is True:
            zplate1 = 2.*com.zshift - zplate1
            zplate2 = 2.*com.zshift - zplate2

        grd.rplate1 = rplate1
        grd.zplate1 = zplate1

        grd.rplate2 = rplate2
        grd.zplate2 = zplate2

        com.nlim = grd.nplate1
        gchange('Comflxgrd')
        com.xlim = grd.rplate1
        com.ylim= grd.zplate1
        rstrike = array(0.)
        zstrike = array(0.)
        flx.findstrike(com.jsptrx[0], rstrike, zstrike)
        grd.rstrike[0] = deepcopy(rstrike)
        grd.zstrike[0] = deepcopy(zstrike)

        com.nlim = grd.nplate2
        gchange('Comflxgrd')
        com.xlim = grd.rplate2
        com.ylim = grd.zplate2
        flx.findstrike(com.jsptrx[1], rstrike, zstrike)
        grd.rstrike[1] = deepcopy(rstrike)
        grd.zstrike[1] = deepcopy(zstrike)
        com.nlim = _nlim
        gchange('Comflxgrd')
        com.xlim = _xlim
        com.ylim = _ylim

        print(grd.rstrike, grd.zstrike)
        grd.grdrun()

        f = plot_flx()
        f.get_axes()[0].plot(grd.rstrike, grd.zstrike, 'go')



    def readplatefile(filename):
        ''' Reads and returns plate '''
        from numpy import array
        lines = []
        with open(filename) as file:
            for line in file:
                line = line.split('#')[0].strip()
                if len(line) > 0:
                    lines.append(line)
        r = []
        z = []
        for point in lines:
            [rbuff, zbuff] = [float(x.lower().replace('d','e')) for x in \
                point.split()]
            r.append(rbuff)
            z.append(zbuff)
        return len(lines), array(r), array(z)


    def uppersn():
        ''' Input file for the UEDGE single-null DIII-D case, shot 160299

        The shot is a set of well-characterized L-mode shots with fwd-B, LSN

        Power flow:
            Pinj=?
                Poh=? MW
            Pwdot=?
            P_rad~?
            P_sep=? MW
            P_div=?
        Boundary condition:
            ncore=3.0e19
            pcoree=pcorei=0.45MW
        '''
        from numpy import ones
        import uedge.contrib.input as setup
        from numpy import pi
        from plot import gridue, plot_plates # Import personal grid plotting function
        from gridgen import plot_efit,plot_flx, plot_poloidal_distribution # Import personal grid visualization tools
        # Input package located in UEDGE/pyscripts/contrib/input


        
        ''' Set up EFIT files and geometry '''
        bbb.mhdgeo=1
        com.aeqdskfname="aeqdsk"
        com.geqdskfname="neqdsk"
        com.geometry="uppersn"

        # Plot read EFIT
    #    plot_efit()

        # I am not 100% sure what these do, but without them we segfault...
        flx.istchkon=1
        flx.dtheta_exclude    = [.75*pi,.50*pi]
        flx.dtheta_overlap_pf = [.05*pi,.01*pi]
        flx.dtheta_overlap_sol=[0.25*pi,0.25*pi]
        
        ''' Set up flux surfaces '''
        # Radial cells
        com.nycore[0] = 6
        com.nysol[0] = 12

    #    com.nycore[0] = 15
    #    com.nysol[0] =16

        # Radial extent of grid in Psi_norm-space
        flx.psi0min1=.95
        flx.psi0min2=.986
        flx.psi0max=1.08

        # Retrieve standard exponentially distributed flux-surfaces
        flx.kymesh=1
        flx.alfcy=2.8
        flx.flxrun()
        plot_flx()

        

        # Manually modify the flux surfaces to conform to the grid nodes
        # Psi_norm of corner: -0.28127680680240885
        # Psi_norm of first corner node: -0.28308801263733624
        # Psi_norm of second corner node: -0.28339239270776795
        # Psi_norm of third corner node: -0.2842713664087835
        # Psi_norm of fourth corner node: -0.2910363664087835
        # Psi_norm of fifth corner node: -0.30365185175200005

        if 1 == 0:
            flx.kymesh=0

            # Set flux surface to hit corner of shelf
            flx.psitop[8] = -0.28308801263733624
            flx.psitop[-9] = -0.28308801263733624

            # Set flux surface to hit corner point
            flx.psitop[12] = -0.28127680680240885
            flx.psitop[-13] = -0.28127680680240885

            dPsi = flx.psitop[12] - flx.psitop[8]
            # Distribute flux-surfaces 6, 7, 8 along the shelf
            flx.psitop[9] = flx.psitop[8] + 0.25*dPsi
            flx.psitop[-10] = flx.psitop[8] + 0.25*dPsi
            
            flx.psitop[10] = flx.psitop[8] + 0.5*dPsi
            flx.psitop[-11] = flx.psitop[8] + 0.5*dPsi
            
            flx.psitop[11] = flx.psitop[8] + 0.75*dPsi
            flx.psitop[-12] = flx.psitop[8] + 0.75*dPsi
            
            # Set a flux-surface at second corner point
            flx.psitop[7] = -0.28339239270776795
            flx.psitop[-8] = -0.28339239270776795

            # Set a flux-surface at the third corner point
            flx.psitop[6] = -0.2842713664087835
            flx.psitop[-7] = -0.2842713664087835

            # Set a flux-surface at the fourt corner point
            # - account for two flux-surfaces in the interior
            #   to maintain equal spacing
            flx.psitop[3] = -0.2910363664087835
            flx.psitop[-4] = -0.2910363664087835

            # Set the interior flux-surfaces
            dPsi = flx.psitop[6] - flx.psitop[3]
            flx.psitop[4] = flx.psitop[3] + 0.33*dPsi
            flx.psitop[-5] = flx.psitop[3] + 0.33*dPsi

            flx.psitop[5] = flx.psitop[3] + 0.66*dPsi
            flx.psitop[-6] = flx.psitop[3] + 0.66*dPsi

            # Finally, set the outermost flux-tube to hit the closest corner
            flx.psitop[0] = -0.30365185175200005
            flx.psitop[-1] = -0.30365185175200005

            # Space the remaining flux-surfaces uniformly
            dPsi = flx.psitop[3] - flx.psitop[0]
            flx.psitop[1] = flx.psitop[0] + 0.33*dPsi
            flx.psitop[-2] = flx.psitop[0] + 0.33*dPsi
            
            flx.psitop[2] = flx.psitop[0] + 0.66*dPsi
            flx.psitop[-3] = flx.psitop[0] + 0.66*dPsi
        

        ''' Poloidal cell distribution '''

        # Set up poloidal cell distribution
        com.nxleg[0]=[8, 12]
        com.nxcore[0]=[10, 10]

    #    com.nxleg[0] = [12, 10]
    #    com.nxcore[0] = [14,14]

        
        # Start by looking/modifying an orthogonal grid
        com.ismmon=0

        # Modify poloidal cell distribution
        grd.kxmesh=2
        grd.slpxt=1.3
        grd.nxgas = [3,5]
    #    grd.nxgas = [3,8]
        grd.dxgas = [3e-3, 2e-3]
        grd.alfx = [0.7, 0.7]
        plot_poloidal_distribution()
        # Visualize grid
    #    flx.flxrun()
    #    grd.grdrun()
    #    gridue(zoom='divertor')

        # Switch to nonorthohonal grid
        grd.iplate=1;
        plate_d3d_usn()
        grd.nsmooth = 1

        
    #    com.ismmon=1
    #    flx.flxrun()
    #    grd.grdrun()
    #    plot_plates(gridue(zoom='divertor'))

    #    com.ismmon=2
    #    flx.flxrun()
    #    grd.grdrun()
    #    plot_plates(gridue(zoom='divertor'))

        # Use 50% surface compression and gradual mixing length
        com.ismmon=3
        grd.wtmesh1 = 0.5
        grd.dmix0 = 1

        flx.flxrun()
        plot_flx()
        grd.grdrun()
        plot_plates(gridue(zoom='divertor').get_axes()[0])

        return

    # mesh construction -- #generate non orthogonal grid
        #plate_d3d_baf1()
    #    bbb.isybdryog=1
    #    com.isnonog = 1

        # Add an additional cell at the X-point?
    #    com.nxxpt=0
    #    grd.nxmod=2
    #    grd.alfxptu=0.5
    #    grd.alfxptl=0.5



