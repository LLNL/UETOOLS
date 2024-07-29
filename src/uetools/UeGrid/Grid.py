class Grid:
    """ Helper class, still only offering plotting routines """
    def __init__(self, case):
        self.get = case.getue
        self.setue = case.setue
        self.plot = GridPlot(self)

    def pick_aeqdskdata(self,geqdsk, ncontour=250, interpres=2000, **kwargs):
        from uedge import com
        from scipy.interpolate import griddata
        from scipy.optimize import fmin
        from numpy import linspace, array, gradient, sum, meshgrid, where
        from matplotlib.pyplot import subplots, ginput, waitforbuttonpress

        f, ax = subplots(figsize=(7,9))
        ax.plot([],[], 'ko')
        ax.plot([],[], 'ro')
        ax.plot([],[], 'bo')
        self.plot.efit(geqdsk, ax=ax, ncontour=ncontour, **kwargs)

        rorig, zorig = self.get('nxefit'), self.get('nyefit')
        xf = lambda r: linspace(0, self.get('xdim'), r)+self.get('rgrid1')
        yf = lambda r: linspace(0, self.get('zdim'), r)-(self.get('zdim')*0.5-self.get('zmid'))
        xo, yo = meshgrid(xf(rorig), yf(zorig))
        fold = self.get('fold').transpose()
        grad = gradient( fold, xo[0,1]-xo[0,0], yo[1,0]-yo[0,0])

        xi, yi = meshgrid(xf(interpres), yf(interpres))
        interp = griddata( (xo.ravel(), yo.ravel()), fold.ravel(), (xi,yi))
        gradinterp = griddata( (xo.ravel(), yo.ravel()), array((grad[0]**2+grad[1]**2)**0.5).ravel(), (xi, yi))

#        c = ax.pcolormesh(xi, yi, gradinterp, cmap='hot_r', vmax=0.1)#, vmin=-0.2, vmax=0.2)
        #c = ax.pcolormesh(xi, yi, (gradinterp[0]**2+gradinterp[1]**2)**0.5, cmap='bwr', vmin=-0.2, vmax=0.2)
#        f.colorbar(c, ax=ax)

        

        print('Manually identify the following points in order')
        print('(Choose by clicking in the figure)')#, undo by right-clicking)')
        textbox = f.get_axes()[0].text(.87, 1.5, "", backgroundcolor='w', zorder=10)
        boxtext = 'Please click on the:\n{}'

        pts = []
        gradmin = []
        psimin = []
        ptlabels = [' - Magnetic axis', ' - Lower X-point', ' - Upper X-point']
        colors = ['b','r','m']
        for i in range(3):
            print(ptlabels[i])
            textbox.set_text(boxtext.format(ptlabels[i])) 
            pt = ginput(1,0)
            nearest = [abs(xf(interpres)-pt[0][0]).argmin(), abs(yf(interpres)-pt[0][1]).argmin()]
            gradmin_nearest = gradinterp[nearest[1], nearest[0]]
            psimin_nearest = interp[nearest[1], nearest[0]]
            gradmin_optimum = gradinterp[
                        nearest[1]-int(interpres/50):nearest[1]+int(interpres/50),  
                        nearest[0]-int(interpres/50):nearest[0]+int(interpres/50)
            ].min()
            y, x = where(gradinterp == gradmin_optimum)
            psimin_optimum = interp[y[0], x[0]]
            if psimin_optimum < psimin_nearest:
                gradmin.append(gradmin_optimum)
                psimin.append(psimin_optimum)
                pts.append([xf(interpres)[x[0]], yf(interpres)[y[0]]])
            else:
                gradmin.append(gradmin_nearest)
                psimin.append(psimin_nearest)
                pts.append([xf(interpres)[nearest[0]], yf(interpres)[nearest[1]]])
            cross = ax.plot(*pts[-1], '+', color=colors[i])
            contour = ax.contour(xf(rorig), yf(zorig), fold, [psimin[-1]], colors=colors[i], 
            linewidths=1.5, linestyles='-') 

            try:
                textbox.set_text(boxtext.format(ptlabels[i+1])) 
            except:
                textbox.set_text("All points defined!") 
            ax.draw(f.canvas.renderer)


#        magx, lxpt, uxpt = ginput(3,0)
#        textbox.set_visible(False)
        
        # Reconstruct EFIT grid
        

        return textbox
    


class GridPlot:
    """ Initial class with routines useful to grid construction """
    
    def __init__(self, grid):
        self.get = grid.get
        self.setue = grid.setue
        pass

    def flx(self, ax=None, surfaces=None):
        ''' Plots flux surfaces as defined in grid setup
            Based on plotflx.bas by Rensink/Rognlien/Porter 

        Keyword arguments:
        ==================
        ax [None] - axis to plot on. If None, creates new figure
        surfaces [None] - range-object with surfaces to plot. If none, 
                            plots all
    
        Returns:
        ========
        matplotlib.pyplot.Figure
        '''
        from Forthon import packageobject
        from matplotlib.pyplot import subplots, Figure, Axes

        packageobject('flx').__getattribute__('flxrun')
        if ax is None:
            f, ax = subplots(figsize=(7,9))
        elif isinstance(ax, Axes):
            pass
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        else:
            raise TypeError("Axes type {} not compatible".format(type(ax)))

        mxs = 2*(self.get('nycore')[0]+self.get('nysol')[0]+self.get('nyout')[0])+2
        if surfaces is None:
            rng = range(mxs+1)
        elif not isinstance(surfaces, range):
            raise TypeError("surfaces must be type 'range'")
        else:
            rng =surfaces
            if rng.stop > mxs:
                rng = range(rng.start, mxs)
            if rng.start<0:
                rng = range(0, rng.stop)
        # Plot vessel if present in EQDSK files
        if self.get('nlim') > 0:
            ax.plot(self.get('xlim'), self.get('ylim'), 'k-', linewidth=3)
            ax.plot(self.get('xlim'), self.get('ylim'), '-', 
                    linewidth=1.5, color='yellow'
            )
        # Plot target plates, if specified
        try:
            ax.plot(self.get('rplate1'), self.get('zplate1'), 'ro-')
        except:
            pass
        try:
            ax.plot(self.get('rplate2'), self.get('zplate2'), 'b.-')
        except:
            pass

        jmin = self.get('jmin')
        jmax = self.get('jmax')
        jsptrx = self.get('jsptrx')
        xcurve = self.get('xcurve')
        ycurve = self.get('ycurve')
        ijumpf = self.get('ijumpf')

        # Plot the flux surfaces within the sepcified range
        for i in rng:
            # Plot SOL flux surfaces for each half-mesh
            if ((i >= jmin[0]-1) and (i < jsptrx[0])) or \
                ((i >= jsptrx[1]) and (i <= jmax[1])):
                ax.plot(xcurve[:,i][abs(xcurve[:,i])>0], 
                        ycurve[:,i][abs(ycurve[:,i])>0], 
                        'k-', linewidth=0.3
                )
            # Plot CORE/PFR flux surfaces for each half-mesh
            elif ((i >= jsptrx[0]) and (i <= jmax[0])) or \
                    ((i >= jmin[1]-1) and (i <= jsptrx[1]+1)):
                ax.plot(xcurve[:,i][abs(xcurve[:,i])>0][:ijumpf[i]], 
                        ycurve[:,i][abs(ycurve[:,i])>0][:ijumpf[i]], 
                        'k-', linewidth=0.3
                )
                ax.plot(xcurve[:,i][abs(xcurve[:,i])>0][ijumpf[i]:], 
                        ycurve[:,i][abs(ycurve[:,i])>0][ijumpf[i]:], 
                        'k-', linewidth=0.3
                )
            for i in [jsptrx[0]-1, jsptrx[1]-1]:
                ax.plot(xcurve[:,i][abs(xcurve[:,i])>0], 
                        ycurve[:,i][abs(ycurve[:,i])>0], 
                        'r-', linewidth=0.5
                )
        ax.set_aspect('equal')
        ax.set_xlabel('Horizontal position [m]')
        ax.set_ylabel('Vertical position [m]')
        return ax.get_figure()


    def efit(self, geqdsk, aeqdsk=None, ax=None, ncontour=80, 
        color='grey', linestyle='solid', sepcolor='k', linewidth=0.5):
        """ Function to plot EFIT contours

        """
        from matplotlib.pyplot import subplots, Figure, Axes
        from copy import deepcopy
        from numpy import linspace
        from uedge import com, flx
        from os.path import exists
        from scipy.interpolate import interp2d
        from Forthon import packageobject
        # Backup original pointers
        oldaeqdskfname = deepcopy(self.get('aeqdskfname'))
        oldgeqdskfname = deepcopy(self.get('geqdskfname'))
        # Set new file paths
        self.setue('geqdskfname', geqdsk)

        # Check whether the aeqdsk file can be located: if not, do not execute aeqdsk()
        if aeqdsk is not None:
            if exists(aeqdsk):
                self.setue('aeqdskfname', aeqdsk)
                packageobject('flx').__getattribute__('aeqdsk')()

        if exists(geqdsk):
            packageobject('flx').__getattribute__('neqdsk')()
        else:
            raise FileNotFoundError(    \
                'EFIT geqdsk file "{}" not found.\nAborting...'.format(\
                geqdsk)
            )

        if ax is None:
            f, ax = subplots(figsize=(7,9))
        elif isinstance(ax, Axes):
            pass
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        else:
            raise TypeError("Axes type {} not compatible".format(type(ax)))

        # Reconstruct EFIT grid
        x = linspace(0, self.get('xdim'), self.get('nxefit'))+self.get('rgrid1')
        y = linspace(0, self.get('zdim'), self.get('nyefit'))-(self.get('zdim')*0.5-self.get('zmid'))
        
        


        fold = self.get('fold').transpose()
        # Interpolate on EFIT grid to find X-points
        interp = interp2d(x, y, fold)

        ax.contour(x, y, fold, ncontour, colors=color, 
            linewidths=linewidth, linestyles=linestyle) 

        rseps2 = self.get('rseps2')
        zseps2 = self.get('zseps2')
        # Check whether the upper X-point exists
        if (x.min() <= rseps2 <= x.max()) and \
            (y.min() <= zseps2 <= y.max()):
                upperxpoint = interp(rseps2, zseps2)
                ax.contour(x, y, fold, [upperxpoint], 
                    colors=sepcolor, linewidths=1, linestyles='solid') 
        
        rseps = self.get('rseps')
        zseps = self.get('zseps')
        # Check whether the lower X-point exists
        if (x.min() <= rseps <= x.max()) and \
            (y.min() <= zseps <= y.max()):
                lowerxpoint = interp(rseps, zseps)
                ax.contour(x, y, fold, [lowerxpoint], 
                    colors=sepcolor, linewidths=1, linestyles='solid') 

        ax.plot(self.get('xlim'), self.get('ylim'), 'k-', linewidth=2)
        try:
            ax.plot(self.get('xlim'), self.get('ylim'), 'k-', linewidth=2)
        except:
            pass
        ax.set_aspect('equal')
        ax.set_xlabel('Horizontal position [m]')
        ax.set_ylabel('Vertical position [m]')
        ax.set_title(self.get('runid')[0].decode('UTF-8').strip())

        # Restore original pointers
        self.setue('aeqdskfname', oldaeqdskfname)
        self.setue('geqdskfname', oldgeqdskfname)

        return ax.get_figure()


