

class DB_1DPlots():

    def tiITsep(self, ylabel="$\mathrm{T_{i,sep}^{IT}}~[eV]$", **kwargs):
        return self.plotITsep(self.get("ti") / 1.602e-19, 
                    ylabel=ylabel, **kwargs)

    def tiOTsep(self, ylabel="$\mathrm{T_{i,sep}^{OT}~ [eV]$", **kwargs):
        return self.plotOTsep(self.get("ti") / 1.602e-19, 
                    ylabel=ylabel, **kwargs)


    def teITsep(self, ylabel="$\mathrm{T_{e,sep}^{IT}}~[eV]$", **kwargs):
        return self.plotITsep(self.get("te") / 1.602e-19, 
                    ylabel=ylabel, **kwargs)

    def teOTsep(self, ylabel="$\mathrm{T_{e,sep}^{OT}}~[eV]$", **kwargs):
        return self.plotOTsep(self.get("te") / 1.602e-19, 
                    ylabel=ylabel, **kwargs)

    def neOTsep(self, ylabel="$\mathrm{n_{e,sep}^{OT}}~[m^{-3}]$",**kwargs):
        return self.plotOTsep(self.get("ne"), 
                    ylabel=ylabel, **kwargs)

    def neITsep(self, ylabel="$\mathrm{n_{e,sep}^{IT}}~[m^{-3}]$", **kwargs):
        return self.plotITsep(self.get("ne"), 
                    ylabel=ylabel, **kwargs)

    def niOTsep(self, species, ylabel=None, **kwargs):
        if ylabel is None:
            ylabel="$\mathrm{{n_{{i={},sep}}^{{OT}}~[m^{{-3}}]}}$".format(species)
        return self.plotOTsep(self.get("ni")[:,:,:,species], 
                    ylabel=ylabel, **kwargs)

    def niITsep(self, species, ylabel=None, **kwargs):
        if ylabel is None:
            ylabel="$\mathrm{{n_{{i={},sep}}^{{IT}} [m^{{-3}}]}}$".format(species)
        return self.plotITsep(self.get("ni")[:,:,:,species], 
                    ylabel=ylabel, **kwargs)

    def ngOTsep(self, species, ylabel=None, **kwargs):
        if ylabel is None:
            ylabel="$\mathrm{{n_{{g={},sep}}^{{OT}}~[m^{{-3}}]}}$".format(species)
        return self.plotOTsep(self.get("ng")[:,:,:,species], 
                    ylabel=ylabel, **kwargs)

    def ngITsep(self, species, ylabel=None, **kwargs):
        if ylabel is None:
            ylabel="$\mathrm{{n_{{g={},sep}}^{{IT}}~[m^{{-3}}]}}$".format(species)
        return self.plotITsep(self.get("ng")[:,:,:,species], 
                    ylabel=ylabel, **kwargs)

    def tgOTsep(self, species, ylabel=None, **kwargs):
        if ylabel is None:
            ylabel="$\mathrm{{T_{{g={},sep}}^{{OT}}~[eV]}}$".format(species)
        return self.plotOTsep(self.get("tg")[:,:,:,species]/1.602e-19, 
                    ylabel=ylabel, **kwargs)

    def tgITsep(self, species, ylabel=None, **kwargs):
        if ylabel is None:
            ylabel="$\mathrm{{T_{{g={},sep}}^{{IT}}~[eV]}}$".format(species)
        return self.plotITsep(self.get("tg")[:,:,:,species]/1.602e-19, 
                    ylabel=ylabel, **kwargs)



    def tiITmax(self, ylabel="$\mathrm{T_{i,max}^{IT}}~[eV]$", **kwargs):
        return self.plotITmax(self.get("ti") / 1.602e-19, 
                    ylabel=ylabel, **kwargs)

    def tiOTmax(self, ylabel="$\mathrm{T_{i,max}^{OT}}~[eV]$", **kwargs):
        return self.plotOTmax(self.get("ti") / 1.602e-19, 
                    ylabel=ylabel, **kwargs)


    def teITmax(self, ylabel="$\mathrm{T_{e,max}^{IT}}~[eV]$", **kwargs):
        return self.plotITmax(self.get("te") / 1.602e-19, 
                    ylabel=ylabel, **kwargs)

    def teOTmax(self, ylabel="$\mathrm{T_{e,max}^{OT}}~[eV]$", **kwargs):
        return self.plotOTmax(self.get("te") / 1.602e-19, 
                    ylabel=ylabel, **kwargs)

    def neOTmax(self, ylabel="$\mathrm{n_{e,max}^{OT}}~[m^{-3}]$",**kwargs):
        return self.plotOTmax(self.get("ne"), 
                    ylabel=ylabel, **kwargs)

    def neITmax(self, ylabel="$\mathrm{n_{e,max}^{IT}}~[m^{-3}]$", **kwargs):
        return self.plotITmax(self.get("ne"), 
                    ylabel=ylabel, **kwargs)

    def niOTmax(self, species, ylabel=None, **kwargs):
        if ylabel is None:
            ylabel="$\mathrm{{n_{{i={},max}}^{{OT}}~[m^{{-3}}]}}$".format(species)
        return self.plotOTmax(self.get("ni")[:,:,:,species], 
                    ylabel=ylabel, **kwargs)

    def niITmax(self, species, ylabel=None, **kwargs):
        if ylabel is None:
            ylabel="$\mathrm{{n_{{i={},max}}^{{IT}} [m^{{-3}}]}}$".format(species)
        return self.plotITmax(self.get("ni")[:,:,:,species], 
                    ylabel=ylabel, **kwargs)

    def ngOTmax(self, species, ylabel=None, **kwargs):
        if ylabel is None:
            ylabel="$\mathrm{{n_{{g={},max}}^{{OT}}~[m^{{-3}}]}}$".format(species)
        return self.plotOTmax(self.get("ng")[:,:,:,species], 
                    ylabel=ylabel, **kwargs)

    def ngITmax(self, species, ylabel=None, **kwargs):
        if ylabel is None:
            ylabel="$\mathrm{{n_{{g={},max}}^{{IT}}~[m^{{-3}}]}}$".format(species)
        return self.plotITmax(self.get("ng")[:,:,:,species], 
                    ylabel=ylabel, **kwargs)

    def tgOTmax(self, species, ylabel=None, **kwargs):
        if ylabel is None:
            ylabel="$\mathrm{{T_{{g={},max}}^{{OT}}~[eV]}}$".format(species)
        return self.plotOTmax(self.get("tg")[:,:,:,species]/1.602e-19, 
                    ylabel=ylabel, **kwargs)

    def tgITmax(self, species, ylabel=None, **kwargs):
        if ylabel is None:
            ylabel="$\mathrm{{T_{{g={},max}}^{{IT}}~[eV]}}$".format(species)
        return self.plotITmax(self.get("tg")[:,:,:,species]/1.602e-19, 
                    ylabel=ylabel, **kwargs)


    def tiOMP(self, ylabel="$\mathrm{T_{i,sep}^{OMP}~[eV]}$", **kwargs):
        return self.plotOMP(self.get("ti")/1.602e-19, ylabel=ylabel, **kwargs)
        
    def teOMP(self, ylabel="$\mathrm{T_{e,sep}^{OMP}~[eV]}$", **kwargs):
        return self.plotOMP(self.get("te")/1.602e-19, ylabel=ylabel, **kwargs)
        
    def tgOMP(self, ylabel=None, species=None, **kwargs):
        if ylabel is None:
            ylabel = "$\mathrm{{T_{{g={},sep}}^{{OMP}}~[eV]}$".format(species)
        return self.plotOMP(self.get("tg")[:,:,:,species]/1.602e-19, ylabel=ylabel, **kwargs)

 
    def niOMP(self, ylabel=None, species=None, **kwargs):
        if ylabel is None:
            ylabel = "$\mathrm{{n_{{i={},sep}}^{{OMP}}~[m^{{-3}}]}}$".format(species)
        return self.plotOMP(self.get("ni")[:,:,:,species], ylabel=ylabel, **kwargs)
        
    def neOMP(self, ylabel="$\mathrm{n_{e,sep}^{OMP}~[m^{-3}]}}$", **kwargs):
        return self.plotOMP(self.get("ne"), ylabel=ylabel, **kwargs)
        
    def ngOMP(self, ylabel=None, species=None, **kwargs):
        if ylabel is None:
            ylabel = "$\mathrm{{n_{{g={},sep}}^{{OMP}}~[m^{{-3}}]}}$".format(species)
        return self.plotOMP(self.get("ng")[:,:,:,species], ylabel=ylabel, **kwargs)

               


    # TODO: generalize to account for different grid sizes
    def plotITsep(self, var, **kwargs):
        return self.plotscan(var, (1, self.iysptrx + 1), **kwargs)

    def plotOTsep(self, var, **kwargs):
        return self.plotscan(var, (-2, self.iysptrx + 1), **kwargs)

    def plotITmax(self, var, inds=(None, None), **kwargs):
        from numpy import max

        return self.plotvar(
            self.sortvalues, max(var[:, 1, slice(*inds)], axis=1), **kwargs
        )

    def plotOTmax(self, var, inds=(None, None), **kwargs):
        from numpy import max

        return self.plotvar(
            self.sortvalues, max(var[:, -2, slice(*inds)], axis=1), **kwargs
        )

    def plotOMP(self, var, ylabel=None, **kwargs):
        return self.plotscan(var, (self.ixmp, self.iysptrx + 1), ylabel=ylabel, **kwargs)

    def plotscan(self, var, location=(), **kwargs):
        for index in location:
            var = var[:, index]
        return self.plotvar(self.sortvalues, var, **kwargs)

    def plotvar(self, xvar, yvar, ax=None, xlabel=None, ylabel=None, 
            marker=".", linestyle="", color="k", **kwargs):
        """Plots yvar as a function of xvar for all cases"""
        from matplotlib.pyplot import subplots, Figure

        if ax is None:
            f, ax = subplots()
        elif ax is Figure:
            ax = f.get_axes()[0]

        ax.plot(xvar, yvar, marker=marker, linestyle=linestyle, color=color, 
            **kwargs)
        if xlabel is None:
            ax.set_xlabel(self.sortlabel)
        else:
            ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        return ax.get_figure()

   
