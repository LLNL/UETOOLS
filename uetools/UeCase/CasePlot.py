# Plotting routines for case objects
from uetools.UePlot import Plot

class Caseplot(Plot):
    def __init__(self):
        super(Caseplot, self).__init__()
        # Calculate distances along targets
        self.otdistance = self.radialdistance(-2)
        self.itdistance = self.radialdistance(1)
        return


    # TODO: implement profile plots


    # TODO: fix this function 
    def radialdistance(self, ind):
        from numpy import cumsum
        # Coordinates of edge of grid
        dr0 = self.get('rm', 2)[ind, 1]
        dz0 = self.get('zm', 2)[ind, 1]
        # Halfway point of cell centers
        drplatecenter = 0.5*(self.get('rm', 2)[ind] + self.get('rm', 4)[ind])
        dzplatecenter = 0.5*(self.get('zm', 2)[ind] + self.get('zm', 4)[ind])
        # Calculate strike point coordinates
        iysptrx = self.get('iysptrx')
        rsp = self.get('rm', 4)[ind, iysptrx] - drplatecenter[iysptrx]
        zsp = self.get('zm', 4)[ind, iysptrx] - dzplatecenter[iysptrx]
        # Normalize to grid edge
        drplatecenter -= dr0
        dzplatecenter -= dz0
        # Calculate cumulative distance along target
        dist = cumsum((drplatecenter**2 + dzplatecenter**2)**0.5)
        # Calculate distance along plate to strike point
        dsep = dist[iysptrx] + (rsp**2 + zsp**2)**0.5
        return dist-dsep

    def it(self, variable, marksep=True, staggered=False, **kwargs):
        # Exchange YYC for working radialdistance
        fig = self.plotprofile(self.itdistance[1:-1], 
            variable[0**staggered, 1:-1], **kwargs)
        # Add Sep location if requested
        if marksep is True:
            fig.get_axes()[0].axvline(0, color='g', linewidth=1)
        return fig

    def ot(self, variable, marksep=True, **kwargs):
        # Exchange YYC for working radialdistance
        fig = self.plotprofile(self.otdistance[1:-1], 
            variable[-2, 1:-1], **kwargs)
        # Add Sep location if requested
        if marksep is True:
            fig.get_axes()[0].axvline(0, color='g', linewidth=1)
        return fig

    def omp(self):
        return

    def imp(self):
        return

    def row(self):
        return

    def neOT(self, s=None, **kwargs):
        return self.ot(self.get('ne', s), **kwargs)

    def teOT(self, s=None, **kwargs):
        return self.ot(self.get('te', s)/self.get('ev'), **kwargs)

    def niOT(self, s=None, **kwargs):
        return self.ot(self.get('ni', s), **kwargs)

    def tiOT(self, s=None, **kwargs):
        return self.ot(self.get('ti', s)/self.get('ev'), **kwargs)

    # Expand the 2D plot list
    def grid(self, **kwargs):
        return self.plotmesh(**kwargs)

    def heatmap(self, var, s=None, **kwargs):
        return self.plotmesh(self.get(var, s), **kwargs)

    def ne2D(self, s=None, **kwargs):
        return self.plotmesh(self.get('ne', s), **kwargs)

    def te2D(self, s=None, **kwargs):
        return self.plotmesh(self.get('te', s)/self.get('ev'), **kwargs)

    def ni2D(self, s=None, **kwargs):
        return self.plotmesh(self.get('ni', s), **kwargs)

    def ti2D(self, s=None, **kwargs):
        return self.plotmesh(self.get('te', s)/self.get('ev'), **kwargs)

    def ng2D(self, s=None, **kwargs):
        self.plotmesh(self.get('ng', s), **kwargs)

    def tg2D(self, s=None, **kwargs):
        return self.plotmesh(self.get('tg', s)/self.get('ev'), **kwargs)

    
