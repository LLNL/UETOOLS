
class Surfaces:
    def __init__(self, case):
        self.get = case.get
        self.getue = case.getue
        self.setue = case.setue
        self.pump_area = []
        
    def plot_pump_areas(self, ax=None, **kwargs):
        from matplotlib.pyplot import subplots, Axes, Figure
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Axes):
            f = ax.get_figure()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
            f = ax.get_figure()
        else:
            raise TypeError(f"ax type {type(ax)} not compatible!")

        for polygon in self.pump_area:
            ax.plot(*polygon.exterior.xy, **kwargs)

    def geometric_pump(self, points, albedo, reference_area=None, albedo_bg=1, **kwargs):
        """ Geometric definition of pumping boundary surfaces
        
        points - list of points defining the pumping surface:
                 intersecting boundary cells are pumped by albedo
        
        """
        from shapely import Polygon
        # Define area to pump
        pump_area = Polygon(points)
        self.pump_area.append(pump_area)
        
        # Define boundary arrays and parameters
        boundaries = {
            'it': {
                'alb': 'alblb',
                'recyc': 'recylb_use'
            },
            'ot': {
                'alb': 'albrb',
                'recyc': 'recyrb_use'
            },
            'pfr': {
                'alb': 'albedoo',
                'recyc': 'recypf_use'
            },
            'wall': {
                'alb': 'albedoi',
                'recyc': 'recywall_use'
            },
        }

        # User-defined albedos
        self.setue('albedo_by_user', 1)
        self.lfs_target_pump(pump_area, albedo, reference_area, albedo_bg, **kwargs)
        
    def lfs_target_pump(self, polygon, albedo, reference_area=None, albedo_bg=1, reset=True):
        from shapely import LineString
        print("Setting LFS target pumping by albedos and recycling")
        alb = self.getue('albrb', copy=True)
        recyc = self.getue('recyrb_use', copy=True)
        if reset:
            # TODO: separate albedo and recyc? Separate species?
            alb[:] = albedo_bg
            recyc[:,0,0] = albedo_bg - 1
        rm = self.getue('rm')
        zm = self.getue('zm')
        intersects = []
        intersect_fraction = []
        for iy in range(1, self.get('ny')+1):
            surf = LineString([(rm[-2,iy,2], zm[-2,iy,2]), (rm[-2,iy,4], zm[-2,iy,4])])
            if polygon.contains(surf):
                intersects.append(iy)
                intersect_fraction.append(1)
                alb[iy] = albedo
                recyc[iy,0,0] = albedo - 1
            elif polygon.intersects(surf):
                frac = polygon.intersection(surf).length / surf.length
                if frac < 1.e-6:
                    continue
                intersects.append(iy)
                intersect_fraction.append(frac)
                alb[iy] = albedo * frac
                recyc[iy,0,0] = (albedo - 1) * frac
        if reference_area is not None:
            total_area = 0
            sx = self.get('sx')
            for i in range(len(intersects)):
                iy = intersects[i]
                frac = intersect_fraction[i]
                total_area += sx[-2, iy]*frac
            frac = reference_area/total_area
            for iy in intersects:
                alb[iy] = 1 - (1-alb[iy])*frac
                recyc[iy,0,0] *= frac
                for s in range(alb.shape[1]):
                    alb[iy, s] = 1 + recyc[iy, 0]
#                alb[iy] = 1 + recyc[iy]
        self.setue('albrb', alb)
        self.setue('recyrb_use', recyc)
                
            
                
