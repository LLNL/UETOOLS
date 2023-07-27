

class Misc():
    def smooth_curve(x, y, s=1, **kwargs):
        from scipy.interpolate import splrep, BSpline
        return x, BSpline(*splrep(x, y, s=s))(x)
        
        


