

class Misc():
    def smooth_curve(self, x, y, s=1, **kwargs):
        from scipy.interpolate import splrep, BSpline
        return x, BSpline(*splrep(x, y, s=s))(x)
        
        
    def get_bottomkeys(self,vardict, keys=None):
        """ Returns a list of keys at the bottom of nested dict """
        for key, var in vardict.items():
            if isinstance(var, dict):
                keys = self.get_bottomkeys(var, keys)
            # Hit bottom of dict
            elif isinstance(key, str): 
                try:
                    keys.append(key)
                except:
                    keys = []
                    keys.append(key)
        return keys


