

class Misc():
    def smooth_curve(self, x, y, s=1, **kwargs):
        """ Returns smoothed x and y data """
        from scipy.interpolate import splrep, BSpline
        return x, BSpline(*splrep(x, y, s=s))(x)
        
        
    def get_bottomkeys(self, vardict, keys=None):
        """ Returns a list of keys at the bottom of nested dict """
        # Iterate through all dictionary entries
        for key, var in vardict.items():
            # If there is a nested dict, traverse down recursively
            if isinstance(var, dict):
                keys = self.get_bottomkeys(var, keys)
            # Otherwise, ensure that the key is a variable, not index
            elif isinstance(key, str): 
                # Assert there is a list to store to, and store variable
                try:
                    keys.append(key)
                except:
                    keys = []
                    keys.append(key)
            # The YAML can define which indices to set: this catches and
            # passes on any such instances
            else:
                continue
        return keys


