

class ConvergeStep():
    def __init__(self):
        return

    def conv_step(self, increment, name, var, ivar=None, stop=None, 
            b0=False, **kwargs):
        from copy import deepcopy
        from numpy import ndarray
        
        if isinstance(self.get(var), ndarray) and (ivar is None):
            print('No index for var specified. Aborting!')
            return
        if stop is None:
            try:
                stop = self.get(var)[ivar]
            except:
                stop = deepcopy(bbb.__getattribute__(var))
            if increment < 0:
                stop /= 5
            else:
                stop *= 5

        try:
            self.get(var)[ivar]
            varstr = '{}[{}]'.format(var,ivar)
        except:
            varstr = var
        while True:    
            _var = self.getue(var, cp=False) 
            if (b0 is True) and (var == 'b0'):
                currval = 1/deepcopy(_var)
                currval += increment
                currval = 1/currval
                self.setue(var, currval)
            else:
                try:
                    _var[ivar] += increment
                    currval = self.getue(var)[ivar]
                except:
                    _var += increment
                    currval = self.getue(var)
                
            _label = self.getue('label', cp=False)
            _label[0] = '{}_{}={:.3e}'.format(name, var, currval)
            print('======================================')
            print('Solving for {}={:.2E}'.format(varstr, currval))
            print('======================================')
            try:
                self.setue('dtreal', kwargs['dtreal'])
            except:
                self.setue('dtreal', 1e-9)
            self.setue('issfon', 1)
            self.setue('ftol', 1e-5)
            self.exmain()
            self.converge(**kwargs)
            if self.getue('iterm') != 1:
                break
            if (increment > 0) and (currval > stop):
                    break
            elif (increment < 0) and (currval < stop):
                    break


    def conv_step_ncore(self, increment, name, iisp=0, stop=None, **kwargs):
        self.conv_step(increment, name, 'ncore', ivar=iisp, stop=None, 
            **kwargs)

    def conv_b0(self, name, increment=0.03, stop=1, **kwargs):
        self.conv_step(increment, name, 'b0', stop=stop, b0=True)


