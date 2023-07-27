
class PostProcessors:
    
    def postproc_helloworld(self):
        print('Hello world!')

    def set_psinormc(self, simagxs, sibdrys):
        self.psinormc =  (self.get('psi')[self.get('ixmp'),:,0]-simagxs)/ \
            (sibdrys-simagxs)

    def set_psinormf(self, simagxs, sibdrys):
        psi = self.get('psi')[self.get('ixmp')]
        self.psinormf = ((0.5*(psi[:,3] + psi[:,4])) -simagxs)/ \
            (sibdrys-simagxs)
