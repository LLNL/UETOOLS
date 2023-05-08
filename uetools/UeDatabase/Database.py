from uetools import Case

class Database():
    def __init__(self, databasename, dbidentifier='_UeDB', sortvar='ne',
        sortlocation='midplane', rerun=False):
        from os import getcwd
        from os.path import isfile, isdir
        self.cwd = getcwd()
        self.cases = {}
        self.dbidentifier = dbidentifier
        self.datbasename = databasename
        self.sortvar = sortvar
        self.rerun = rerun
        self.sortlocation = sortlocation
        if isfile(databasename):
            self.restore(databasename)
        else:
            self.create_database(databasename)
        # TODO: Store commonly used grid locations
        # TODO: Account for different grids
        self.ixmp = self.get('ixmp')[0]
        self.iysptrx = self.get('iysptrx')[0]

        # Make sort location more advanced
        if self.sortlocation == 'midplane':
            self.sortlocation = (self.ixmp, self.iysptrx+1)
        self.scanvar = self.get(sortvar)[:, self.sortlocation[0], 
            self.sortlocation[1]]
        # Sort here
        # Save if requested

    def top(self):
        self.chdir(self.cwd)

    def chdir(self, path):
        from os import chdir
        chdir(path)

    def is_case(self, file):
        from h5py import File
        try:
            f = File(file, 'r')
            ret = True
            # Check that necessary groups are present in file
            for entry in ['centered', 'staggered', 'restore', 'grid', 'setup']:
                if entry not in f.keys():
                    ret = False
            f.close()
            return ret
        except:
            return False

    def save(self, savename):
        """ Save the database """
        print('TBD')
    
    def restore(self, restorename):
        """ Restore a database """
        print('TBD')

    def sort(self, variable, location, increasing = True):
        """ Sorts cases according to variable at location """
        print('TBD')

    def create_database(self, path):
        from os import walk
        from os.path import abspath
        from tqdm import tqdm
        path = abspath(path)
        self.chdir(path)
        createdb = []
        for parent, dirs, files in walk('.'):
            subdir = parent.split('/')[-1]
            databases = [db for db in files if self.is_case('{}/{}'.format(\
                parent, db))]
            if self.rerun is False:
                # HDF5s identified
                if len(databases) == 1:
                    # Verify all necessary data is present
                    self.cases['{}/{}'.format(subdir, databases[0].replace(\
                        '.hdf5', ''))] = Case('{}/{}'.format(parent, 
                        databases[0]), inplace=True)
                elif len(databases) > 1:
                    for db in databases:
                        self.cases['{}/{}'.format(subdir, db.replace('.hdf5', 
                        ''))] = Case('{}/{}'.format(parent, db), inplace=True)
                # No database found, store location where input is
                # Changing dirs while executing walk breaks the call
                elif 'input.yaml' in files:
                    createdb.append(parent)
            else:
                if 'input.yaml' in files:
                    createdb.append(parent)
        # Now, create and read the files
        if len(createdb)>0:
            for newdbfolder in tqdm(createdb):
                subdir = newdbfolder.split('/')[-1]
                self.chdir(path)
                self.chdir(newdbfolder)
                Case('input.yaml', verbose=False).save('{}{}.hdf5'.format(subdir, 
                    self.dbidentifier))
                self.cases['{}/{}'.format(subdir, '{}{}'.format(subdir, 
                    self.dbidentifier))] = Case('{}{}.hdf5'.format(subdir, 
                    self.dbidentifier), inplace=True)
                self.chdir(path)
        self.top() 

    def get(self, var, **kwargs):
        """ Returns an array of var for all cases """
        # TODO: Supress verbose output
        from numpy import array
        ret = []
        for key, case in self.cases.items():
            ret.append(case.get(var, **kwargs))
        return array(ret)


    def plotscan(self, var, location, **kwargs):
        self.plotvar(self.scanvar, self.get(var)[:, location[0], location[1]],
            **kwargs)
        
        
        return

    def plotvar(self, xvar, yvar, ax=None, **kwargs):
        """ Plots yvar as a function of xvar for all cases """
        from matplotlib.pyplot import subplots
        if ax is None:
            f, ax = subplots()

        ax.plot(xvar, yvar, **kwargs)

        return ax.get_figure() 

    
    def plot_series(self):
        """ Returns a series of figures to scroll through """
        print('TBD')
    
    def animation(self):
        """ Creates an animation from a series of figures """ 
        print('TBD')
