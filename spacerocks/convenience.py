from astropy.table import Table
import pandas as pd
import copy

class Convenience:

    def __len__(self):
        '''
        This method allows you to use the len() function on a SpaceRocks object.
        '''
        return len(self.name)


    def __getitem__(self, idx):
        '''
        This method allows you to index a SpaceRocks object.
        '''
        p = copy.copy(self)
        for attr in self.__dict__.keys():
            setattr(p, attr, getattr(self, attr)[idx])

        return p

    def astropy_table(self):
        '''
        Write the rocks to an astropy table. This can handle units, though
        it is generally less elegant than pandas.
        '''
        return Table(self.__dict__)


    def pandas_df(self):
        '''
        Write the rocks to a pandas dataframe. Pandas can't handle astropy
        units (yet), so if you want to keep units intact you'll have to use
        an Astropy Table. to_pandas()
        '''
        return self.astropy_table().to_pandas()
