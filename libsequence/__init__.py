__version__='0.2.2'

from ._libsequence import *

def get_includes():
    """
    Returns absolute path to location of libsequence headers
    """
    import os
    import libsequence
    return os.path.dirname(libsequence.__file__)+'/src/libsequence'

class Windows:
    """
    An iterable list of sliding windows created from a :class:`libsequence.PolyTable`
    """
    def __init__(self, pt, window_size, step_len, starting_pos = 0.,  ending_pos = 1):
        self.windows=[]
        if isinstance(pt,_libsequence.SimData):
            # from .windows_cpp import SimDataWindows
            temp = SimDataWindows(pt,window_size,step_len,starting_pos,ending_pos)
        else:
            # from .windows_cpp import PolySitesWindows
            temp = PolySitesWindows(pt,window_size,step_len,starting_pos,ending_pos)
        for i in range(len(temp)):
            self.windows.append(temp[i])
        temp=None
    def __iter__(self):
        return iter(self.windows)
    def __next__(self):
        return next(self.windows)
    def __getitem__(self,i):
        return self.windows[i]
    def __len__(self):
        return len(self.windows)
