import re
import numpy as np
from io import BytesIO as StringIO
import pandas as pd


def load_array(xgrowstring, onlytiles=False):
    tiles = np.genfromtxt(
        StringIO(
            (re.sub(r'(\[|\])', '',   # remove [ and ]
             re.sub(r'; \.\.\.', '',  # remove ; ... at line ends
                    xgrowstring))
             ).encode()),
        skip_header=4,                # remove first lines
        skip_footer=1,                # remove end junk
        dtype='uint'
        )
    if onlytiles:
        return tiles
    data = pd.Series(
        np.genfromtxt(StringIO((xgrowstring.split('\n')[2]).encode()))[1:-1],
        index=['gmc', 'gse', 'k', 'time', 'tiles', 'mismatches', 'events',
               'perimeter', 'g', 'dgbonds'])
    return {'data': data, 'tiles': tiles}


def load_trace(s):
    data = pd.DataFrame(
        np.genfromtxt(StringIO(s.encode())),
        columns=['gmc', 'gse', 'k', 'time', 'tiles', 'mismatches', 'events',
                 'perimeter', 'g', 'dgbonds'])
    return data


def load_data(s):
    data = pd.Series(
        np.genfromtxt(StringIO(s.encode())),
        index=['gmc', 'gse', 'k', 'time', 'tiles', 'mismatches', 'events',
               'perimeter', 'g', 'dgbonds'])
    return data


def show_array(a, ts, **kwargs):
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from alhambra.tilestructures import xcolors  # FIXME:  put in here!
    mcolors = {n: tuple(z / 255.0 for z in eval(x[3:]))
               for n, x in xcolors.items()}
    cmap = colors.ListedColormap(['black'] + [mcolors[x['color']]
                                              for x in ts['tiles']])
    try:
        plt.imshow(
            a['tiles'], cmap=cmap, vmin=0, vmax=len(ts['tiles']), **kwargs)
    except:
        plt.imshow(a, cmap=cmap, vmin=0, vmax=len(ts['tiles']), **kwargs)
