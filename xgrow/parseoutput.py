import re
from io import BytesIO as StringIO
import numpy as np
import pandas as pd
import pkg_resources

rgbv = pkg_resources.resource_stream(__name__, 'rgb.txt')
xcolors = {" ".join(y[3:]): "rgb({},{},{})".format(y[0], y[1], y[2])
           for y in [x.decode().split() for x in rgbv]}
del rgbv


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


def show_array(a, ts, emptycolor='black', **kwargs):
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    
    mcolors = {n: tuple(z / 255.0 for z in eval(x[3:]))
               for n, x in xcolors.items()}
    cmap = colors.ListedColormap([emptycolor] + [mcolors[x['color']]
                                              for x in ts['tiles']])
    try:
        return plt.imshow(
            a['tiles'], cmap=cmap, vmin=0, vmax=len(ts['tiles']), **kwargs)
    except:
        return (plt.imshow(a, cmap=cmap, vmin=0, vmax=len(ts['tiles']), **kwargs), cmap)
