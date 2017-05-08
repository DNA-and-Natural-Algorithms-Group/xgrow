import re
import numpy as np
import string
from io import BytesIO as StringIO
import pandas as pd

def load_array(xgrowstring,onlytiles=False):
    tiles = np.genfromtxt(
        StringIO(                                       # genfromtxt needs an io object
            (re.sub(r'(\[|\])','',                       # remove [ and ]
                    re.sub(r'; \.\.\.','',xgrowstring))).encode()),   # remove ; ... at line ends
        skip_header=4,                                  # remove first lines
        skip_footer=1,                                  # remove end junk
        dtype='uint'
        )
    if onlytiles:
        return tiles
    data = pd.Series(
        np.genfromtxt(StringIO((xgrowstring.split('\n')[2]).encode()))[1:-1],
        index = ['gmc','gse','k','time','tiles','mismatches','events','perimeter','g','dgbonds']
        )
    return {'data': data, 'tiles': tiles}

def load_trace(s):
    data = pd.DataFrame(
    np.genfromtxt(StringIO(s.encode())),
        columns = ['gmc','gse','k','time','tiles','mismatches','events','perimeter','g','dgbonds'])
    return data

def load_data(s):
    data = pd.Series(
    np.genfromtxt(StringIO(s.encode())),
        index = ['gmc','gse','k','time','tiles','mismatches','events','perimeter','g','dgbonds'])
    return data
