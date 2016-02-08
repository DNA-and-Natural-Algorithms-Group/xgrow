import re
import numpy as np
import string
from io import StringIO
import pandas as pd

def loadflake(xgrowstring,onlytiles=False):
    tiles = np.genfromtxt(
        StringIO(                                       # genfromtxt needs an io object
            re.sub(r'(\[|\])','',                       # remove [ and ]
                re.sub(r'; \.\.\.','',xgrowstring))),   # remove ; ... at line ends
        skip_header=4,                                  # remove first lines
        skip_footer=1,                                  # remove end junk
        dtype='uint'
        )
    if onlytiles:
        return tiles
    data = pd.Series(
        np.genfromtxt(StringIO(xgrowstring.split('\n')[2]))[1:-1],
        index = ['gmc','gse','k','time','tiles','mismatches','events','perimeter','g','dgbonds']
        )
    return {'data': data, 'tiles': tiles}
