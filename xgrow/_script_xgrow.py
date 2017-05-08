import os
import sys
import re
from collections import namedtuple
from xgrow import xgrow

TF_KEYS = ['wander', 'pause', 'movie', 'zero_bonds']

parsed_traditional = {
    'no_fission': lambda x: ('fission', 'off'),
    'fission': lambda x: ('fission', 'on'),
    'chunk_fission': lambda x: ('fission', 'chunk'),
    'importfile': lambda x: ('importfile', x.group(3)
                             if x.group(3) else 'xgrow_export_output') }

parsed_new = {
    'w': lambda x: ('window',not bool(x.group(1))),
    'nw': lambda x: ('window',bool(x.group(1))) }

def main():
    """\
The new Python xgrow wrapper, designed to run both older xgrow tilesets with
xgrow command line arguments, and stxg tilesets originally run by xgrow-wrap.
Note that this does *not* currently run Alhambra tilesets: that functionality
will likely be moved into Alhambra itself.
"""
    # Why doesn't this use argparse or click? Because, in order to maintain
    # backward compatibility, it can't... What we want is the following:
    # 1. Tileset must come before non-dashed options, but otherwise can be anywhere.
    # 2. Dashed options can go before or after the tileset.
    # 3. Dashed options can take split or = form.
    # 4. Undashed options must take xgrow (=) form.

    raw_args = sys.argv[1:]
    args = {}
    pos = 0
    tilepath = None
    while pos < len(raw_args):
        # First try to parse as an argument:
        m = re.match(r'-{1,2}(no[-_])?([^=\s]+)(=)?(\S+)?',raw_args[pos])
        if m and m.group(2) in parsed_new.keys():
            args.__setitem__(*parsed_new[m.group(2)](m))
            pos += 1
        elif m and m.group(4):
            args[m.group(2)] = m.group(4)
            pos += 1
        elif m:
            if m.group(2) not in TF_KEYS:
                args[m.group(2)] = raw_args[pos+1]
                pos += 2
            else:
                if m.group(1):
                    args[m.group(2)] = False
                else:
                    args[m.group(2)] = True
                pos += 1
        elif not tilepath:
            tilepath = raw_args[pos]
            pos += 1
        else:
            m = re.match(r'([^=\s]+)(=)?(\S+)?',raw_args[pos])
            if m.group(1) in parsed_traditional.keys():
                args.__setitem__(*parsed_traditional[m.group(1)](m))
            elif m.group(3):
                args[m.group(1)]=m.group(3)
            else:
                args[m.group(1)]=True
            pos += 1

    fd = open(tilepath,'r')
    for l in fd:
        if l[0]=="%":
            continue
        else:
            if l[0:4]=='tile' or l[0:3]=='num':
                filetype='xgrow'
            else:
                filetype='stxg'
                from xgrow import stxg
            break
    fd.seek(0)

    if filetype=='xgrow':
        xgrow.run_old( fd.read(), args )
    else:
        xgrow.run( yaml.load(fd), params )
    
if __name__ == '__main__':
    main()
