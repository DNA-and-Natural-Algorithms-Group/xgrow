#!/usr/bin/env python

# xgrow-wrap is a wrapper for Xgrow that allows it to use stxg
# (structured xgrow) yaml files for tilesets instead of the usual
# xgrow tile files. It has essentially no other options at the moment,
# and relies on the stxg.py library for stxg support.

# Created in 2013 by Constantine Evans <cge@dna.caltech.edu>

import os
import yaml
import xgrowutils.stxg as stxg
import tempfile
import sys
import argparse

def main():
    
    p = argparse.ArgumentParser(description="The xgrow-wrap wrapper for xgrow, running stxg and tilesetdesigner/yaml files!")

    p.add_argument("-p","--perfect", action="store_true", help="for tilesetdesigner/yaml input, use perfect sticky ends",default=None)
    p.add_argument("-e","--energetics", action="store_false", help="for tilesetdesigner/yaml input, use sequence-dependent sticky ends", dest="perfect")
    p.add_argument("-r","--rotate", action="store_true", help=" tilesetdesigner/yaml input, use sequence-dependent sticky ends", default=False)
    p.add_argument("-f","--norotate", action="store_false", help=" tilesetdesigner/yaml input, use sequence-dependent sticky ends", dest="rotate")
    p.add_argument("-x","--xgrowpath", help="path to xgrow (optional, can also use XGROW_DIR env variable", dest="xgrowpath")
    p.add_argument("-c","--daoecomb", action="store_true", help="if using energetics, use the combined daoe energetics algorithms",default=None)
    p.add_argument("-t","--temperature", help="if using daoecomb, set the temperature for energetics calcs", default=37)
    p.add_argument("tileset", help="an stxg or tilesetdesigner/yaml file to run")
    p.add_argument("xgrow_options", nargs=argparse.REMAINDER, help="options to pass to xgrow")

    args = p.parse_args()

    # Load the stxg file.
    try:
        stream = open(args.tileset)
        xg = yaml.load( stream )
    except IOError as e:
        sys.stderr.write("Error opening tileset file '{0}': {1}\n".format(e.filename, e.strerror))
        stream.close()
        return 1
    except yaml.scanner.ScannerError as e:
        sys.stderr.write("Error loading tileset file '{0}' at line {1}: {2}\n".format(
            argv[1], e.problem_mark.line, e.problem
            ))
        stream.close()
        return 1
    finally:
        stream.close()
        
    # Look for xgrow.
    paths = [ os.path.curdir ]
    if args.xgrowpath:
        paths.append( args.xgrowpath )
    if 'XGROW_DIR' in os.environ:
        paths.append( os.environ['XGROW_DIR'] )
    paths += os.environ['PATH'].split(os.pathsep)
    
    xgrowpath = None
    for path in paths:
        testpath = os.path.join(path.strip('"'),'xgrow')
        if os.path.isfile(testpath):
            xgrowpath = testpath
            break
    if not xgrowpath:
        sys.stderr.write("Can't find xgrow. The xgrow binary should be in the current folder\n"+
                         ", your path, or in the directory specified by XGROW_DIR.\n")
        return 1
    gseavgparam = ""        
    # Now convert from stxg to xgrow.
    try:
        tilestring = stxg.to_xgrow(xg)
    except KeyError:
        sys.stderr.write("Assuming a YAML tileset definition file.\n")
        if args.perfect == None: # The user didn't specify, so let's make an educated guess
            if 'ends' in xg.keys() and 'fseq' in xg['ends'][0].keys():
                sys.stderr.write("Found end sequences. Using them.\n")
                args.perfect=False
            else:
                args.perfect=True
                sys.stderr.write("Assuming perfect ends.\n")
        enc = None
        if args.daoecomb:
            from stickydesign.energetics_daoe import energetics_daoe
            enc = energetics_daoe(temperature=float(args.temperature))
        xgg = stxg.from_yaml_endadj(xg,perfect=args.perfect,rotate=args.rotate,energetics=enc)


        if 'gse_calc_avg' in xgg['xgrowargs'].keys():
            gseavgparam = "Gse={}".format(xgg['xgrowargs']['gse_calc_avg'])
            sys.stderr.write("Found gse_calc_avg of {}.\n".format(xgg['xgrowargs']['gse_calc_avg']))

        tilestring = stxg.to_xgrow(xgg)
    except Exception as e:
        sys.stderr.write("Error converting stxg file to xgrow: {0}\n".format(repr(e)))
        return 1
        
    # Now create a temporary file, and run xgrow.
    
    temp = tempfile.NamedTemporaryFile(delete = False, mode='w')
    temp.write( tilestring )
    temp.close()
    
    ret = os.system( xgrowpath + " " + temp.name + " " + gseavgparam + " " + " ".join(args.xgrow_options))
    
    
    if ret != 0:
        sys.stderr.write("Xgrow failed, giving return code {0}.\n".format(ret)+
                         "Tile file has been kept at path {0}\n".format(temp.name))
        return 2
    else:
        os.remove(temp.name)
        
    

def print_help():
    print("""
The STXG Xgrow Wrapper
======================

Usage: xgrow-wrap [stxg-file] [xgrow-arguments]

xgrow-wrap is a wrapper for xgrow that uses stxg (structured xgrow, essentially
yaml) files instead of normal xgrow tile files. The first argument should be
the stxg file. All arguments after this are passed directly to xgrow without
modification, and, like xgrow files, will override arguments in the stxg file.

For a description of arguments to xgrow, use xgrow --

xgrow-wrap searches for xgrow in the current directory, then any directory
specified by the XGROW_DIR environment variable, and then your PATH.

Send questions or comments to Constantine Evans <cge@dna.caltech.edu>.

Last updated in 2013. This is version 0.0.1.
    """)

if __name__ == '__main__':
    main()
