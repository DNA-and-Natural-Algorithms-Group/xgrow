# Structured Xgrow libary for reading, writing and converting Structured Xgrow
# (STXG) files.  If you're just manipulating stxg files, you probably don't
# need this (just use yaml), unless you want to use the dump function, which
# produces somewhat cosmetically better output. The library is primarily useful
# for converting to and from xgrow tile files.
#
# If you just want to run xgrow with stxg files, consider using xgrow-wrap.

import yaml
import re
from io import StringIO  # Fixme: stupid Python 2/3 hack.
import datetime

version = 0.5

# Option names and their argument types.
keyvaloptions = [
    ('block', int), ('size', int), ('rand', str), ('k', float), ('Gmc', float),
    ('Gse', float), ('Gas', float), ('Gam', float), ('Gae', float),
    ('Gfc', float), ('T', float), ('blast_rate_alpha', float),
    ('blast_rate_beta', float), ('blast_rate_gamma', float), ('seed', str),
    ('update_rate', int), ('tracefile', str), ('untiltiles', str),
    ('tmax', float), ('emax', int), ('smax', int), ('smin', int),
    ('untiltilescount', str), ('clean_cycles', int), ('error_radius', float),
    ('datafile', str), ('arrayfile', str), ('exportfile', str),
    ('importfile', str), ('min_strength', float)
]
keyopts = [x[0] for x in keyvaloptions]

# Options that need special handling.
specialoptions = ['stoic', 'nowindow', 'doubletile', 'vdoubletile', 'fission']

# Options that take true/fale values.
truefalseoptions = ['pause', 'wander', 'periodic']


# Classes for controlling flow/block style of various yaml elements. These are
# used in dump()
class blockseq(dict):
    pass


def blockseq_rep(dumper, data):
    return dumper.represent_mapping(
        u'tag:yaml.org,2002:seq', data, flow_style=False)


class blockmap(dict):
    pass


def blockmap_rep(dumper, data):
    return dumper.represent_mapping(
        u'tag:yaml.org,2002:map', data, flow_style=False)


class flowmap(dict):
    pass


def flowmap_rep(dumper, data):
    return dumper.represent_mapping(
        u'tag:yaml.org,2002:map', data, flow_style=True)


yaml.add_representer(blockseq, blockseq_rep)
yaml.add_representer(flowmap, flowmap_rep)
yaml.add_representer(blockmap, blockmap_rep)


def dump(tileset, *xargs, **pargs):
    """stxg.dump( tileset, ... ) is a wrapper for yaml.dump that slightly
    prettifies the output taking into account what tilesets usually look
    like. It passes all arguments onto yaml.dump.
    """

    # Make tile and bond sections flow-style
    tileset['tiles'] = [flowmap(x) for x in tileset['tiles']]
    tileset['bonds'] = [flowmap(x) for x in tileset['bonds']]

    # If xgrowargs is there, make it block-style
    if 'xgrowargs' in tileset.keys():
        tileset['xgrowargs'] = blockmap(tileset['xgrowargs'])

    return yaml.dump(tileset, *xargs, **pargs)


def load(stream, *xargs, **pargs):
    """
    stxg.load( stream, ... ) is a wrapper for yaml.load that currently does
    nothing else but call yaml.load
    """
    return yaml.load(stream, *xargs, **pargs)


def to_xgrow(stxg, stream=None):
    """
    Given an stxg structure/dict, create an xgrow tile file, and return it
    as a string. Alternatively, if the stream= parameter is given a stream,
    write the file to it instead.
    """
    # If there's no stream to write to, create a StringIO to write to:
    if stream:
        xgrowf = stream
    else:
        xgrowf = StringIO()

    # Add some comments:

    for n, tile in enumerate(stxg['tiles']):
        tileedges = [str(x) for x in tile['edges']]
        for e in tileedges:
            if e != '0' and e not in [x['name'] for x in stxg['bonds']]:
                stxg['bonds'].append({'name': e, 'strength': 1})
    # Start writing the tileset
    try:
        xgrowf.write("num tile types=%d\n" % len(stxg['tiles']))
        xgrowf.write("num binding types=%d\n" % len(stxg['bonds']))
    except KeyError as e:
        raise ValueError("Missing section {0}.".format(e.message))

    bondhasname = ['name' in x.keys() for x in stxg['bonds']]

    if all(bondhasname):
        xgrowf.write("binding type names={ %s }\n" % " ".join(
            [x['name'] for x in stxg['bonds']]))
    elif any(bondhasname):
        raise ValueError("Some bonds have names, but not all. " +
                         "All bonds must be either named or unnamed.")

    xgrowf.write("tile edges={\n")

    for n, tile in enumerate(stxg['tiles']):
        xgrowf.write("{ %s }" % " ".join([str(x) for x in tile['edges']]))
        if 'stoic' in tile:
            xgrowf.write("[%g]" % tile['stoic'])
        if 'color' in tile:
            xgrowf.write("(%s)" % tile['color'])
        xgrowf.write("   % ")
        if 'name' in tile:
            xgrowf.write(tile['name'])
        xgrowf.write(" (tile #%d)\n" % n)

    xgrowf.write("}\n")

    xgrowf.write("binding strengths={ %s }\n" % " ".join(
        [str(x['strength']) for x in stxg['bonds']]))

    if 'glues' in stxg:
        for x1, x2, g in [tuple(x) for x in stxg['glues']]:
            # dereference
            if type(x1) != int:
                try:
                    x1 = [n for n, z in enumerate(stxg['bonds'])
                          if z['name'] == x1][0] + 1
                except IndexError:
                    print(x2)
                    raise ValueError(x2)
            if type(x2) != int:
                try:
                    x2 = [n for n, z in enumerate(stxg['bonds'])
                          if z['name'] == x2][0] + 1
                except IndexError:
                    print(x2)
                    raise ValueError(x2)
            xgrowf.write("g(%d,%d)=%g\n" % (x1, x2, float(g)))

    if 'xgrowargs' in stxg.keys():
        for key, val in stxg['xgrowargs'].items():
            if key == 'doubletiles':
                for x1, x2 in val:
                    # dereference
                    if type(x1) != int:
                        try:
                            x1 = [n for n, z in enumerate(stxg['tiles'])
                                  if z['name'] == x1][0] + 1
                        except IndexError:
                            print(x1)
                            raise ValueError(x1)
                    if type(x2) != int:
                        try:
                            x2 = [n for n, z in enumerate(stxg['tiles'])
                                  if z['name'] == x2][0] + 1
                        except IndexError:
                            print(x2)
                            raise ValueError(x2)
                    xgrowf.write("doubletile=%d,%d\n" % (x1, x2))
            elif key == 'vdoubletiles':
                for x1, x2 in val:
                    # dereference
                    if type(x1) != int:
                        try:
                            x1 = [n for n, z in enumerate(stxg['tiles'])
                                  if z['name'] == x1][0] + 1
                        except IndexError:
                            print(x1)
                            raise ValueError(x1)
                    if type(x2) != int:
                        try:
                            x2 = [n for n, z in enumerate(stxg['tiles'])
                                  if z['name'] == x2][0] + 1
                        except IndexError:
                            print(x2)
                            raise ValueError(x2)
                    xgrowf.write("vdoubletile=%d,%d\n" % (x1, x2))
            else:
                xgrowf.write("%s=%s\n" % (key, str(val)))

    xgrowf.write("% Tileset created by stxg.py version {0} on {1}\n".format(
        version, str(datetime.date.today())))

    if stream:
        return None
    else:
        return xgrowf.getvalue()


def from_xgrow(xgst):
    """
    Given a string of an xgrow tile file, convert it to stxg.
    
    Note that the conversion may require some work afterward.
    """

    ntre = re.compile(r"^\s*num tile types=(\d+)", re.MULTILINE)
    nbre = re.compile(r"^\s*num binding types=(\d+)", re.MULTILINE)
    bsre = re.compile(r"^\s*binding strengths\s*=\s*{([^}]+)}", re.MULTILINE)
    tere = re.compile(
        r"^\s*tile edges\s*=[\s\n]*{((?:[^{}]*{[^{}]+}[^{}]*)+)}",
        re.MULTILINE)
    bnre = re.compile(r"^\s*binding type names\s*=\s*{([^}]+)}", re.MULTILINE)
    glre = re.compile(r"^\s*g\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*=\s*(\S+)",
                      re.MULTILINE)

    xgs = re.sub(r"%[^\n]\n", "\n", xgst)

    m = ntre.search(xgs)
    if not m:
        #  warning.warn("Expected number of tile types.")
        pass
    ntiles = int(m.group(1))

    m = nbre.search(xgs)
    if not m:
        #  warning.warn("Expected number of binding types.")
        pass
    nbind = int(m.group(1))

    bsstring = bsre.search(xgs).group(1)
    tstring = tere.search(xgs).group(1)
    bnmatch = bnre.search(xgs)
    if bnmatch:
        bnstring = bnmatch.group(1)
        bnames = bnstring.split()

    bstrengths = bsstring.split()
    tilevals = re.findall(
        r"{\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*}(?:\[(\S+)\])?(?:\((\S+)\))?",
        tstring)

    tiles = []
    for tv in tilevals:
        tile = {}
        tile['edges'] = []
        for e in tv[0:4]:
            try:
                e = int(e)
            except:
                pass
            tile['edges'].append(e)
        if tv[4] != '':
            tile['stoic'] = float(tv[4])
        if tv[5] != '':
            tile['color'] = tv[5]
        tiles.append(tile)

    if bnmatch:
        assert len(bnames) == len(bstrengths)
        bonds = [{'name': name,
                  'strength': float(strength)}
                 for name, strength in zip(bnames, bstrengths)]
    else:
        bonds = [{'strength': float(strength)} for strength in bstrengths]

    assert len(bonds) == nbind
    assert len(tiles) == ntiles

    ts = {'tiles': tiles, 'bonds': bonds}

    glues = {}
    m = glre.findall(xgs)
    if m:
        for x1, x2, s in m:
            if bnmatch:
                glues[(bnames[int(x1) - 1], bnames[int(x2) - 1])] = float(s)
            else:
                glues[int(x1), int(x2)] = float(s)
        ts['glues'] = glues

    xgrowargs = {}
    for opt, opttype in keyvaloptions:
        m = re.findall(r"^\s*" + opt + r"\s*=\s*(\S+)", xgs, re.MULTILINE)
        if len(m) > 0:
            xgrowargs[opt] = opttype(m[-1])

        # fission
    m = re.findall(r"^\s*(chunk_fission|fission|no_fission)", xgs,
                   re.MULTILINE)
    if len(m) > 0:
        if m[-1] == 'chunk_fission':
            xgrowargs['fission'] = 'chunk'

        elif m[-1] == 'fission':
            xgrowargs['fission'] = 'on'
        else:
            xgrowargs['fission'] = 'off'

    # doubletile
    m = re.findall(r"^\s*doubletile\s*=\s*(\d+)\s*,\s*(\d+)", xgs,
                   re.MULTILINE)
    if m:
        xgrowargs['doubletiles'] = [[int(x1), int(x2)] for x1, x2 in m]

    if re.search(r"^\s*-nw", xgs, re.MULTILINE):
        xgrowargs['window'] = False

    ts['xgrowargs'] = xgrowargs
    return ts


def revcomp(seq):
    import stickydesign as sd
    seq = seq.lower()
    return "".join(reversed([sd.wc[nt] for nt in seq]))
