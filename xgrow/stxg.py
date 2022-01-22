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
    ("block", int),
    ("size", int),
    ("rand", str),
    ("k", float),
    ("Gmc", float),
    ("Gse", float),
    ("Gas", float),
    ("Gam", float),
    ("Gae", float),
    ("Gfc", float),
    ("T", float),
    ("blast_rate_alpha", float),
    ("blast_rate_beta", float),
    ("blast_rate_gamma", float),
    ("seed", str),
    ("update_rate", int),
    ("tracefile", str),
    ("untiltiles", str),
    ("tmax", float),
    ("emax", int),
    ("smax", int),
    ("smin", int),
    ("untiltilescount", str),
    ("clean_cycles", int),
    ("error_radius", float),
    ("datafile", str),
    ("arrayfile", str),
    ("exportfile", str),
    ("importfile", str),
    ("min_strength", float),
    ("window", bool),
    ("fission", str),
    ("pause", bool),
    ("wander", bool),
    ("periodic", bool),
]
keyopts = [x[0] for x in keyvaloptions]

# Options that need special handling.
specialoptions = ["stoic"]  # FIXME


def from_xgrow(xgst):
    """
    Given a string of an xgrow tile file, convert it to stxg.

    Note that the conversion may require some work afterward.
    """

    ntre = re.compile(r"^\s*num tile types=(\d+)", re.MULTILINE)
    nbre = re.compile(r"^\s*num binding types=(\d+)", re.MULTILINE)
    bsre = re.compile(r"^\s*binding strengths\s*=\s*{([^}]+)}", re.MULTILINE)
    tere = re.compile(
        r"^\s*tile edges\s*=[\s\n]*{((?:[^{}]*{[^{}]+}[^{}]*)+)}", re.MULTILINE
    )
    bnre = re.compile(r"^\s*binding type names\s*=\s*{([^}]+)}", re.MULTILINE)
    glre = re.compile(r"^\s*g\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*=\s*(\S+)", re.MULTILINE)

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
        r"{\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*}(?:\[(\S+)\])?(?:\((\S+)\))?", tstring
    )

    tiles = []
    for tv in tilevals:
        tile = {}
        tile["edges"] = []
        for e in tv[0:4]:
            try:
                e = int(e)
            except:
                pass
            tile["edges"].append(e)
        if tv[4] != "":
            tile["stoic"] = float(tv[4])
        if tv[5] != "":
            tile["color"] = tv[5]
        tiles.append(tile)

    if bnmatch:
        assert len(bnames) == len(bstrengths)
        bonds = [
            {"name": name, "strength": float(strength)}
            for name, strength in zip(bnames, bstrengths)
        ]
    else:
        bonds = [{"strength": float(strength)} for strength in bstrengths]

    assert len(bonds) == nbind
    assert len(tiles) == ntiles

    ts = {"tiles": tiles, "bonds": bonds}

    glues = {}
    m = glre.findall(xgs)
    if m:
        for x1, x2, s in m:
            if bnmatch:
                glues[(bnames[int(x1) - 1], bnames[int(x2) - 1])] = float(s)
            else:
                glues[int(x1), int(x2)] = float(s)
        ts["glues"] = glues

    xgrowargs = {}
    for opt, opttype in keyvaloptions:
        m = re.findall(r"^\s*" + opt + r"\s*=\s*(\S+)", xgs, re.MULTILINE)
        if len(m) > 0:
            xgrowargs[opt] = opttype(m[-1])

        # fission
    m = re.findall(r"^\s*(chunk_fission|fission|no_fission)", xgs, re.MULTILINE)
    if len(m) > 0:
        if m[-1] == "chunk_fission":
            xgrowargs["fission"] = "chunk"

        elif m[-1] == "fission":
            xgrowargs["fission"] = "on"
        else:
            xgrowargs["fission"] = "off"

    # doubletile
    m = re.findall(r"^\s*doubletile\s*=\s*(\d+)\s*,\s*(\d+)", xgs, re.MULTILINE)
    if m:
        xgrowargs["doubletiles"] = [[int(x1), int(x2)] for x1, x2 in m]

    if re.search(r"^\s*-nw", xgs, re.MULTILINE):
        xgrowargs["window"] = False

    ts["xgrowargs"] = xgrowargs
    return ts
