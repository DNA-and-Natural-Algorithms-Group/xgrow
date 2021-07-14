from dataclasses import dataclass
import re
from typing import Any, Optional, Union
import numpy as np
from io import BytesIO
import pandas as pd


@dataclass
class XgrowOutput:
    data: Optional[pd.Series] = None
    tiles: Optional[np.ndarray] = None


def load_array_file(fn: str, onlytiles: bool = False) -> XgrowOutput:

    s = open(fn, "rb").read()
    s = re.sub(
        rb"(\[|\])",
        b"",  # remove [ and ]
        re.sub(rb"; \.\.\.", b"", s),  # remove ; ... at line ends
    )

    tiles = np.genfromtxt(
        BytesIO(s),
        skip_header=4,  # remove first lines
        skip_footer=1,  # remove end junk
        dtype="uint",
    )
    if onlytiles:
        return tiles
    data = pd.Series(
        np.genfromtxt(BytesIO((s.split(b"\n")[2])))[:-1],
        index=[
            "gmc",
            "gse",
            "k",
            "time",
            "tiles",
            "mismatches",
            "events",
            "perimeter",
            "g",
            "dgbonds",
        ],
    )
    return XgrowOutput(data, tiles)


def load_trace_file(fn: str) -> pd.DataFrame:
    data = pd.DataFrame(
        np.genfromtxt(fn),
        columns=[
            "gmc",
            "gse",
            "k",
            "time",
            "tiles",
            "mismatches",
            "events",
            "perimeter",
            "g",
            "dgbonds",
        ],
    )
    return data


def load_data_file(fn: str) -> pd.Series:
    data = pd.Series(
        np.genfromtxt(fn),
        index=[
            "gmc",
            "gse",
            "k",
            "time",
            "tiles",
            "mismatches",
            "events",
            "perimeter",
            "g",
            "dgbonds",
        ],
    )
    return data


def show_array(a: np.ndarray, ts: dict[str, Any], **kwargs: dict[str, Any]):
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from .xcolors import xcolors  # FIXME:  put in here!

    mcolors = {n: tuple(z / 255.0 for z in eval(x[3:])) for n, x in xcolors.items()}
    cmap = colors.ListedColormap(
        ["black"] + [mcolors[x["color"]] for x in ts["tiles"]]  # type: ignore
    )
    try:
        plt.imshow(a["tiles"], cmap=cmap, vmin=0, vmax=len(ts["tiles"]), **kwargs)
    except KeyError:
        plt.imshow(a, cmap=cmap, vmin=0, vmax=len(ts["tiles"]), **kwargs)
