from dataclasses import dataclass
import re
from typing import TYPE_CHECKING, Any, Optional, Union, Dict
import numpy as np
from io import BytesIO, StringIO
import pandas as pd

if TYPE_CHECKING:
    from .tileset import TileSet


@dataclass
class XgrowOutput:
    data: Optional[pd.Series] = None
    tiles: Optional[np.ndarray] = None

    def show_array(self, ts: "TileSet", **kwargs: Dict[str, Any]):
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        from .xcolors import mcolors  # FIXME:  put in here!

        out = StringIO()

        _, tilenums = ts.to_xgrow(out, return_tilenums=True)
        back = {v: k for k, v in tilenums.items()}
        vmax = max(back.keys())

        tilecolors = {t.name: t.color for t in ts.tiles}

        cmap = colors.ListedColormap(
            ["black"] + [mcolors.get(tilecolors.get(t, None), "gray") for i, t in back.items()]  # type: ignore
        )
        return plt.imshow(self.tiles, cmap=cmap, vmin=0, vmax=vmax, **kwargs)


def load_array_file(fn: str) -> XgrowOutput:

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
