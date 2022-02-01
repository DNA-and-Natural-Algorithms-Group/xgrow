from pathlib import Path
import pytest
import glob
from xgrow.parseoutput import XgrowOutput

from xgrow.tileset import TileSet

import numpy as np

import xgrow


def test_untiltiles():
    ts = TileSet.from_yaml("examples/barish-perfect.yaml")

    out: XgrowOutput = xgrow.run(
        ts,
        outputopts="array",
        T=12,
        untiltiles=["ZG", "ZE"],
        size=128,
        emax=1000,
        window=False,
        importfile="examples/tallrect.seed",
    )
    _, tilenums = ts.to_xgrow(return_tilenums=True)

    assert np.any(out.tiles == tilenums["ZG"])
    assert np.any(out.tiles == tilenums["ZE"])
    assert out.data.tiles == 652
