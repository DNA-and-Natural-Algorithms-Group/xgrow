from pathlib import Path
import pytest
import glob

from xgrow.tileset import TileSet

import xgrow


@pytest.mark.parametrize("file", glob.glob("examples/*.yaml"))
def test_load_examples(file, tmp_path: Path):
    ts = TileSet.from_yaml(file)

    assert ts == TileSet.from_dict(ts.to_dict())


@pytest.mark.parametrize("file", glob.glob("examples/*.yaml"))
def test_run_examples(file):
    ts = TileSet.from_yaml(file)

    if ts.xgrowargs.importfile is not None:
        return
        # ts.xgrowargs.importfile == str("examples" / Path(ts.xgrowargs.importfile))

    out = xgrow.run(ts, outputopts="array", pause=False, tmax=1000, window=False)

    out.show_array(ts)


@pytest.mark.parametrize("file", glob.glob("examples/xgrow-format/sierp*"))
def test_old_examples(file):
    xgrow.run_old(
        open(file).read(),
        extraparams={"pause": "False", "tmax": "1000", "window": "False"},
    )
