from __future__ import annotations
import os
import tempfile
from .xgrow_subprocess import _DEFAULT_SIZE, run_old, Literal, OutputOpts, XgrowOutput
from .tileset import TileSet
from typing import TYPE_CHECKING, Sequence, overload, Tuple, Dict, Any
from typing_extensions import TypeAlias
import subprocess

if TYPE_CHECKING:
    import pandas as pd
    PossibleXgrowOutputs: TypeAlias = "XgrowOutput | pd.Series | pd.DataFrame"

@overload
def run(
    tileset: Dict[str, Any] | TileSet,
    extraparams: Dict[str, Any] | None = None,
    outputopts: None = None,
    *,
    process_info: Literal[True],
    **kwargs: Dict[str, Any],
) -> Tuple[None, subprocess.CompletedProcess[str]]:
    ...


@overload
def run(
    tileset: Dict[str, Any] | TileSet,
    extraparams: Dict[str, Any] | None = None,
    *,
    outputopts: OutputOpts,
    process_info: Literal[True],
    **kwargs: Dict[str, Any],
) -> Tuple[PossibleXgrowOutputs, subprocess.CompletedProcess[str]]:
    ...


@overload
def run(
    tileset: Dict[str, Any] | TileSet,
    extraparams: Dict[str, Any] | None = None,
    outputopts: None = None,
    process_info: Literal[False] = False,
    **kwargs: Dict[str, Any],
) -> None:
    ...


@overload
def run(
    tileset: Dict[str, Any] | TileSet,
    extraparams: Dict[str, Any] | None = None,
    *,
    outputopts: Literal["data"],
    process_info: Literal[False] = False,
    **kwargs: Dict[str, Any],
) -> pd.Series:
    ...


@overload
def run(
    tileset: Dict[str, Any] | TileSet,
    extraparams: Dict[str, Any] | None = None,
    *,
    outputopts: Literal["array"],
    process_info: Literal[False] = False,
    **kwargs: Dict[str, Any],
) -> XgrowOutput:
    ...


@overload
def run(
    tileset: Dict[str, Any] | TileSet,
    extraparams: Dict[str, Any] | None = None,
    *,
    outputopts: Literal["trace"],
    process_info: Literal[False] = False,
    **kwargs: Dict[str, Any],
) -> pd.DataFrame:
    ...


@overload
def run(
    tileset: Dict[str, Any] | TileSet,
    extraparams: Dict[str, Any] | None = None,
    *,
    outputopts: OutputOpts,
    process_info: Literal[False] = False,
    **kwargs: Dict[str, Any],
) -> PossibleXgrowOutputs:
    ...


def run(
    tileset: Dict[str, Any] | TileSet,
    extraparams: Dict[str, Any] | None = None,
    outputopts: OutputOpts | Sequence[OutputOpts] | None = None,
    process_info: bool = False,
    **kwargs: Dict[str, Any],
) -> Any:
    """Given a tileset (class or dict), and a dictionary of extra parameters,
    run xgrow, potentially with particular managed output options.  This
    replaces the xgrow-running code in xgrow_wrap and in xgrow_parallel.

    Parameters
    ==========

    tileset: a tileset to run, as an stxg-format dictionary.

    extraparams: a dictionary of options to xgrow, as in the xrgrowargs
    key in the tileset.

    outputopts: either a string or list, specifying output options.  If a
    string, one of 'final', 'array', or 'trace', corresponding to 'datafile',
    'arrayfile' and 'tracefile', respectively.  If a list of multiple, then
    do those.  These will manage the output, and return the data in usable
    form.

    ui: if True, then show the xgrow ui.  If false, then suppress the ui.  If
    suppressing the ui, the tileset or extraparams should include a terminating
    option, like smax or emax; otherwise, this function may never halt.

    Output
    ======

    If no outputopts, then nothing is returned.

    If outputopt is a string, the function returns the requested output.

    If outputopt is a list of strings, then the function returns a dict of
    outputs with keys of the same name.
    """

    if outputopts is None:
        outputopts = tuple()

    if extraparams is None:
        extraparams = {}

    if not isinstance(tileset, TileSet):
        tileset = TileSet.from_dict(tileset)

    return tileset.run(extraparams=extraparams, outputopts=outputopts, process_info=process_info, **kwargs) # type: ignore
