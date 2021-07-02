from __future__ import annotations
import os
from typing import Any, Literal, Optional, Sequence, Tuple, overload
import tempfile
import subprocess
from subprocess import CompletedProcess
import re
from .tileset import TileSet


_XGROW_BINARY = os.path.join(os.path.dirname(__file__), '_xgrow')
_XGROW_BINARY = re.sub(r" ", r"\ ", _XGROW_BINARY)


def _process_outputs(outputs: dict[OutputOpts, Any]) -> None:
    from . import parseoutput
    for key in outputs.keys():
        if key == 'array':
            outputs[key] = parseoutput.load_array(outputs[key])
        elif key == 'trace':
            outputs[key] = parseoutput.load_trace(outputs[key])
        elif key == 'data':
            outputs[key] = parseoutput.load_data(outputs[key])
        else:
            raise ValueError("Output type {} is unknown.".format(key), key)


def run_raw(args: list[str], process_info=False) -> subprocess.CompletedProcess:
    """
    Dealing with the binary path, just run xgrow, with an arbitrary `argstring`.
    This is primarily designed for situations where (a) the user wants to
    manage everything themself, or (b) another function is handling parsing,
    and calls this to run xgrow with an argstring it constructs.

    Parameters
    ==========

    argstring: a raw argument string for xgrow.

    process_info: if True, pipes stdout and stderr, returning them in the
    CompletedProcess.  If False, then leaves them as is.

    Returns a subprocess CompletedProcess.
    """
    if process_info:
        return subprocess.run([_XGROW_BINARY] + args,
                              shell=False,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              encoding="utf-8")
    else:
        return subprocess.run([_XGROW_BINARY] + args,
                              shell=False)


OutputOpts = Literal['data', 'array', 'trace']


@overload
def run_old(tilestring: str, extraparams: dict = {},
            outputopts: OutputOpts = None,
            process_info: Literal[True] = True) -> Tuple[Any, CompletedProcess]:
    ...


@overload
def run_old(tilestring: str, extraparams: dict = {},
            outputopts: Sequence[OutputOpts] = None,
            process_info: Literal[True] = True) -> Tuple[dict[str, Any],
                                                         CompletedProcess]:
    ...


@overload
def run_old(tilestring: str, extraparams: dict = {},
            outputopts: None = None,
            process_info: Literal[True] = True) -> Tuple[None, CompletedProcess]:
    ...


@overload
def run_old(tilestring: str, extraparams: dict = {},
            outputopts: OutputOpts = None,
            process_info: Literal[False] = False) -> Any:
    ...


@overload
def run_old(tilestring: str, extraparams: dict = {},
            outputopts: Sequence[OutputOpts] = None,
            process_info: Literal[False] = False) -> dict[str, Any]:
    ...


@overload
def run_old(tilestring: str, extraparams: dict = {},
            outputopts: None = None,
            process_info: Literal[False] = False) -> None:
    ...


@overload
def run_old(tilestring: str, extraparams: dict = {},
            outputopts: OutputOpts | Sequence[OutputOpts] | None = None,
            process_info: bool = False) -> Any:
    ...


def run_old(tilestring: str, extraparams: dict = {},
            outputopts: OutputOpts | Sequence[OutputOpts] | None = None,
            process_info: bool = False) -> Any:
    """
    Given an old xgrow tileset definition (as a string), a dict of parameters,
    and possible output options, run xgrow, handling file creation.

    Parameters
    ==========

    tilestring: a string containing an old xgrow tileset definition, which
    will be put in a temporary file and passed directly to xgrow.

    extraparams: a dict of key, value pairs of options for xgrow itself.  There
    is no special handling of these.  Output options should instead go in
    outputopts, which will handle output file creation and reading.

    outputopts: either a string or list, specifying output options.  If a
    string, one of 'data', 'array', or 'trace', corresponding to 'datafile',
    'arrayfile' and 'tracefile', respectively.  If a list of multiple, then
    do those.  These will manage the output, and return the data raw from xgrow
    as strings.

    process_info: if True, add a 'process_info' key to the return dictionary,
    with the subprocess.CompletedProcess instance for the xgrow run.  This
    includes the stdout and stderr, and so on.
    """

    # Create necessary temp files:
    tileset_file = tempfile.NamedTemporaryFile(
        delete=False, mode='w', newline='\n')

    if outputopts is None:
        outputopts = []
    elif isinstance(outputopts, str):
        outputopts = [outputopts]  # type: ignore
    output_files = {output_type: tempfile.mktemp(prefix=f"tmp_xgrow_{output_type}")
                    for output_type in outputopts}  # type: ignore

    tileset_file.write(tilestring)
    tileset_file.close()

    args = [tileset_file.name]
    for param, val in extraparams.items():
        if val is True:
            args.append("{}".format(param))
        else:
            args.append("{}={}".format(param, val))

    args += ["{}file={}".format(output_type, output_file)
             for output_type, output_file
             in output_files.items()]

    ret = run_raw(args, process_info=process_info)
    if ret.returncode != 0:
        raise Exception(
            "Xgrow failed with return code {}.".format(ret.returncode),
            ret,
            args,
            tileset_file.name,
            {t: f for t, f in output_files.values()})

    os.unlink(tileset_file.name)

    output = {}
    for output_type, output_file in output_files.items():
        with open(output_file, 'r') as output_file_reopened:
            output[output_type] = output_file_reopened.read()
        os.unlink(output_file)

    _process_outputs(output)

    if not output:
        output = None
    if len(outputopts) == 1:    # type: ignore
        output = output[outputopts[0]]    # type: ignore

    if not process_info:
        return output
    else:
        return output, ret


@overload
def run(tileset: dict | TileSet,
        extraparams: dict,
        outputopts: None = None,
        ui: Optional[bool] = None,
        process_info: Literal[True] = True) -> Tuple[None, subprocess.CompletedProcess]:
    ...


@overload
def run(tileset: dict | TileSet,
        extraparams: dict,
        outputopts: OutputOpts = None,
        ui: Optional[bool] = None,
        process_info: Literal[True] = True) -> Tuple[Any, subprocess.CompletedProcess]:
    ...


@overload
def run(tileset: dict | TileSet,
        extraparams: dict,
        outputopts: Sequence[OutputOpts] = None,
        ui: Optional[bool] = None,
        process_info: Literal[True] = True) -> Tuple[dict, subprocess.CompletedProcess]:
    ...


@overload
def run(tileset: dict | TileSet, extraparams: dict = {},
        outputopts: None = None,
        ui: Optional[bool] = None, process_info: Literal[False] = False) -> None:
    ...


@overload
def run(tileset: dict | TileSet, extraparams: dict = {},
        outputopts: OutputOpts = None,
        ui: Optional[bool] = None, process_info: Literal[False] = False) -> Any:
    ...


@overload
def run(tileset: dict | TileSet, extraparams: dict = {},
        outputopts: Sequence[OutputOpts] = None,
        ui: Optional[bool] = None,
        process_info: Literal[False] = False) -> dict[str, Any]:
    ...


def run(tileset: dict | TileSet, extraparams: dict = {},
        outputopts: OutputOpts | Sequence[OutputOpts] | None = None,
        ui: Optional[bool] = None, process_info: bool = False) -> Any:
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

    if not isinstance(tileset, TileSet):
        tileset = TileSet.from_dict(tileset)

    tileset.xgrowargs.merge(extraparams)

    if ui is not None:
        tileset.xgrowargs.window = ui

    return run_old(tileset.to_xgrow(), outputopts=outputopts, process_info=process_info)
