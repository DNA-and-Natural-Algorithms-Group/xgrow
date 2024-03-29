from __future__ import annotations
import os
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Literal,
    Sequence,
    SupportsIndex,
    Tuple,
    cast,
    overload,
)
import tempfile
import subprocess
from subprocess import CompletedProcess
import re
from typing_extensions import TypeAlias

from xgrow.parseoutput import XgrowOutput

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd
    PossibleXgrowOutputs: TypeAlias = "XgrowOutput | pd.Series | pd.DataFrame"
    from .tileset import XgrowArgs

OutputOpts = Literal["data", "array", "trace"]

_DEFAULT_SIZE = 64

_XGROW_BINARY = re.sub(r" ", r"\ ", os.path.join(os.path.dirname(__file__), "_xgrow"))


def _process_outputs(output_files: Dict[OutputOpts, str]) -> Dict[OutputOpts, Any]:
    if not output_files:
        return {}

    from . import parseoutput

    outputs: Dict[OutputOpts, Any] = {}
    for key in output_files.keys():
        if key == "array":
            outputs[key] = parseoutput.load_array_file(output_files[key])
        elif key == "trace":
            outputs[key] = parseoutput.load_trace_file(output_files[key])
        elif key == "data":
            outputs[key] = parseoutput.load_data_file(output_files[key])
        else:
            raise ValueError(f"Output type {key} is unknown.", key)
    return outputs


def run_raw(
    args: list[str], process_info: bool = False, quiet: bool = False
) -> subprocess.CompletedProcess[str]:
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

    quiet: if True and process_info is False, then ignores stdout and stderr.

    Returns a subprocess CompletedProcess.
    """
    if process_info:
        return subprocess.run(
            [_XGROW_BINARY] + args,
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding="utf-8",
            check=True,
        )
    elif quiet:
        return subprocess.run(
            [_XGROW_BINARY] + args,
            shell=False,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            encoding="utf-8",
            check=True,
        )
    return subprocess.run(
        [_XGROW_BINARY] + args, shell=False, encoding="utf-8", check=True
    )


@overload
def run_old(
    tilestring: str,
    extraparams: Dict[str, Any] | None = None,
    *,
    outputopts: OutputOpts,
    process_info: Literal[True],
) -> Tuple[XgrowArgs, CompletedProcess[str]]:
    ...


# @overload
# def run_old(
#     tilestring: str,
#     extraparams: Dict[str, Any] | None = None,
#     *,
#     outputopts: Sequence[OutputOpts],
#     process_info: Literal[True],
# ) -> Tuple[Dict[str, PossibleXgrowOutputs], CompletedProcess[str]]:
#     ...


@overload
def run_old(
    tilestring: str,
    extraparams: Dict[str, Any] | None = None,
    outputopts: None = None,
    *,
    process_info: Literal[True],
) -> Tuple[None, CompletedProcess[str]]:
    ...


@overload
def run_old(
    tilestring: str,
    extraparams: Dict[str, Any] | None = None,
    *,
    outputopts: OutputOpts,
    process_info: Literal[False] = False,
) -> PossibleXgrowOutputs:
    ...


# @overload
# def run_old(
#     tilestring: str,
#     extraparams: Dict[str, Any] | None = None,
#     *,
#     outputopts: Sequence[OutputOpts],
#     process_info: Literal[False] = False,
# ) -> Dict[str, Any]:
#     ...


@overload
def run_old(
    tilestring: str,
    extraparams: Dict[str, Any] | None = None,
    outputopts: None = None,
    process_info: Literal[False] = False,
    quiet: bool = False
) -> None:
    ...


@overload
def run_old(
    tilestring: str,
    extraparams: Dict[str, Any] | None = None,
    outputopts: OutputOpts | Sequence[OutputOpts] | None = None,
    process_info: bool = False,
    quiet: bool = False
) -> PossibleXgrowOutputs | Dict[str, PossibleXgrowOutputs]:
    ...


def run_old(
    tilestring: str,
    extraparams: Dict[str, Any] | None = None,
    outputopts: OutputOpts | Sequence[OutputOpts] | None = None,
    process_info: bool = False,
    quiet: bool = False
) -> PossibleXgrowOutputs | Dict[str, PossibleXgrowOutputs]:
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

    if extraparams is None:
        extraparams = {}

    # Create necessary temp files:
    with tempfile.NamedTemporaryFile(
        delete=False, mode="w", newline="\n"
    ) as tileset_file:

        if isinstance(outputopts, str):
            outputopts = [outputopts]  # type: ignore
        if outputopts is None:
            outputopts = []
        output_files: Dict[OutputOpts, str] = {
            output_type: tempfile.mktemp(prefix=f"tmp_xgrow_{output_type}")
            for output_type in outputopts
        }  # type: ignore

        tileset_file.write(tilestring)

    args = [tileset_file.name]
    for param, val in extraparams.items():
        if val is True:
            args.append(f"{param}")
        else:
            args.append(f"{param}={val}")

    args += [
        f"{output_type}file={output_file}"
        for output_type, output_file in output_files.items()
    ]

    ret = run_raw(args, process_info=process_info, quiet=quiet)
    if ret.returncode != 0:
        raise Exception(
            f"Xgrow failed with return code {ret.returncode}.",
            ret,
            args,
            tileset_file.name,
            list(output_files.values()),
        )

    os.unlink(tileset_file.name)

    output = _process_outputs(output_files)
    for output_file in output_files.values():
        os.unlink(output_file)

    if len(outputopts) == 1:  # type: ignore
        output = output[outputopts[0]]  # type: ignore

    if not process_info:
        return output if output is not None else None
    return (output if output is not None else None), ret


# @overload
# def run(
#     tileset: Dict[str, Any] | TileSet,
#     extraparams: Dict[str, Any],
#     *,
#     outputopts: Sequence[OutputOpts],
#     process_info: Literal[True],
#     **kwargs: Dict[str, Any],
# ) -> Tuple[Dict[str, PossibleXgrowOutputs], subprocess.CompletedProcess[str]]:
#     ...

# @overload
# def run(
#     tileset: Dict[str, Any] | TileSet,
#     extraparams: Dict[str, Any] | None = None,
#     *,
#     outputopts: Sequence[OutputOpts],
#     process_info: Literal[False] = False,
#     **kwargs: Dict[str, Any],
# ) -> Dict[str, PossibleXgrowOutputs]:
#     ...
