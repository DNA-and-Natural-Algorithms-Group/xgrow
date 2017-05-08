import sys
import os
import pkg_resources
import tempfile
import subprocess
import copy
from . import stxg
from . import parseoutput

_XGROW_BINARY = pkg_resources.resource_filename(__name__,'_xgrow')

def _process_outputs( outputs ):
    for key in outputs.keys():
        if key == 'array':
            outputs[key] = parseoutput.load_array(outputs[key])
        elif key == 'trace':
            outputs[key] = parseoutput.load_trace(outputs[key])
        elif key == 'data':
            outputs[key] = parseoutput.load_data(outputs[key])
            

def run_raw( argstring: str, process_info=False ):
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
        return subprocess.run(_XGROW_BINARY + " " + argstring,
                              shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              encoding="utf-8")
    else:
        return subprocess.run(_XGROW_BINARY + " " + argstring,
                              shell=True)
        

def run_old( tilestring: str, extraparams: dict, outputopts=None, process_info=False ):
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
    string, one of 'final', 'array', or 'trace', corresponding to 'datafile',
    'arrayfile' and 'tracefile', respectively.  If a list of multiple, then
    do those.  These will manage the output, and return the data raw from xgrow
    as strings.

    process_info: if True, add a 'process_info' key to the return dictionary,
    with the subprocess.CompletedProcess instance for the xgrow run.  This
    includes the stdout and stderr, and so on.
    """

    # Create necessary temp files:
    tileset_file = tempfile.NamedTemporaryFile(delete = False, mode='w')

    if outputopts is None:
        outputopts = []
    elif isinstance(outputopts, str):
        outputopts = [outputopts]
    output_files = { output_type: tempfile.NamedTemporaryFile(delete=False)
                     for output_type in outputopts }

    tileset_file.write(tilestring)
    tileset_file.close()

    # Close all output files.  We'll reopen them later.  This is important,
    # apparently, for Windows compatibility.
    for output_type, output_file in output_files.items():
        output_files[output_type]=output_file.name
        output_file.close()
    
    paramstring = ""
    for param, val in extraparams.items():
        if val is True:
            paramstring += " {}".format(param)
        else:
            paramstring += " {}={}".format(param,val)

    outputstring = " ".join( "{}file={}".format(output_type,output_file)
                             for output_type, output_file
                             in output_files.items() )
            
    argstring = tileset_file.name + " " + \
                paramstring + " " + \
                outputstring

    ret = run_raw( argstring, process_info=process_info )
    if ret.returncode != 0:
        raise Exception(
            "Xgrow failed with return code {}.".format(ret.returncode),
            ret,
            argstring,
            tileset_file.name,
            { t: f.name for t,f in outputstring.values() } )

    os.unlink(tileset_file.name)
    
    output = {}
    for output_type, output_file in output_files.items():
        with open(output_file,'r') as output_file_reopened:
            output[output_type] = output_file_reopened.read()
        os.unlink(output_file)

    _process_outputs(output)
        
    if process_info:
        output['process_info'] = ret
            
    if not output:
        output = None

    return output
    

def run( tileset: dict, extraparams: dict, outputopts=None, ui=False):
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
    do those.  These will manage the output, and return the data in usable form.

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

    tileset_copy = copy.deepcopy(tileset)
    tileset_copy['xgrowargs'].update(extraparams)
    if ui:
        tileset_copy['xgrowargs']['window']=True
    else:
        tileset_copy['xgrowargs']['window']=False
    tilestring = stxg.to_xgrow(tileset_copy)
    
    return run_old(tilestring, {}, outputopts)
