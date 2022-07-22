# Cluster running code for Xgrow using STXG input.

import logging as log
import pandas as pd
import numpy as np
from io import StringIO
from . import parseoutput as xgo


def _run_xgrow(params):
    # imports go here because the nodes might not have imported yet.
    import subprocess
    import tempfile
    import os
    import yaml

    xg, paramstring, outtype = params
    # Move to correct directory on cluster. This needs to be done each time
    # because the nodes need it. In this case, there is a '.xgrowcluster' in
    # the node's user's home folder has the config for each node yaml-encoded.
    with open(os.environ["HOME"] + "/.xgrowcluster") as configfile:
        config = yaml.load(configfile)
    os.chdir(config["xgrowpath"])
    from . import stxg

    with tempfile.NamedTemporaryFile(
        dir=config["temppath"], delete=False, mode="w"
    ) as tilefile:
        tilestring = stxg.to_xgrow(xg)
        tilefile.write(tilestring)

    with tempfile.NamedTemporaryFile(mode="r", dir=config["temppath"]) as datafile:
        z = subprocess.check_call(
            "./xgrow {0} -nw {1}={2} {3}".format(
                tilefile.name, outtype, datafile.name, paramstring
            ),
            shell=True,
        )
        data = datafile.read()

    os.remove(tilefile.name)

    return data


class XgrowResult:
    def __init__(self, res, outtype):
        self.res = res
        self.ot = outtype

    def ready(self):
        return self.res.ready()

    def result(self):

        if self.ot == "final":
            return pd.DataFrame(
                np.loadtxt(StringIO(u"".join(self.res.result()))),
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
        elif self.ot == "array":
            return [xgo.loadflake(x) for x in self.res.result()]
        elif self.ot == "trace":
            dframes = [
                pd.DataFrame(
                    np.loadtxt(StringIO(u"".join(x))),
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
                for x in self.res.result()
            ]
            return dframes


class XgrowCluster:
    """
    A class controlling an IPCluster used for Xgrow simulations.
    """

    def __init__(self, pool):
        self.pool = pool
        self.client = pool.client
        log.debug("Connected to IPCluster with {0} nodes.".format(len(self.client.ids)))

    def do_sims(self, tilestring, paramstring, outtype, num=50):
        """
        Runs a set of simulations with the same parameters on the cluster,
        and return the datafile output for each.

        * tilestring: the tileset file, as a string.
        * paramstring: options to be appended to the xgrow command line.
          -nw and datafile are included by default.
        * num (=50 default): number of runs to run.
        """
        outparam = {"final": "datafile", "trace": "tracefile", "array": "arrayfile"}[
            outtype
        ]
        commands = num * [(tilestring, paramstring, outparam)]

        output = self.pool.map(_run_xgrow, commands)
        return XgrowResult(output, outtype)


class SimRunSingleMismatch:
    """A class for sets of simulations with single mismatch interactions, varying
    the strength of the interaction for each set.
    """

    def __init__(
        self,
        xgcluster,
        baseset,
        extraparams,
        num=50,
        defvals=[0.0, 0.025, 0.05, 0.075, 0.1, 0.2, 0.4, 0.5, 0.8],
    ):
        """baseset is a tileset file as a string, in a way such that glue info can be
        appended to it (no params at the end).  extraparams are all the extra
        parameters

        """
        self.num = num
        self.baseset = baseset
        self.extraparams = extraparams
        self.xgc = xgcluster
        self.res = {}
        self.defvals = defvals

    def runfunction(self, g1, g2, strength, n=50):
        log.debug("In RunFunction")
        return self.xgc.do_sims(
            self.baseset + "\ng({},{})={}\n".format(g1, g2, strength),
            self.extraparams,
            num=n,
        )

    def add(self, name, g1, g2, vals=None):
        "Add a new simulation set, with interactions between g1 and g2"
        self.res[name] = {}
        if not vals:
            vals = self.defvals
        for gval in vals:
            self.res[name][gval] = self.runfunction(g1, g2, gval, n=self.num)

    @property
    def progress(self):
        """Fractional progress of run. Note that this doesn't use the
        AsyncResult.progress method because in normal ipython 0.13 / 0.14 this
        can sometimes hang for 10 seconds owing to a timeout problem.

        """
        done = 0
        total = 0
        for x in self.res.values():
            log.info("{} / {}".format(done, total))
            for y in x.values():
                total += len(y)
                done += len(y) - len(set(y.msg_ids).intersection(y._client.outstanding))
        return (1.0 * done / total, done, total)

    def map(self, function):
        return dict(
            (n, np.array([[k, function(x)] for k, x in self.items()])) for n in self.res
        )
