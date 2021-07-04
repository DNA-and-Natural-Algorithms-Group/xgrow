import sys
import re
from typing import List
import yaml
from xgrow import xgrow

__all__: List[str] = []

TF_KEYS = ["wander", "pause", "movie", "zero_bonds", "window"]

parsed_traditional = {
    "no_fission": lambda x: ("fission", "off"),
    "fission": lambda x: ("fission", "on"),
    "chunk_fission": lambda x: ("fission", "chunk"),
    "importfile": lambda x: (
        "importfile",
        x.group(3) if x.group(3) else "xgrow_export_output",
    ),
}

parsed_new = {
    "w": lambda x: ("window", not bool(x.group(1))),
    "nw": lambda x: ("window", bool(x.group(1))),
}


def print_help():
    print(
        """
xgrow [NEW OPTIONS] TILESET_FILE [NEW+OLD OPTIONS]

The new Python xgrow wrapper, designed to run both older xgrow tilesets with
xgrow command line arguments, and stxg tilesets originally run by xgrow-wrap.
Note that this does *not* currently run Alhambra tilesets: that functionality
will likely be moved into Alhambra itself.

TILESET_FILE can be either an stxg file, or an older Xgrow-type file.

Options to this script can be provided in two different ways:

    - as normal command line options, using -KEY VALUE or -KEY=VALUE for options
      expecting a value, or -KEY or -no-KEY for true/false options.

    - for backwards compatibility, as old Xgrow options without a dash, but only
      after TILESET_FILE has been specified.  This means that this wrapper of
      xgrow is completely compatible with command lines for the actual Xgrow,
      with the exception of "xgrow --" to show help.

A full list of options can be printed with xgrow --help-options, and is also in the
documentation.  In general, any option that can be used in an stxg or Xgrow file can
be used here, though there are also the shorthands -w and -nw for -window and -no-window.
    """
    )


def print_options():
    print(
        """
These options are currently taken from xgrow itself.  They need to be updated.

  block=  display block size, 1...10
  size=   field side length (power-of-two) [default 256]
  rand=   random number seed
  T=      threshold T (relative to Gse) for irreversible Tile Assembly Model
  k=      hybridization rate constant (/sec)
  Gmc=    initiation free energy  (units kT)
  Gse=    interaction free energy per binding
  Gmch=   initiation free energy  for hydrolyzed units
  Gseh=   interaction free energy for hydrolyzed units
  Ghyd=   free energy of hydrolysis
  Gas=    activation energy for spontaneous hydrolysis
  Gam=    activation energy for mismatched sticky ends
  Gae=    activation energy for unmatched sticky ends
  Gah=    activation energy for hydrolyzed neighbors
  Gao=    delta a. e. for output vs input-triggers hydrolysis
  Gfc=    log concentration of flakes (otherwise no depletion)
  stoic=n[s]            set stoichiometry of tile n to relative value s
  anneal=g,t            anneal Gse from g to given Gse with time constant t
  linan=h,s,C0,Cfin,dt  do a linear anneal, where delta G at each temperature is
                        governed by (Delta) h and (Delta) s in kcals per mol.
                        Go from C0 > Cfin changing temps 1 degree every dt secs
                        (in incremenets of 0.1C).  Ignores Gse value.
  seed=i,j,n            seed tile type n at position i,j
  tinybox=V             use dynamic flakes in a box of volume V (in liters).
  addflakes=i,j,n:N@Gfc simulate N separate flakes
  stripe=o[:p,w]*       width w stripe with p errors, offset o
  wander                wandering `seed' designation
  fission               can tile be removed if two flakes result?
  no_fission             the answer is no [default]
  chunk_fission         allow pairs & 2x2 blocks to dissociate as one (implies fission)
  blast_rate_alpha      square kxk holes are blasted with this per-tile rate for 1x1 [default=0]
                         (rate relative to tile addition, i.e. scaled by total concentration)
  blast_rate_beta        rate scales as 1/k^beta [default=4]
  blast_rate_gamma       rate also scales as exp(-gamma*(k-1)) [default=0]
  zero_bonds            can tiles be added if they bond with 0 strength?
  no_zero_bonds          the answer is no [default]
  min_strength=         set minimum bond energy below which attachments are considered incorrect
  periodic              periodic boundary conditions
  -linear               simulate linear A B tiles, write errs > stdout
  -nw                   no X window (only if ?max set)
  update_rate=          update display every so-many events
  tracefile=            append datafile info (see below) EVERY so-many events
  movie                 export MATLAB-format flake array information EVERY so-many events
  tmax=                 quit after time t has passed
  emax=                 quit after e events have occurred
  smax=                 quit when the fragment or total size of fragments is size s
  mmax=                quit when there are at least mm mismatches by xgrow's interpretation
  smin=                 quit when the fragment or total size of fragments goes to or below size s
  fsmax=                 quit when a single fragment reaches size s
  untiltiles=           quit when all (numbered) tiles in the comma-delineated list are in the assembly
  clean_cycles=         at end, remove how many layers of weakly attached tiles?
                        [default=0]
  clean_X=              for cleaning, minimal ratio of off-rate to on-rate [default=1.0]
  fill_cycles=          at end, add how many layers of strongly attached tiles?
                        [default=0]
  fill_X=               for filling, minimal ratio of off-rate to on-rate [default=1.0]
  error_radius=         when writing to file, #mismatches counts only those for which
                        all surounding tiles are present (after clean/fill) [default=0]
  repair_unique_T=      alternative clean/fill, called Rx: remove mismatches, fill in interior sites
                        if there is a unique strength-T tile, then fill in by strongest tile
  datafile=             append Gmc, Gse, ratek, time, size, #mismatched se, events, perimeter, dG, dG_bonds for each flake
  arrayfile=            output MATLAB-format flake array information on exit (after cleaning)
  exportfile=           on-request output of MATLAB-format flake array information
                        [defaults to 'xgrow_export_output']
  importfile=FILENAME   import all flakes from FILENAME.
  importfile            import all flakes from xgrow_export_output.
  pause                 start in paused state; wait for user to request simulation to start.
  testing               run automated tests instead of a simulation.
    """
    )


def main():
    # Why doesn't this use argparse or click? Because, in order to maintain
    # backward compatibility, it can't... What we want is the following:
    # 1. Tileset must come before non-dashed options, but otherwise can be anywhere.
    # 2. Dashed options can go before or after the tileset.
    # 3. Dashed options can take split or = form.
    # 4. Undashed options must take xgrow (=) form.

    raw_args = sys.argv[1:]
    args = {}
    pos = 0
    tilepath = None
    while pos < len(raw_args):
        # First try to parse as an argument:
        m = re.match(r"-{1,2}(no[-_])?([^=\s]+)(=)?(\S+)?", raw_args[pos])
        if m and m.group(2) in ["h", "help"]:
            print_help()
            sys.exit(0)
        if m and m.group(2) in ["help-options"]:
            print_options()
            sys.exit(0)
        if m and m.group(2) in parsed_new.keys():
            args.__setitem__(*parsed_new[m.group(2)](m))
            pos += 1
        elif m and m.group(4):
            args[m.group(2)] = m.group(4)
            pos += 1
        elif m:
            if m.group(2) not in TF_KEYS:
                args[m.group(2)] = raw_args[pos + 1]
                pos += 2
            else:
                if m.group(1):
                    args[m.group(2)] = False
                else:
                    args[m.group(2)] = True
                pos += 1
        elif not tilepath:
            tilepath = raw_args[pos]
            pos += 1
        else:
            m = re.match(r"([^=\s]+)(=)?(\S+)?", raw_args[pos])
            if not m:
                print(f"Can't parse {raw_args[pos]}.")
                sys.exit(1)
            if m.group(1) in parsed_traditional.keys():
                args.__setitem__(*parsed_traditional[m.group(1)](m))
            elif m.group(3):
                args[m.group(1)] = m.group(3)
            else:
                args[m.group(1)] = True
            pos += 1

    if not tilepath:
        print("No tileset file specified!\nUse -h or --help for help.")
        sys.exit(1)

    fd = open(tilepath, "r")
    for line in fd:
        if line[0] == "%":
            continue
        else:
            if line[0:4] == "tile" or line[0:3] == "num":
                fd.seek(0)
                xgrow.run_old(fd.read(), args)
                break
            else:
                fd.seek(0)
                xgrow.run(yaml.safe_load(fd), args)
                break


if __name__ == "__main__":
    main()
