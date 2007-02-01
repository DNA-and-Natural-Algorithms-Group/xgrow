#
# January 2007
#  this file will use the actual weave (as in the one used in experiments
#  it will also use a more realistic gmc = 13
#  and a larger range of Taos (closer to 2)
#

if [ -n "$1" ]; then num_trials=$1; else num_trials=1;  fi

gse=13;
gse_max=7.2;
gse_min=6.5100;
gmc=13;

smax=3000;
tmax=2000;
xgrow="/home/wrightc/SURF06/xgrow/xgrow";
tileset="/home/wrightc/SURF06/xgrow/zztilesets/weave";
datafile="/home/wrightc/SURF06/xgrow/zzscripts/final/errorrates";
seed="/home/wrightc/SURF06/xgrow/zzseeds/weave";


echo "running 6wide weave @num_trials times from gse=$gse to $gse_min";
echo $pn;

x=0
while [[ $x -eq "0" ]]; do
  for ((trial=1; trial <= num_trials; trial++)); do
    echo "---- round $trial of $num_trials. gse=$gse ----"

pn="nice $xgrow $tileset/np.tiles Gmc=$gmc Gse=$gse datafile=$datafile/pn size=512 smax=$smax tmax=$tmax importfile=$seed clean_cycles=1 -nw";
p2="nice $xgrow $tileset/p2.tiles Gmc=$gmc Gse=$gse datafile=$datafile/p2 size=512 smax=$smax tmax=$tmax importfile=$seed clean_cycles=1 -nw";
p3="nice $xgrow $tileset/p3.tiles Gmc=$gmc Gse=$gse datafile=$datafile/p3 size=512 smax=$smax tmax=$tmax importfile=$seed clean_cycles=1 -nw";
p4="nice $xgrow $tileset/p4.tiles Gmc=$gmc Gse=$gse datafile=$datafile/p4 size=512 smax=$smax tmax=$tmax importfile=$seed clean_cycles=1 -nw";

    $pn

    $p2

    $p3
    $p3

    $p4
    $p4

  done
  gse=$(echo "$gse - 0.4" | bc)
  x=$(echo "$gse_min > $gse" | bc)
done


# add some more sample density to T near 2

gse=8.8
x=0
while [[ $x -eq "0" ]]; do
  for ((trial=1; trial <= num_trials; trial++)); do
    echo "--=- round $trial of $num_trials. gse=$gse -=--"

pn="nice $xgrow $tileset/np.tiles Gmc=$gmc Gse=$gse datafile=$datafile/pn size=512 smax=$smax tmax=$tmax importfile=$seed clean_cycles=1 -nw";
p2="nice $xgrow $tileset/p2.tiles Gmc=$gmc Gse=$gse datafile=$datafile/p2 size=512 smax=$smax tmax=$tmax importfile=$seed clean_cycles=1 -nw";
p3="nice $xgrow $tileset/p3.tiles Gmc=$gmc Gse=$gse datafile=$datafile/p3 size=512 smax=$smax tmax=$tmax importfile=$seed clean_cycles=1 -nw";
p4="nice $xgrow $tileset/p4.tiles Gmc=$gmc Gse=$gse datafile=$datafile/p4 size=512 smax=$smax tmax=$tmax importfile=$seed clean_cycles=1 -nw";

    $pn

    $p2

    $p3
    $p3
    $p3

    $p4
    $p4
    $p4

  done
  gse=$(echo "$gse - 0.4" | bc)
  x=$(echo "$gse_min > $gse" | bc)
done

