
# usage: run_np.bash trials 
# Gse_min Gse_max Eps_min Eps_max
#
#  delta_Gse = 0.5
#  delta_Eps = 0.03

if [ -n "$1" ]; then num_trials=$1; else num_trials=1;  fi


gse=10.6;
gse_min=6.0;
gse_step=.5;
gmc=10;

smax=2000

xgrow="/home/wrightc/SURF06/xgrow/xgrow"
tileset="/home/wrightc/SURF06/xgrow/zztilesets"
datafile="/home/wrightc/SURF06/xgrow/zzscripts/data/06_28_06_2"
seed="/home/wrightc/SURF06/xgrow/zzseeds/6nomelt-0000-"

echo
echo "Running 6wide np 2p 3p 4p $num_trials times from gse = $gse to $gse_min in $gse_step increments"
echo

while [[ $gse > $gse_min ]]; do

for ((trial=1; trial <= num_trials; trial++)); do
 echo
 echo "Begining round $trial. gse = $gse"

 echo
 echo "np"
 nice $xgrow $tileset/6np_nomelt.tiles Gmc=$gmc Gse=$gse datafile="$datafile/error_np" size=512 smax=$smax importfile=$seed clean_cycles=1 -nw

 echo
 echo "2p"
 nice $xgrow $tileset/6p2_nomelt.tiles Gmc=$gmc Gse=$gse datafile="$datafile/error_2p" size=512 smax=$smax importfile=$seed clean_cycles=1 -nw

 echo
 echo "3p"
 nice $xgrow $tileset/6p3_nomelt_error.tiles Gmc=$gmc Gse=$gse datafile="$datafile/error_3p" size=512 smax=$smax importfile=$seed clean_cycles=1 -nw

 echo
 echo "4p"
 nice $xgrow $tileset/6p4_nomelt.tiles Gmc=$gmc Gse=$gse datafile="$datafile/error_4p" size=512 smax=$smax importfile=$seed clean_cycles=1 -nw

done

gse=$(echo "$gse - $gse_step" | bc) 

done
