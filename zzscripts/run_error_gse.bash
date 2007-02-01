
# usage: run_np.bash trials 
# Gse_min Gse_max Eps_min Eps_max
#
#  delta_Gse = 0.5
#  delta_Eps = 0.03

if [ -n "$1" ]; then num_trials=$1; else num_trials=1;  fi


gse=6.5300;
gse_max=7.2;
gmc=13;

smax=2000;
tmax=2000;
xgrow="/home/wrightc/SURF06/xgrow/xgrow";
tileset="/home/wrightc/SURF06/xgrow/zztilesets";
datafile="/home/wrightc/SURF06/xgrow/zzscripts/data/error";
seed="/home/wrightc/SURF06/xgrow/zzseeds/6nomelt-0000-";

echo
echo "Running 6wide np 2p 3p 4p $num_trials times from gse = $gse to $gse_min in $gse_step increments"
echo

x=0
while [[ $x -eq "0" ]]; do

for ((trial=1; trial <= num_trials; trial++)); do
 
 echo
 echo "Begining round $trial. gse = $gse"

 echo "---np---"
# nice $xgrow $tileset/6np_nomelt.tiles Gmc=$gmc Gse=$gse datafile="$datafile/1_np" size=512 smax=$smax tmax=$tmax  importfile=$seed clean_cycles=1 -nw

 echo "---2p---"
# nice $xgrow $tileset/6p2_nomelt.tiles Gmc=$gmc Gse=$gse datafile="$datafile/1_2p" size=512 smax=$smax tmax=$tmax importfile=$seed clean_cycles=2 -nw

 echo "---3p---"
 nice $xgrow $tileset/6p3_nomelt_error.tiles Gmc=$gmc Gse=$gse datafile="$datafile/1_3p" size=512 smax=$smax tmax=$tmax importfile=$seed clean_cycles=3 -nw

 echo "---4p---"
 nice $xgrow $tileset/6p4_nomelt.tiles Gmc=$gmc Gse=$gse datafile="$datafile/1_4p" size=512 smax=$smax tmax=$tmax importfile=$seed clean_cycles=4 -nw


 echo "---3p---"
 nice $xgrow $tileset/6p3_nomelt_error.tiles Gmc=$gmc Gse=$gse datafile="$datafile/1_3p" size=512 smax=$smax tmax=$tmax importfile=$seed clean_cycles=3 -nw

 echo "---4p---"
 nice $xgrow $tileset/6p4_nomelt.tiles Gmc=$gmc Gse=$gse datafile="$datafile/1_4p" size=512 smax=$smax tmax=$tmax importfile=$seed clean_cycles=4 -nw

done

gse=$(echo "$gse * 1.01" | bc) 
x=$(echo  "$gse > $gse_max" | bc)

done
