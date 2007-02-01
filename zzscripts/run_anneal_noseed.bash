
# usage: run_np.bash trials 
# Gse_min Gse_max Eps_min Eps_max
#
#  delta_Gse = 0.5
#  delta_Eps = 0.03

if [ -n "$1" ]; then num_trials=$1; else num_trials=1;  fi


gse=10;
g=5;
gmc=10;
smax=3000;

xgrow="/home/wrightc/SURF06/xgrow/xgrow";
tileset="/home/wrightc/SURF06/xgrow/zztilesets";
datafile="/home/wrightc/SURF06/xgrow/zzscripts/data/anneal";
#seed="/home/wrightc/SURF06/xgrow/zzseeds/6nomelt-0000-";

echo
echo "Running 6wide np 2p 3p 4p $num_trials times, $t hours to $tmax hours, from gse = $g to $gse"
echo


for ((t=10; t <=50; t=t+5)); do
 
 for ((trial=1; trial <= num_trials; trial++)); do
 
  echo
  echo "Begining round $trial. t = $t"

  tsec=$(echo "$t * 60" | bc); 
  tend=$(echo "tsec * 5" | bc);

  echo "---np---"
  nice $xgrow $tileset/6np_nomelt.tiles Gmc=$gmc Gse=$gse datafile="$datafile/a2_ns_np" size=512 smax=$smax tmax=$tend clean_cycles=1 anneal=$g,$tsec -nw

  echo "---2p---"
  nice $xgrow $tileset/6p2_nomelt.tiles Gmc=$gmc Gse=$gse datafile="$datafile/a2_ns_2p" size=512 smax=$smax  tmax=$tend clean_cycles=2 anneal=$g,$tsec -nw


  echo "---3p---"
  nice $xgrow $tileset/6p3_nomelt_error.tiles Gmc=$gmc Gse=$gse datafile="$datafile/a2_ns_3p" size=512 smax=$smax  tmax=$tend clean_cycles=3 anneal=$g,$tsec -nw


  echo "---4p---"
  nice $xgrow $tileset/6p4_nomelt.tiles Gmc=$gmc Gse=$gse datafile="$datafile/a2_ns_4p" size=512 smax=$smax  tmax=$tend clean_cycles=4 anneal=$g,$tsec -nw


 done

done
