# usage: run_np.bash trials 
# Gse_min Gse_max Eps_min Eps_max
#
#  delta_Gse = 0.5
#  delta_Eps = 0.03

if [ -n "$1" ]; then num_trials=$1; else num_trials=1;  fi
# if [ -n "$2" ]; then gse_min=$(echo "$2+0" | bc); else gse_min=7.5;  fi
# if [ -n "$3" ]; then gse_max=$(echo "$3+0" | bc); else gse_max=$gse_min;  fi
# if [ -n "$4" ]; then eps_min=$(echo "$4+0" | bc); else eps_min=.15; fi
# if [ -n "$5" ]; then eps_max=$(echo "$5+0" | bc); else eps_max=$eps_min; fi

echo "Running 6ide-2bit-proofreading $num_trials times"
# Gse = $gse_min ... $gse_max Eps = $eps_min ... $eps_max"
echo

gse=5.4;
gmc=10;

# gse=$gse_min;
# eps=$eps_min;
# deps=.03;
# dgse=.5;

for ((trial=1; trial <= num_trials; trial++)); do
 echo
 echo "Begining round $trial."

# gse_stop=0;
# while [[ $gse_stop -eq "0" ]]; do
#  eps_stop=0;
#  while [[ $eps_stop -eq "0" ]]; do
#   gmc=$(echo "2*$gse - $eps" | bc)
#   echo "sierpinski2x2 Gmc=$gmc Gse=$gse"

   echo "-0000-"

   nice xgrow /home/wrightc/SURF06/xgrow/zztilesets/6p2TileFile.tiles -nw Gmc=$gmc Gse=$gse datafile=/home/wrightc/SURF06/xgrow/zzscripts/data/2p0000 size=512 smax=3000 importfile=/home/wrightc/SURF06/xgrow/zztilesets/6seed-0000-

#   echo "-0011-"

#   nice xgrow /home/wrightc/SURF06/xgrow/zztilesets/6p2TileFile.tiles -nw Gmc=$gmc Gse=$gse datafile=/home/wrightc/SURF06/xgrow/zzscripts/data/2p0011 size=512 smax=3000 importfile=/home/wrightc/SURF06/xgrow/zztilesets/6seed-0011-
 
#   echo "-1100-"

#   nice xgrow /home/wrightc/SURF06/xgrow/zztilesets/6p2TileFile.tiles -nw Gmc=$gmc Gse=$gse datafile=/home/wrightc/SURF06/xgrow/zzscripts/data/2p1100 size=512 smax=3000 importfile=/home/wrightc/SURF06/xgrow/zztilesets/6seed-1100-

#    echo "-1111-"

#   nice xgrow /home/wrightc/SURF06/xgrow/zztilesets/6p2TileFile.tiles -nw Gmc=$gmc Gse=$gse datafile=/home/wrightc/SURF06/xgrow/zzscripts/data/2p1111 size=512 smax=3000 importfile=/home/wrightc/SURF06/xgrow/zztilesets/6seed-1111-


#   eps=$(echo "$eps + $deps" | bc) 	
#   eps_stop=$(echo "$eps > $eps_max" | bc)
#   if [[ $eps_stop -eq "1" ]]; then eps=$eps_min; fi
#  done
#  gse=$(echo "$gse + $dgse" | bc) 	
#  gse_stop=$(echo "$gse > $gse_max" | bc)
#  if [[ $gse_stop -eq "1" ]]; then gse=$gse_min; fi
# done

done


