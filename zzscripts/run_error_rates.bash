
# usage: run_np.bash trials 
# Gse_min Gse_max Eps_min Eps_max
#
#  delta_Gse = 0.5
#  delta_Eps = 0.03

if [ -n "$1" ]; then num_trials=$1; else num_trials=1;  fi


gse=5;
gse_max=10;
gse_step=1;
gmc=10;

echo
echo "Running 6wide np 2p 3p 4p $num_trials times"
echo

while [[ $gse -lt $gse_max ]] do

for ((trial=1; trial <= num_trials; trial++)); do

 echo
 echo "Begining round $trial. gse = $gse"

 echo
 echo "np"

   nice /home/wrightc/SURF06/xgrow/xgrow /home/wrightc/SURF06/xgrow/zztilesets/6np_nomelt.tiles Gmc=$gmc Gse=$gse datafile=/home/wrightc/SURF06/xgrow/zzscripts/data/06_26_06/error_np size=512 smax=3000 importfile=/home/wrightc/SURF06/xgrow/zzseeds/6nomelt-0000- clean_cycles=1 -nw

 echo
 echo "2p"

   nice /home/wrightc/SURF06/xgrow/xgrow /home/wrightc/SURF06/xgrow/zztilesets/6p2_nomelt.tiles Gmc=$gmc Gse=$gse datafile=/home/wrightc/SURF06/xgrow/zzscripts/data/06_26_06/error_2p size=512 smax=3000 importfile=/home/wrightc/SURF06/xgrow/zzseeds/6nomelt-0000- clean_cycles=2 -nw
 
 echo
 echo "3p"

   nice /home/wrightc/SURF06/xgrow/xgrow /home/wrightc/SURF06/xgrow/zztilesets/6p3_nomelt_error.tiles Gmc=$gmc Gse=$gse datafile=/home/wrightc/SURF06/xgrow/zzscripts/data/06_27_06/error_3p size=512 smax=3000 importfile=/home/wrightc/SURF06/xgrow/zzseeds/6nomelt-0000- clean_cycles=3 -nw

 echo
 echo "4p"

   nice /home/wrightc/SURF06/xgrow/xgrow /home/wrightc/SURF06/xgrow/zztilesets/6p4_nomelt.tiles Gmc=$gmc Gse=$gse datafile=/home/wrightc/SURF06/xgrow/zzscripts/data/06_27_06/error_4p size=512 smax=500 importfile=/home/wrightc/SURF06/xgrow/zzseeds/6nomelt-0000- clean_cycles=4 -nw

done

gse = (echo "$gse + $gse_step" | bc) 

done
