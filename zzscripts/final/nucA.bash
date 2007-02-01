#
# January 2007
#  this file will measure our simulated nucleation Tao. it uses GMC
#  consistent with other simulations
#
#  this file will use the actual weave (as in the one used in experiments
#  it will also use a more realistic gmc = 13
#  and a larger range of Taos (closer to 2)
#

if [ -n "$1" ]; then num_trials=$1; else num_trials=1;  fi


### seed ###

gse=13;
gmc=13;

smax=80;
#tmax=2000;
#xgrow="/home/wrightc/SURF06/xgrow/xgrow/";
xgrow="/research/stow/xgrow-2.1/xgrow";
tileset="/home/wrightc/SURF06/xgrow/zztilesets/weave";
datafile="/home/wrightc/SURF06/xgrow/zzscripts/final/nucA/seed";
seed="/home/wrightc/SURF06/xgrow/zzseeds/weave";

pn="nice $xgrow $tileset/np.tiles Gmc=$gmc datafile=$datafile/np size=512 smax=$smax importfile=$seed linan=59.4,0.177,40,12,3600 -nw";
p2="nice $xgrow $tileset/p2.tiles Gmc=$gmc datafile=$datafile/p2 size=512 smax=$smax importfile=$seed linan=59.4,0.177,40,12,3600 -nw";
p3="nice $xgrow $tileset/p3.tiles Gmc=$gmc datafile=$datafile/p3 size=512 smax=$smax importfile=$seed linan=59.4,0.177,40,12,3600 -nw";
p4="nice $xgrow $tileset/p4.tiles Gmc=$gmc datafile=$datafile/p4 size=512 smax=$smax importfile=$seed linan=59.4,0.177,40,12,3600 -nw";


echo "running nucleationA SEED @num_trials times"


#for ((trial=1; trial <= num_trials; trial++)); do
#  echo "...round $trial..."
#  $pn
#  $p2
#  $p3
#  $p4
#done



### no seed ###

gse=4;
fsmax=80;
size=64;
tmax=20000;
datafile="/home/wrightc/SURF06/xgrow/zzscripts/final/nucA/noseed";

pn="nice $xgrow $tileset/np.tiles Gmc=$gmc Gse=$gse datafile=$datafile/np size=$size fsmax=$fsmax tmax=$tmax linan=59.4,0.177,40,12,3600 tinybox=9.183e-18 -nw";
p2="nice $xgrow $tileset/p2.tiles Gmc=$gmc Gse=$gse datafile=$datafile/p2 size=$size fsmax=$fsmax tmax=$tmax linan=59.4,0.177,40,12,3600 tinybox=9.183e-18 -nw";
p3="nice $xgrow $tileset/p3.tiles Gmc=$gmc Gse=$gse datafile=$datafile/p3 size=$size fsmax=$fsmax tmax=$tmax linan=59.4,0.177,40,12,3600 tinybox=9.183e-18 -nw";
p4="nice $xgrow $tileset/p4.tiles Gmc=$gmc Gse=$gse datafile=$datafile/p4 size=$size fsmax=$fsmax tmax=$tmax linan=59.4,0.177,40,12,3600 tinybox=9.183e-18 -nw";


echo "running nucleationA NO SEED @num_trials times from gse=$gse to $gse_min"
echo $pn
#$pn

for ((trial=1; trial <= num_trials; trial++)); do
  echo "...round $trial..."
  $pn
  $p2
  $p3
  $p4
done


