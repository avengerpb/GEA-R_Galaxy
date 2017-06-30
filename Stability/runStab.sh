input=$1;
cp $input input.csv
lsmean=$2;
expdes=$3;
traits=$4;
staparms=$5;
selected=$6;
report=$7;
adapt=$8;
plot=$9;
plot2=${10};
tmp=$RANDOM

{

echo 'datos1 = read.csv("input.csv")'
echo 'ExpDes1 <- c("'$expdes'")'
echo 'BiplotFormat1 <- c("pdf")'
echo 'traits1 <- c("'$traits'")'
echo 'LSMeans_par1 <-c("'$lsmean'")'
echo 'sta.parms1 <-c('$staparms')'
echo 'tmp1 <- c("'$tmp'")'
echo 'selected1 <- c("'$selected'")'
echo 'source("/home/dquoc/galaxy/tools/GEA-R/Stability/Stability.R")'

echo 'STABILITY()'
} >run.R

Rscript ./run.R
dname='Output'
mv $dname*$tmp*/Stab*  $report
mv $dname*$tmp*/Adapt* $adapt
mv $dname*$tmp*/PlotCV* $plot
mv $dname*$tmp*/PlotEb* $plot2
