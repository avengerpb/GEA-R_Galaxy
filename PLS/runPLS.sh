input=$1;
cp $input input.csv
cov=$2;
cp $cov cov.csv
expdes=$3;
traits=$4;
typecov=$5;
selected=$6;
report=$7;
biplot=$8;
tmp=$RANDOM

{

echo 'datos = read.csv("input.csv")'
echo 'file_nameCOV = read.csv("cov.csv")'
echo 'ExpDes <- c("'$expdes'")'
echo 'BiplotFormat <- c("pdf")'
echo 'traits <- c("'$traits'")'
echo 'typecov <- c("'$typecov'")'
echo 'tmp <- c("'$tmp'")'
echo 'selected <- c("'$selected'")'
echo 'source("/home/dquoc/galaxy/tools/GEA-R/PLS/PLS.R")'

echo 'PLSR()'
} >run.R

Rscript ./run.R
dname='OutputPLS'
mv $dname*$tmp*/Final*  $report
mv $dname*$tmp*/Biplot* $biplot








