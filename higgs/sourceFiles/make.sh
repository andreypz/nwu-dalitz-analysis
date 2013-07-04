#!/bin/sh

path=$1
fname=$2
user="andreypz"
pnfs="/pnfs/cms/WAX/11/store/user"
sample=${path:26}
ls -1 -d  $pnfs/$user/$sample/*.root > $fname 

r1="/pnfs/"
r2="dcap://cmsdca1.fnal.gov:24140/pnfs/fnal.gov/usr/"

sed -i "s|$r1|$r2|g" $fname

