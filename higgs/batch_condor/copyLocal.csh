#!/bin/csh

set dir  = $1
set username = $2

ls /pnfs/cms/WAX/11/store/user/$username/$dir/nuTuple_*.root  > DataFiles.txt

if (! -d ~/nobackup/$dir) then
   mkdir ~/nobackup/$dir
endif

set count = 1
set nFiles = `cat DataFiles.txt | wc -l`

while($count <= $nFiles)
   set filepath = `cat DataFiles.txt | head -n $count | tail -1`
   set filename = `cat DataFiles.txt | head -n $count | tail -1 | cut -d "/" -f 11`
   echo "Copying store/user/$username/$dir/$filename to ~/nobackup/$dir/"
   dccp $filepath ~/nobackup/$dir/.  
   @ count ++
end 

rm DataFiles.txt
