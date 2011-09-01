#!/bin/csh

set source = $1
set count = 1
set nJobs = `cat $1 | wc -l`

while($count <= $nJobs)
   cp condorSettings batchSubmit
   set arg = `cat $1 | head -n $count | tail -1`
   sed -i "s/ARGS/$arg/g" batchSubmit
   condor_submit batchSubmit
   rm batchSubmit
   @ count ++
end

