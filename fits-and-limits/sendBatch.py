#!/usr/bin/env python
import sys

sys.path.append("../batch-condor")
import BatchMaster as b
cfg = b.JobConfig

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog subdir")
parser.add_option("-j", dest="jobs",   default=2, help="Number od jobs")
parser.add_option("-t", dest="trials", default=5, help="Number of trials per job")

(options, args) = parser.parse_args()


if len(args) < 1:
    parser.print_usage()
    exit(1)
    
''' Specify parameters '''
executable  = 'batchJob.sh'

version    = args[0]
outputPath = version+'-bias'

inputSamples = []
inputSamples.append(cfg('None', 'None', int(options.jobs), str(options.trials) +' 120'))
inputSamples.append(cfg('None', 'None', int(options.jobs), str(options.trials) +' 125'))
inputSamples.append(cfg('None', 'None', int(options.jobs), str(options.trials) +' 130'))
inputSamples.append(cfg('None', 'None', int(options.jobs), str(options.trials) +' 135'))
inputSamples.append(cfg('None', 'None', int(options.jobs), str(options.trials) +' 140'))
inputSamples.append(cfg('None', 'None', int(options.jobs), str(options.trials) +' 145'))
inputSamples.append(cfg('None', 'None', int(options.jobs), str(options.trials) +' 150'))


batcher = b.BatchMaster(inputSamples, outputPath, shortQueue = False, stageDir = '../StageBatch_'+version,
                        executable = executable, prefix = '', bias=True)
print "Submitting to batch?"
batcher.submit_to_batch()
            
