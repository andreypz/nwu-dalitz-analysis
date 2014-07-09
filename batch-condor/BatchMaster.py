import sys, os, subprocess, fileinput, math, tempfile, datetime

class JobConfig():
  '''Class for storing configuration for each dataset'''
  def __init__(self, dataName = 'TEST', inDir = '', nJobs = 1, arguments = 'TEST 0 muon'):
    self._dataName  = dataName
    self._inDir     = inDir
    self._nJobs     = nJobs
    self._args      = arguments

class BatchMaster():
    '''A tool for submitting batch jobs'''
    def __init__(self, configList, outDir = '.', shortQueue = False, stageDir = 'outputBatch', executable = 'execBatch.sh', prefix = 'test', bias=False):
      self._prefix  = prefix
      self._current    = os.path.abspath('.')
      self._stageDir   = stageDir
      self._outDir     = outDir
      self._configList = configList
      self._executable = executable
      self._shortQueue = shortQueue
      self._bias      = bias


    def get_current_time(self):
      '''
      Returns a string of the current time with the format
      '''

      now = datetime.datetime.now()
      currentTime = '{0:02d}{1:02d}{2:02d}_{3:02d}{4:02d}{5:02d}'.format(now.year, now.month, now.day, now.hour, now.minute, now.second)
      return currentTime


    def make_directory(self, filePath, clear = True):
      '''
      Create save path in case it doesnt already exist
      '''
      print "Creating directory", filePath
      if not os.path.exists(filePath):
        print "mkdir", filePath
        os.system('mkdir -p '+filePath)

      if clear and len(os.listdir(filePath)) != 0:
        os.system('rm '+filePath+'/*')


    def split_jobs(self, directory, nJobs):
      '''
      Split jobs by dataset
      '''

      fileList = os.listdir(directory)
      nFiles = len(fileList)

      # Split files to requested number.  Cannot exceed
      # the number of files being run over.
      if nJobs > nFiles:
        nJobs = nFiles

      nFilesPerJob = int(math.ceil(float(nFiles)/float(nJobs)))
      fileSplit = [fileList[i:i+nFilesPerJob] for i in range(0, len(fileList), nFilesPerJob)]

      return fileSplit


    def make_batch_lpc(self, cfg, count, sourceFiles):
        '''
        Prepares for submission to lpc.  Does the following:

        1. Generates input_files.txt with files to run over
        2. Write batch configuration file
        4. ????
        5. Poor house :(
        '''
        if not self._bias:
          ## make file with list of inputs ntuples for the analyzer
          input = open('{0}/input_{1}_{2}.txt'.format(self._stageDir, cfg._dataName, str(count+1)), 'w')

          path = cfg._inDir
          for i,source in enumerate(sourceFiles):
            input.write(path+'/'+source+'\n')
          input.close()
        else:
          input =None

        ## Writing the batch config file
        batch_tmp = open('{0}/.batch_tmp_{1}_{2}'.format(self._stageDir, cfg._dataName, count+1), 'w')
        batch_tmp.write('Universe              = vanilla\n')

        ## TESTING CONDORG ## --> Doesn't seem to work :(
        #batch_tmp.write('Universe              = globus\n')
        #batch_tmp.write('globusscheduler       = cmsosgce.fnal.gov/jobmanager-condor\n')
        if self._bias:
          batch_tmp.write('Arguments             = {0} {1} {2}\n'.format(self._outDir, count, cfg._args))
          print "Creating condor cfg for:",self._outDir, count+1,cfg._args
        else:
          batch_tmp.write('Arguments             = {0} {1} {2} {3}\n'.format(self._outDir, count+1, cfg._dataName, cfg._args))

        batch_tmp.write('Executable            = {0}\n'.format(self._executable))
        batch_tmp.write('Should_Transfer_Files = YES\n')
        batch_tmp.write('WhenToTransferOutput  = ON_EXIT\n')
        if self._bias:
          batch_tmp.write('Transfer_Input_Files  = {0}/source.tar.gz\n'.format(self._stageDir))
        else:
          batch_tmp.write('Transfer_Input_Files  = {0}/source.tar.gz, {0}/input_{1}_{2}.txt\n'.format(self._stageDir, cfg._dataName,  str(count+1)))
          batch_tmp.write('Transfer_Output_Files = out_synch_{0}_{1}.txt, hhhh_{2}_{3}.root\n'.format(
              cfg._dataName,  str(count+1), cfg._dataName,  str(count+1)))

        batch_tmp.write('Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000\n')
        batch_tmp.write('Notification = Never\n')
        if self._shortQueue:
          batch_tmp.write('+LENGTH = "SHORT"\n')
        batch_tmp.write('Output    = '+self._stageDir+'/report_'+cfg._dataName+'_'+str(count)+'_$(Cluster)_$(Process).stdout\n')
        batch_tmp.write('Error     = '+self._stageDir+'/report_'+cfg._dataName+'_'+str(count)+'_$(Cluster)_$(Process).stderr\n')
        batch_tmp.write('Log       = '+self._stageDir+'/report_'+cfg._dataName+'_'+str(count)+'_$(Cluster)_$(Process).log   \n')
        batch_tmp.write('Queue\n')
        batch_tmp.close()


    def submit_to_batch(self, bSystem = 'lpc'):
        '''
        Submits batch jobs to batch.  Currently only works
        for lpc batch system, but should be updated for more
        general use
        '''

        self._stageDir  = self._stageDir + '_' + self.get_current_time()
        self._outDir    = self._outDir + '/' + self._prefix
        self.make_directory(self._stageDir, clear=False)
        self.make_directory(self._outDir, clear=False)
        os.system('cp {0} {1}'.format(self._executable, self._stageDir))
        os.chdir(self._stageDir)


        print "Making a tar-ball", self._stageDir
        if not self._bias:
          os.system('tar czf {0}/source.tar.gz ../zgamma ../testAnalyzer ../src ../interface ../plugins ../data 2> /dev/null'.format(self._stageDir))
        else:
          os.system('tar czf {0}/source.tar.gz ../fits-and-limits  2> /dev/null'.format(self._stageDir))


        if bSystem=='lpc' or bSystem=='tier3':
          for cfg in self._configList:
            if self._bias:
              for i in xrange(cfg._nJobs):
                print 'Submitting: ', i
                self.make_batch_lpc(cfg, i, None)
                subprocess.call('condor_submit .batch_tmp_{0}_{1}'.format(cfg._dataName, i+1), shell=True)
            else:
              sourceFiles = self.split_jobs(cfg._inDir, cfg._nJobs)

              for i, source in enumerate(sourceFiles):
                self.make_batch_lpc(cfg, i, source)
                subprocess.call('condor_submit .batch_tmp_{0}_{1}'.format(cfg._dataName, i+1), shell=True)
                # subprocess.call('condor_submit {0}/.batch_tmp_{1}_{2}'.format(self._stageDir, cfg._dataName, i+1), shell=True)
        else:
          print "Non LPC or Tier3?  Sorry, this is not supported yet"




