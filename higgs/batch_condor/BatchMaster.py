import os, subprocess, fileinput, math, tempfile

class JobConfig():
    '''Class for storing configuration for each dataset'''
    def __init__(self, dataName = 'TEST', inDir = '/uscms/home/andreypz/nobackup/TEST', nJobs = 1, arguments = 'TEST 0 muon 2011A', selection = 'muon'):
        self._dataName  = dataName
        self._inDir     = inDir
        self._nJobs     = nJobs
        self._args      = arguments
        self._selection = selection 
        
    def GetMembers(self, verbose = False):
        '''Returns data members of class instance'''
        if verbose:
            print 'dataset: %s \t number of jobs: %i \t triggers: %s' % (self._dataName, self._nJobs, self._triggers)
        return (self._dataName, self._inDir, self._nJobs, self._triggers)

class BatchMaster():
    '''A tool for submitting batch jobs'''
    def __init__(self, configList, outDir = '.', shortQueue = False):
        self._current    = os.path.abspath('.')
        self._outDir     = outDir
        self._configList = configList
        self._executable = 'batchJob.csh'
        self._shortQueue = shortQueue
    
    def MakeDirectory(self, filePath, clear = True):
        '''Create save path in case it doesn't already exist'''
        if not os.path.exists(filePath):
            os.system('mkdir -p '+filePath)
        elif clear:
            os.system('rm -rf '+filePath+'/*')

    def SplitJobs(self, directory, nJobs):
        '''Split jobs by dataset'''
        fileList = os.listdir(directory)
        nFiles = len(fileList)
        
        if nJobs > nFiles:
            nJobs = nFiles
        nFilesPerJob = int(math.ceil(float(nFiles)/nJobs))
        fileSplit = [fileList[i:i+nFilesPerJob] for i in range(0, len(fileList), nFilesPerJob)]

        #print fileList
        return fileSplit

    def MakeExecutable(self, config, sourceFiles, count):
        '''Writes config into executable'''
        infile   = open(self._executable, 'r')
        exec_tmp = tempfile.NamedTemporaryFile(prefix = config._dataName+'_'+config._selection+'_', delete=False)

        for i, line in enumerate(infile.readlines()):
            if i != 6:
                exec_tmp.write(line)
            else:
                path = ''
                if config._inDir[:5] == '/pnfs':
                    path = 'dcap://cmsdca1.fnal.gov:24140/pnfs/fnal.gov/usr'+config._inDir[5:]
                else:
                    path = config._inDir

                for file in sourceFiles:
                    exec_tmp.write('echo '+repr(path+'/'+file)+' >> input.txt\n')
                    
        exec_tmp.seek(0)
        infile.close()
                
        return exec_tmp

    def MakeBatchConfig(self, config, count, exec_tmp, sourceFiles):
        '''Write batch configuration file'''
        batch_tmp = tempfile.NamedTemporaryFile(prefix = 'TEST', delete=False)

        batch_tmp.write('Arguments  = %s %s %s %s %s\n' % (self._current, self._outDir+'/'+config._selection, str(count+1), config._dataName, config._args))
        batch_tmp.write('Executable            = %s\n' % exec_tmp.name)
        batch_tmp.write('Should_Transfer_Files = YES\n')
        batch_tmp.write('WhenToTransferOutput  = ON_EXIT\n')
        batch_tmp.write('Transfer_Input_Files  = \n')
        batch_tmp.write('Universe              = vanilla\n')
        batch_tmp.write('Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000\n')
        if self._shortQueue:
            batch_tmp.write('+LENGTH               = "SHORT"\n')
        batch_tmp.write('Output                = reports/report_$(Cluster)_$(Process).stdout\n')
        batch_tmp.write('Error                 = reports/report_$(Cluster)_$(Process).stderr\n')
        batch_tmp.write('Log                   = reports/report_$(Cluster)_$(Process).log   \n')
        batch_tmp.write('Queue\n')
            
        batch_tmp.seek(0)
        return batch_tmp


    def CreateWorkingDir(self, codedir):
        '''Creating local working directory which will be used for copying all the files'''
        self.MakeDirectory(self._outDir, clear=False)
        code_copy_dir = self._outDir+'/'+codedir
        self.MakeDirectory(code_copy_dir, clear=True)

        for d in ["src", "plugins", "data", "higgs"]:
            self.MakeDirectory(code_copy_dir+'/'+d, clear=True)
            if   d == "data":
                os.system("cp -r "+self._outDir+"/../../"+d+"/* "+code_copy_dir+"/"+d)
            elif d =="plugins":
                self.MakeDirectory(code_copy_dir+'/'+d+'/lib', clear=False)
                os.system("cp "+self._outDir+"/../../"+d+"/lib/*.so "+code_copy_dir+"/"+d+"/lib")
                os.system("cp "+self._outDir+"/../../"+d+"/*.cc "+code_copy_dir+"/"+d)
                os.system("cp "+self._outDir+"/../../"+d+"/*.h "+code_copy_dir+"/"+d)
            elif d == "src":
                os.system("cp "+self._outDir+"/../../"+d+"/*.cc "+code_copy_dir+"/"+d)
                os.system("cp "+self._outDir+"/../../"+d+"/*.h "+code_copy_dir+"/"+d)
            else:
                os.system("cp "+self._outDir+"/../../"+d+"/*.C "+code_copy_dir+"/"+d)
                os.system("cp "+self._outDir+"/../../"+d+"/*.h "+code_copy_dir+"/"+d)
                os.system("cp "+self._outDir+"/../../"+d+"/*.py "+code_copy_dir+"/"+d)

        
    def SubmitToLPC(self):
        '''Submits batch jobs to lpc batch'''
        for cfg in self._configList:
            self.MakeDirectory(self._outDir+'/'+cfg._selection, clear=False)
            self.MakeDirectory(self._outDir+'/'+cfg._selection+'/printouts', clear=False)
            sourceFiles = self.SplitJobs(cfg._inDir, cfg._nJobs)
            for i, source in enumerate(sourceFiles):
                exec_tmp  = self.MakeExecutable(cfg, source, i)
                batch_tmp = self.MakeBatchConfig(cfg, i, exec_tmp, source)
                subprocess.call('condor_submit ' + batch_tmp.name , shell=True)
                exec_tmp.close()
                batch_tmp.close()
