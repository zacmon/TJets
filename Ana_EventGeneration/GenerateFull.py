import os
import time
import threading
import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--type",      type=str,  default="ww",              help="Type of the physical process you wish to generate.  Possible values are [\"jj\"=parton jets, \"ww\"=boosted W bosons], \"tt\"=boosted top quarks]")
parser.add_argument("--nevents",   type=int,  default=1000,              help="The number of events to create in a single job")
parser.add_argument("--jobstart",  type=int,  default=0,                 help="When you want to generate multiple jobs, this is the lower job number in the list of [jobStart,jobEnd] jobs.")
parser.add_argument("--jobend",    type=int,  default=1,                 help="When you want to generate multiple jobs, this is the upper job number in the list of [jobStart,jobEnd] jobs.")
parser.add_argument("--recompile", default=False,  action="store_true",  help="Flag that will recompile the program before running")
parser.add_argument("--pythia",    default=False,  action="store_true",  help="Flag that will run the generation of Pythia events")
parser.add_argument("--ntuple",    default=False,  action="store_true",  help="Flag that will run the NTupler over pre-made Pythia particle ntuples")
parser.add_argument("--debug",     default=False,  action="store_true",  help="Flag that will turn on the debugging printout")
args = parser.parse_args()   

arg_process   = args.type
arg_nevents   = args.nevents
arg_jobstart  = args.jobstart
arg_jobend    = args.jobend
arg_recompile = args.recompile
arg_pythia    = args.pythia
arg_ntuple    = args.ntuple
arg_debug = "false"
if args.debug:
    arg_debug  = "true"

print "Arguments     :"
print "arg_process   :",arg_process 
print "arg_nevents   :",arg_nevents 
print "arg_jobstart  :",arg_jobstart
print "arg_jobend    :",arg_jobend  
print "arg_recompile :",arg_recompile
print "arg_pythia    :",arg_pythia   
print "arg_ntuple    :",arg_ntuple   
print "arg_debug     :",arg_debug   

#making output directory for storage
outputdirPythia="GenPythia/"
outputdirPythia+=time.strftime("%Y%m%d")
outputdirPythia+="/"
os.system("mkdir -p "+outputdirPythia)

#making output directory for storage
outputdirNtuple="GenNTuple/"
outputdirNtuple+=time.strftime("%Y%m%d")
outputdirNtuple+="/"
os.system("mkdir -p "+outputdirNtuple)


#######################
#EVENT GENERATION
#######################
if arg_pythia:
    #compiling code freshly
    if arg_recompile:
        os.system("source compile_pythia.sh")

    #loop over and produce jobs from jobStart to jobEnd
    for iJob in range(arg_jobstart,arg_jobend):

        pythiaFilename  = outputdirPythia+"pythia_"+arg_process+"_"+str(iJob)+".root"
    
        runcommand = "./MakeNTupleFromPythia.exe   "+arg_process+"   "+pythiaFilename+"   "+str(arg_nevents)+"   "+arg_debug

        print "Running Command : ",runcommand

        os.system(runcommand)
    
        print "Finished : ",runcommand
    
    
#######################
#JET NTUPLING
#######################
if arg_ntuple:
    #compiling code freshly
    if arg_recompile:
        os.system("source compile_ntupler.sh")

    #loop over and produce jobs from jobStart to jobEnd
    for iJob in range(arg_jobstart,arg_jobend):

        pythiaFilename  = outputdirPythia+"pythia_"+arg_process+"_"+str(iJob)+".root"
        ntupleFilename  = outputdirNtuple+"ntuple_"+arg_process+"_"+str(iJob)+".root"
    
        runcommand = "./NTupler.exe   "+arg_process+"   "+pythiaFilename+"   "+ntupleFilename+"   "+arg_debug

        print "Running Command : ",runcommand

        os.system(runcommand)
    
        print "Finished : ",runcommand

