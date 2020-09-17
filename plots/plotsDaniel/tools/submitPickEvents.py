import os
import subprocess

eventList = []
with open("events.txt","r") as f:
    while True:
        l = f.readline()
        if l =='': break
        eventList.append(l.replace('\n',''))

#runsAndEras = {\
#    '2016B':    {'range': (272007,   275376), 'dataset':'/DoubleMuon/Run2016B-03Feb2017_ver2-v2/MINIAOD'},
#    '2016C':    {'range': (275657,   276283), 'dataset':'/DoubleMuon/Run2016C-03Feb2017-v1/MINIAOD'},
#    '2016D':    {'range': (276315,   276811), 'dataset':'/DoubleMuon/Run2016D-03Feb2017-v1/MINIAOD'},
#    '2016E':    {'range': (276831,   277420), 'dataset':'/DoubleMuon/Run2016E-03Feb2017-v1/MINIAOD'},
#    '2016F':    {'range': (277772,   278808), 'dataset':'/DoubleMuon/Run2016F-03Feb2017-v1/MINIAOD'},
#    '2016G':    {'range': (278820,   280385), 'dataset':'/DoubleMuon/Run2016G-03Feb2017-v1/MINIAOD'},
#    '2016Hv2':  {'range': (280919,   284037), 'dataset':'/DoubleMuon/Run2016H-03Feb2017_ver2-v1/MINIAOD'},   
#    '2016Hv3':  {'range': (284037,   284044), 'dataset':'/DoubleMuon/Run2016H-03Feb2017_ver3-v1/MINIAOD'},
#    }

#runsAndEras = {\
#    '2016B':    {'range': (272007,   275376), 'dataset':'/SingleElectron/Run2016B-03Feb2017_ver2-v2/MINIAOD'},
#    '2016C':    {'range': (275657,   276283), 'dataset':'/SingleElectron/Run2016C-03Feb2017-v1/MINIAOD'},
#    '2016D':    {'range': (276315,   276811), 'dataset':'/SingleElectron/Run2016D-03Feb2017-v1/MINIAOD'},
#    '2016E':    {'range': (276831,   277420), 'dataset':'/SingleElectron/Run2016E-03Feb2017-v1/MINIAOD'},
#    '2016F':    {'range': (277772,   278808), 'dataset':'/SingleElectron/Run2016F-03Feb2017-v1/MINIAOD'},
#    '2016G':    {'range': (278820,   280385), 'dataset':'/SingleElectron/Run2016G-03Feb2017-v1/MINIAOD'},
#    '2016Hv2':  {'range': (280919,   284037), 'dataset':'/SingleElectron/Run2016H-03Feb2017_ver2-v1/MINIAOD'},
#    '2016Hv3':  {'range': (284037,   284044), 'dataset':'/SingleElectron/Run2016H-03Feb2017_ver3-v1/MINIAOD'},
#    }

runsAndEras = {\
    '2016B':    {'range': (272007,   275376), 'dataset':'/DoubleMuon/Run2016B-07Aug17_ver2-v2/MINIAOD'},
    '2016C':    {'range': (275657,   276283), 'dataset':'/DoubleMuon/Run2016C-07Aug17-v1/MINIAOD'},
    '2016D':    {'range': (276315,   276811), 'dataset':'/DoubleMuon/Run2016D-07Aug17-v1/MINIAOD'},
    '2016E':    {'range': (276831,   277420), 'dataset':'/DoubleMuon/Run2016E-07Aug17-v1/MINIAOD'},
    '2016F':    {'range': (277772,   278808), 'dataset':'/DoubleMuon/Run2016F-07Aug17-v1/MINIAOD'},
    '2016G':    {'range': (278820,   280385), 'dataset':'/DoubleMuon/Run2016G-07Aug17-v1/MINIAOD'},
    '2016H':    {'range': (280919,   284037), 'dataset':'/DoubleMuon/Run2016H-07Aug17-v1/MINIAOD'},   
    }

runsAndEras = {\
    '2017B':    {'range': (297020,   299329), 'dataset':'/DoubleMuon/Run2017B-17Nov2017-v1/MINIAOD'},
    '2017C':    {'range': (299337,   302029), 'dataset':'/DoubleMuon/Run2017C-17Nov2017-v1/MINIAOD'},
    '2017D':    {'range': (302030,   303434), 'dataset':'/DoubleMuon/Run2017D-17Nov2017-v1/MINIAOD'},
    '2017E':    {'range': (303435,   304826), 'dataset':'/DoubleMuon/Run2017E-17Nov2017-v1/MINIAOD'},
    '2017F':    {'range': (304911,   306462), 'dataset':'/DoubleMuon/Run2017F-17Nov2017-v1/MINIAOD'},
    }

runsAndEras = {\
    '2018A':    {'range': (315252, 316995  ), 'dataset':''},
    '2018B':    {'range': (316998, 319312  ), 'dataset':''},
    '2018C':    {'range': (319313, 320393  ), 'dataset':''},
    '2018D':    {'range': (320394, 325273  ), 'dataset':'/DoubleMuon/Run2018D-PromptReco-v2/MINIAOD'},
    }

#runsAndEras = {\
#    '2016': {'range': (0, 1e10), 'dataset':'/TTTo2L2Nu_noSC_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'},
#    }

for i,ev in enumerate(eventList):
    run_str = ev.replace(':','_')
    run = int(ev.split(':')[0])
    print run
    dataset = False
    for era in runsAndEras:
        if run >= runsAndEras[era]['range'][0] and run <= runsAndEras[era]['range'][1]:
            dataset = runsAndEras[era]['dataset']
    if dataset:
        print "Event %s should be in datset %s"%(ev, dataset)
    else: print ":("
    oFile = "picks/event_%s"%i
    cmd = ["edmPickEvents.py", dataset, ev, "--output", oFile, "--runInteractive"]
    print cmd
    print ' '.join(cmd)
    subprocess.call(cmd)

