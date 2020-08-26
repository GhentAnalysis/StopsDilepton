import ROOT, os
ROOT.gROOT.SetBatch(True)
from RootTools.core.standard             import *
from Analysis.Tools.helpers import writeObjToFile

TTLep_pow_16 = Sample.fromDirectory("TTLep_pow_16", "/afs/hephy.at/data/cms07/nanoTuples/stops_2016_nano_v0p24/dilep/TTLep_pow") 
TTLep_pow_17 = Sample.fromDirectory("TTLep_pow_17_ext", "/afs/hephy.at/data/cms06/nanoTuples/stops_2017_nano_v0p23/dilep/TTLep_pow_ext")
TTLep_pow_18 = Sample.fromDirectory("TTLep_pow_18", "/afs/hephy.at/data/cms07/nanoTuples/stops_2018_nano_v0p24/dilep/TTLep_pow")

mc = [TTLep_pow_16, TTLep_pow_17, TTLep_pow_18]
#mc = [TTLep_pow_17]

for sample in mc:
    print "At", sample.name
    #sample.reduceFiles(to=1)
    print sample.chain.GetEntries()
    sample.chain.GetEntry(0)
    res = sample.chain.GetTree().CopyTree("dl_mt2ll>100","")
    print res.GetEntries()
    writeObjToFile( sample.name+'_mt2ll100.root', res )
