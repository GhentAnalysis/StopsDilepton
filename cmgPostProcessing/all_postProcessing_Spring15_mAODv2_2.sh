#!/bin/sh

python cmgPostProcessing.py --skim=$1 $2  --samples=TTJets_SingleLeptonFromTbar,TTJets_SingleLeptonFromTbar_ext &
python cmgPostProcessing.py --skim=$1 $2  --samples=TTJets_SingleLeptonFromT,TTJets_SingleLeptonFromT_ext &
python cmgPostProcessing.py --skim=$1 $2  --samples=TTLep_pow &
