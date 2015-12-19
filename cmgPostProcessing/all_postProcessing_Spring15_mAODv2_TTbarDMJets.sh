#!/bin/sh
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi1_Mphi10" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi1_Mphi20" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi1_Mphi50" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi1_Mphi100" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi1_Mphi200" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi1_Mphi300" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi1_Mphi500" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi10_Mphi10" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi10_Mphi50" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi10_Mphi100" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi50_Mphi50" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi50_Mphi200" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi50_Mphi300" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi150_Mphi200" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi150_Mphi500" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi150_Mphi1000" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_pseudoscalar_Mchi500_Mphi500" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_scalar_Mchi1_Mphi10" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_scalar_Mchi1_Mphi20" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_scalar_Mchi1_Mphi100" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_scalar_Mchi1_Mphi200" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_scalar_Mchi1_Mphi300" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_scalar_Mchi1_Mphi500" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_scalar_Mchi1_Mphi1000" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_scalar_Mchi10_Mphi10" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_scalar_Mchi10_Mphi100" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_scalar_Mchi50_Mphi50" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_scalar_Mchi50_Mphi200" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_scalar_Mchi50_Mphi300" &
nohup krenew -t -K 10 -- bash -c "python cmgPostProcessing.py --skim=$1 $2 --samples=TTbarDMJets_scalar_Mchi500_Mphi500" &
##python cmgPostProcessing.py --skim=dilepTiny  --samples=TTbarDMJets_scalar_Mchi1_Mphi50 &
##python cmgPostProcessing.py --skim=dilepTiny  --samples=TTbarDMJets_scalar_Mchi10_Mphi50 &
##python cmgPostProcessing.py --skim=dilepTiny  --samples=TTbarDMJets_scalar_Mchi150_Mphi200 &
