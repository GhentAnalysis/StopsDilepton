# Signal region

#python systematicVariation.py --era Run2016 --plot_directory v0p23 --reweightPU Central --selection lepSel-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --signal T2tt #--dpm #--overwrite
#python systematicVariation.py --era Run2017 --plot_directory v0p23 --reweightPU Central --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --signal T2tt --add #--dpm #--overwrite
#python systematicVariation.py --era Run2018 --plot_directory v0p23 --reweightPU VUp --selection lepSel-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --signal T2tt --add #--dpm #--overwrite

python systematicVariation.py --era Run2016 --plot_directory v0p23 --reweightPU Central --selection lepSel-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-mt2ll100 --signal T2tt --add  
#python systematicVariation.py --era Run2017 --plot_directory v0p23 --reweightPU Central --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-mt2ll100 --signal T2tt --add #--dpm #--overwrite
#python systematicVariation.py --era Run2018 --plot_directory v0p23 --reweightPU VUp --selection lepSel-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-mt2ll100 --signal T2tt --add #--dpm #--overwrite

#python systematicVariation.py --dpm --add --era Run2016 --signal T2tt --plot_directory v0p23 --reweightPU Central --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-mt2ll100 --newMetSigPlots #--overwrite
#python systematicVariation.py --dpm --era Run2017 --signal T2tt --plot_directory v0p23 --reweightPU Central --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-mt2ll100 --add --newMetSigPlots #--overwrite
#python systematicVariation.py --dpm --era Run2018 --signal T2tt --plot_directory v0p23 --reweightPU VUp --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-mt2ll100 --newMetSigPlots --add #--overwrite
#
#python systematicVariation.py --dpm --era Run2016 --signal T2tt --plot_directory v0p23 --reweightPU Central --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --newMetSigPlots --add #--overwrite
#python systematicVariation.py --dpm --era Run2017 --signal T2tt --plot_directory v0p23 --reweightPU Central --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add --newMetSigPlots #--overwrite
#python systematicVariation.py --dpm --era Run2018 --signal T2tt --plot_directory v0p23 --reweightPU VUp --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --newMetSigPlots --add #--overwrite

#
#
## CR plots 02.28.2020
#python systematicVariation.py --era Run2016 --plot_directory v0p22 --reweightPU Central --selection lepSel-POGMetSig12-njet01-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-allZ
#python systematicVariation.py --era Run2017 --plot_directory v0p22 --reweightPU Central --selection lepSel-badEEJetVeto-POGMetSig12-njet01-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-allZ --add
#python systematicVariation.py --era Run2018 --plot_directory v0p22 --reweightPU VUp     --selection lepSel-POGMetSig12-njet01-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-allZ --add
#
## DY Check 02.24.2020
#python systematicVariation.py --era Run2016 --plot_directory v0p22 --reweightPU Central --normalize --selection lepSel-njet2p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ 
#python systematicVariation.py --era Run2017 --plot_directory v0p22 --reweightPU Central --normalize --selection lepSel-badEEJetVeto-njet2p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ --add 
#python systematicVariation.py --era Run2018 --plot_directory v0p22 --reweightPU VUp     --normalize --selection lepSel-njet2p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ --add   
#
#python systematicVariation.py --era Run2016 --plot_directory v0p22 --reweightPU Central --normalize --selection lepSel-POGMetSigTo12-njet2p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ 
#python systematicVariation.py --era Run2017 --plot_directory v0p22 --reweightPU Central --normalize --selection lepSel-POGMetSigTo12-badEEJetVeto-njet2p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ --add 
#python systematicVariation.py --era Run2018 --plot_directory v0p22 --reweightPU VUp     --normalize --selection lepSel-POGMetSigTo12-njet2p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ --add   

## DY chi-squared - Check 01.21.2020
#python systematicVariation.py --era Run2016 --plot_directory v0p23 --reweightPU Central  --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ  --normalize #--dpm --overwrite  --small
#python systematicVariation.py --era Run2017 --plot_directory v0p23 --reweightPU Central  --selection lepSel-badEEJetVeto-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ --add  --normalize #--dpm --overwrite  --small
#python systematicVariation.py --era Run2018 --plot_directory v0p23 --reweightPU VUp  --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ --add --normalize #--dpm --overwrite  

## DY CR tail study 01_29_2020
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --beta noScaling --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --beta noScaling --selection lepSel-badEEJetVeto-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --beta noScaling --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ
#
#
## Sidebands
#
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --beta overflowBin #--dpm #--overwrite 
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add --beta overflowBin #--dpm #--overwrite
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add --beta overflowBin #--dpm #--overwrite
#
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --normalize --selection lepSel-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-mt2ll100 --signal T2tt #--dpm --overwrite
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --normalize --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-mt2ll100 --signal T2tt --add #--dpm --overwrite
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --normalize --selection lepSel-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-mt2ll100 --signal T2tt --add #--dpm --overwrite
#
## Sidebands
#
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ #--dpm #--overwrite 
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add #--dpm #--overwrite
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --normalize --selection lepSel-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-mt2ll100 --signal T2tt --add #--dpm --overwrite
#
## Sidebands
#
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ #--dpm #--overwrite 
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add #--dpm #--overwrite
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add #--dpm #--overwrite
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add --beta overflowBin #--dpm #--overwrite
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add --beta overflowBin #--dpm #--overwrite
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add --beta overflowBin #--dpm #--overwrite
#
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ #--dpm #--overwrite
#
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add #--dpm #--overwrite
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add #--dpm #--overwrite
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add #--dpm #--overwrite
#
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ-mt2ll100 --add --normalize #--dpm #--overwrite
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-mt2ll100 --add --normalize #--dpm #--overwrite
#
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-badEEJetVeto-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add --beta overflowBin #--dpm #--overwrite
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-badEEJetVeto-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add --beta overflowBin #--dpm #--overwrite
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add --beta overflowBin #--dpm #--overwrite
#
#
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-badEEJetVeto-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add #--dpm #--overwrite
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-badEEJetVeto-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add #--dpm #--overwrite
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add #--dpm #--overwrite
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-badEEJetVeto-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add --beta overflowBin #--dpm #--overwrite
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-badEEJetVeto-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add --beta overflowBin #--dpm #--overwrite
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add --beta overflowBin #--dpm #--overwrite
#
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ --add #--dpm #--overwrite
#
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-badEEJetVeto-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add #--dpm #--overwrite
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-badEEJetVeto-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add #--dpm #--overwrite
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --scaling mc --variation_scaling --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add #--dpm #--overwrite
#
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ-mt2ll100 --add --normalize #--dpm #--overwrite
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-mt2ll100 --add --normalize #--dpm #--overwrite
#
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --scaling mc --variation_scaling --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add --beta overflowBin #--dpm #--overwrite
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --scaling mc --variation_scaling --selection lepSel-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add --beta overflowBin #--dpm #--overwrite
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --scaling mc --variation_scaling --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add --beta overflowBin #--dpm #--overwrite
#
#
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --scaling mc --variation_scaling --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add #--dpm #--overwrite
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --scaling mc --variation_scaling --selection lepSel-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add #--dpm #--overwrite
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --scaling mc --variation_scaling --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ --add #--dpm #--overwrite
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --scaling mc --variation_scaling --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add --beta overflowBin #--dpm #--overwrite
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --scaling mc --variation_scaling --selection lepSel-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add --beta overflowBin #--dpm #--overwrite
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --scaling mc --variation_scaling --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add --beta overflowBin #--dpm #--overwrite
#
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --scaling mc --variation_scaling --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ --add #--dpm #--overwrite
#
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --scaling mc --variation_scaling --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add #--dpm #--overwrite
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --scaling mc --variation_scaling --selection lepSel-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add #--dpm #--overwrite
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --scaling mc --variation_scaling --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --add #--dpm #--overwrite
#
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ-mt2ll100 --add --normalize #--dpm #--overwrite
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-mt2ll100 --add --normalize #--dpm #--overwrite
#
#
#
######################################################################################################################################################################################
#
## reproduced with signal 02.24.2020
### fig:metSig-presel-comb (EventSelection.tex)
#python systematicVariation.py --signal T2tt --era Run2016 --plot_directory v0p22 --reweightPU Central --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --beta noScaling
#python systematicVariation.py --signal T2tt --era Run2017 --plot_directory v0p22 --reweightPU Central --selection lepSel-badEEJetVeto-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --beta noScaling --add
#python systematicVariation.py --signal T2tt --era Run2018 --plot_directory v0p22 --reweightPU VUp     --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1 --beta noScaling --add


#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --normalize --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ
#
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --normalize --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ #--dpm #--overwrite
#
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --normalize --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --normalize --selection lepSel-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1
#python systematicVariation.py --era Run2016 --plot_directory v0p19 --reweightPU Central --normalize --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1
#
##2017
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --normalize --selection lepSel-badEEJetVeto-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --normalize --selection lepSel-badEEJetVeto-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --normalize --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ
#
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --normalize --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ
#
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --normalize --selection lepSel-badEEJetVeto-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --normalize --selection lepSel-badEEJetVeto-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1
#python systematicVariation.py --era Run2017 --plot_directory v0p19 --reweightPU Central --normalize --selection lepSel-badEEJetVeto-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1
#
##2018
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --normalize --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --normalize --selection lepSel-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --normalize --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1-onZ
#
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --normalize --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-onZ
#
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --normalize --selection lepSel-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --normalize --selection lepSel-POGMetSigTo12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1
#python systematicVariation.py --era Run2018 --plot_directory v0p19 --reweightPU VUp --normalize --selection lepSel-POGMetSig12-njet2p-btag0-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1



