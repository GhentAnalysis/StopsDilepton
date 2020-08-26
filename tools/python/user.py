import os

runOnGentT2 = True

if os.environ['USER'] in ['schoef', 'rschoefbeck', 'schoefbeck']:
    # Where you store cmg output
    cmg_directory      = "/scratch/rschoefbeck/cmgTuples/80X_1l_9"
    # Where postprocessed data goes 
    data_output_directory      = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
    #data_output_directory      = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
    postprocessing_output_directory = "/afs/hephy.at/data/cms07/nanoTuples/"
    # Where you store the data
    data_directory      = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
    # Where the plots go
    plot_directory      = "/afs/hephy.at/user/r/rschoefbeck/www/StopsDilepton"
    #plot_directory      = "/afs/cern.ch/work/s/schoef/www/"

    cache_dir = "/afs/hephy.at/data/cms01/stopsDilepton/signals/caches/"
    # Analysis result files
    analysis_results        =  '/afs/hephy.at/data/cms05/StopsDileptonLegacy/results/v9/'
    dpm_directory           = '/dpm/oeaw.ac.at/home/cms/store/user/schoef/'
    cern_proxy_certificate  = '/afs/cern.ch/user/s/schoef/private/.proxy'
    # directory with veto lists
    veto_lists = "/afs/hephy.at/data/rschoefbeck01/StopsDilepton/vetoLists/"
    # 715 release for limit calculation 
    combineReleaseLocation = '/afs/hephy.at/work/r/rschoefbeck/CMS/tmp/CMSSW_7_1_5/src'
    runOnGentT2 = False
    MVA_preprocessing_directory  = '/afs/hephy.at/work/g/gungersback/StopsDilepton/MVA_preprocessing'
    MVA_model_directory          = '/afs/hephy.at/work/g/gungersback/StopsDilepton/MVA_models'

if os.environ['USER'] in ['robert.schoefbeck']:
    #data_output_directory      = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
    postprocessing_output_directory = "/mnt/hephy/cms/robert.schoefbeck/StopsDileptonLegacy/nanoTuples"
    # Where the plots go
    plot_directory      = "/mnt/hephy/cms/robert.schoefbeck/www/StopsDileptonLegacy/plots"
    # Analysis result files
    analysis_results        =  '/mnt/hephy/cms/robert.schoefbeck/StopsDileptonLegacy/results/v9'
    cern_proxy_certificate  = '/users/robert.schoefbeck/.private/.proxy'
    cache_dir               = "/mnt/hephy/cms/robert.schoefbeck/StopsDileptonLegacy/caches"
    # directory with veto lists
    #veto_lists = "/afs/hephy.at/data/rschoefbeck01/StopsDilepton/vetoLists/"
    # 715 release for limit calculation 
    #combineReleaseLocation = '/afs/hephy.at/work/r/rschoefbeck/CMS/tmp/CMSSW_7_1_5/src'
    runOnGentT2 = False
    #MVA_preprocessing_directory  = '/afs/hephy.at/work/g/gungersback/StopsDilepton/MVA_preprocessing'
    #MVA_model_directory          = '/afs/hephy.at/work/g/gungersback/StopsDilepton/MVA_models'

if os.environ['USER'] in ['phussain']:
    # Where you store cmg output
    cmg_directory      = "/scratch/rschoefbeck/cmgTuples/80X_1l_9"
    # Where postprocessed data goes 
    data_output_directory      = "/afs/hephy.at/data/cms02/nanoTuples/"
    # Where you store the data
    data_directory      = "/afs/hephy.at/data/dspitzbart03/nanoTuples/"
    # Where the plots go
    plot_directory      = "/afs/hephy.at/user/p/phussain/www/stopsDilepton/"
    private_results_directory     = "/afs/hephy.at/data/cms01/"
    cache_dir = "/afs/hephy.at/data/cms01/stopsDilepton/signals/caches/"
    #plot_directory      = "/afs/cern.ch/work/d/dspitzba/www/stopsDilepton/"
    # Analysis result files
    #analysis_results        = '/afs/hephy.at/data/dspitzbart02/StopsDileptonLegacy/results/v1/' #Path to analysis results
    #dpm_directory           = '/dpm/oeaw.ac.at/home/cms/store/user/dspitzba/'
    #cern_proxy_certificate  = '/afs/cern.ch/user/d/dspitzba/private/.proxy'
    #postprocessing_output_directory = "/afs/hephy.at/data/dspitzbart03/nanoTuples/"
    postprocessing_output_directory = "/afs/hephy.at/data/cms02/nanoTuples/"
    # directory with veto lists
    analysis_results        = '/afs/hephy.at/data/cms05/StopsDileptonLegacy/results/v7/'
    #analysis_results        = '/afs/hephy.at/work/p/phussain/StopsDileptonLegacy/results/v4/'
    #analysis_results        = '/afs/hephy.at/work/p/phussain/StopsDileptonLegacy/results/vpseudomet/'
    veto_lists = "/afs/hephy.at/data/rschoefbeck01/StopsDilepton/vetoLists/"
    MVA_preprocessing_directory  = '/afs/hephy.at/work/g/gungersback/StopsDilepton/MVA_preprocessing'
    MVA_model_directory          = '/afs/hephy.at/work/g/gungersback/StopsDilepton/MVA_models'
    # 715 release for limit calculation 
    combineReleaseLocation = '/afs/hephy.at/work/d/dspitzbart/stops/CMSSW_7_4_7/src'
    runOnGentT2 = False

if os.environ['USER'] in ['dspitzbart', 'dspitzba']:
    # Where postprocessed data goes 
    data_output_directory      = "/afs/hephy.at/data/cms02/nanoTuples/"
    # Where you store the data
    data_directory      = "/afs/hephy.at/data/cms02/nanoTuples/"
    # Where the plots go
    plot_directory      = "/afs/hephy.at/user/d/dspitzbart/www/stopsDileptonLegacy/"
    # Analysis result files
    analysis_results        = '/afs/hephy.at/data/cms05/StopsDileptonLegacy/results/v8/' #Path to analysis results
    #analysis_results        = '/afs/hephy.at/data/cms05/StopsDileptonLegacy/results/v6/' #Path to analysis results
    dpm_directory           = '/dpm/oeaw.ac.at/home/cms/store/user/dspitzba/'
    cern_proxy_certificate  = '/afs/cern.ch/user/d/dspitzba/private/.proxy'
    #postprocessing_output_directory = "/afs/hephy.at/data/dspitzbart03/nanoTuples/"
    cache_dir = "/afs/hephy.at/data/cms01/stopsDilepton/signals/caches/"
    postprocessing_output_directory = "/afs/hephy.at/data/cms03/nanoTuples/"
    if 'cern' in os.getenv("HOSTNAME"):
      postprocessing_output_directory = "/eos/home-d/dspitzba/nanoTuples/"
      analysis_results        = '/afs/cern.ch/work/d/dspitzba/StopsDilepton/results/80X_v30'
      plot_directory      = "/eos/home-d/dspitzba/plots/"

if os.environ['USER'] in ['mdoppler']:
    # Where you store cmg output
    cmg_directory      = "/scratch/rschoefbeck/cmgTuples/80X_1l_9"
    # Where postprocessed data goes 
    data_output_directory      = "/afs/hephy.at/data/cms05/nanoTuples/"
    postprocessing_output_directory = "/afs/hephy.at/data/cms05/nanoTuples/"
    # Where you store the data
    data_directory      = "/afs/hephy.at/data/dspitzbart01/nanoTuples/"
    # Where the plots go
    plot_directory      = "/afs/hephy.at/user/m/mdoppler/www/stopsDileptonLegacy/"
    #plot_directory      = "/afs/cern.ch/work/d/dspitzba/www/stopsDilepton/"
    # Analysis result files
    analysis_results        = '/afs/hephy.at/data/cms05/StopsDileptonLegacy/results/v7/'
    #dpm_directory           = '/dpm/oeaw.ac.at/home/cms/store/user/dspitzba/'
    #cern_proxy_certificate  = '/afs/cern.ch/user/d/dspitzba/private/.proxy'
    # directory with veto lists
    veto_lists = "/afs/hephy.at/data/rschoefbeck01/StopsDilepton/vetoLists/"
    cache_dir = "/afs/hephy.at/data/cms01/stopsDilepton/signals/caches/"
    MVA_preprocessing_directory  = '/afs/hephy.at/work/g/gungersback/StopsDilepton/MVA_preprocessing'
    MVA_model_directory          = '/afs/hephy.at/work/g/gungersback/StopsDilepton/MVA_models'
    # 715 release for limit calculation 
    combineReleaseLocation = '/afs/hephy.at/work/d/dfreisinger/cms/CMSSW_7_4_7/src'
    if 'cern' in os.getenv("HOSTNAME"):
      analysis_results        = '/afs/cern.ch/work/d/dspitzba/StopsDilepton/results/80X_v30'
      combineReleaseLocation = '/afs/cern.ch/work/d/dspitzba/CMS/stop/CMSSW_7_4_7/src'
    runOnGentT2 = False

if os.environ['USER'] in ['llechner']:
    # Where postprocessed data goes 
    data_output_directory      = "/afs/hephy.at/data/cms01/nanoTuples/"
    # Where you store the data
    data_directory      = "/afs/hephy.at/data/cms01/nanoTuples/"
    # Where the plots go
    plot_directory      = "/afs/hephy.at/user/l/llechner/www/stopsDileptonLegacy/"
    # Analysis result files
    analysis_results        = '/afs/hephy.at/data/llechner01/StopsDileptonLegacy/results/v1/' #Path to analysis results
    dpm_directory           = '/dpm/oeaw.ac.at/home/cms/store/user/llechner/'
    cache_dir = "/afs/hephy.at/data/cms01/stopsDilepton/signals/caches/"
    cern_proxy_certificate  = '/afs/cern.ch/user/l/llechner/private/.proxy'
    #postprocessing_output_directory = "/afs/hephy.at/data/dspitzbart03/nanoTuples/"
    postprocessing_output_directory = "/afs/hephy.at/data/cms01/nanoTuples/"
    if 'cern' in os.getenv("HOSTNAME"):
      postprocessing_output_directory = "/eos/home-d/llechner/nanoTuples/"
      analysis_results        = '/afs/cern.ch/work/l/llechner/StopsDilepton/results/80X_v30'
      plot_directory      = "/eos/home-d/llechner/plots/"

