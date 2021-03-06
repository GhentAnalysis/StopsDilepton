#Standard imports
import numpy as np
import math

from optparse import OptionParser
parser = OptionParser()
#parser.add_option("--year",       dest="year",                  default=2016, type="int",    action="store",      help="Which year?")
parser.add_option('--logLevel',   dest="logLevel",              default='INFO',              action='store',      help="log level?", choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'])

(options, args) = parser.parse_args()

# Logging
import Analysis.Tools.logger as logger
logger  = logger.get_logger(options.logLevel, logFile = None)
import StopsDilepton.tools.logger as logger
logger = logger.get_logger(options.logLevel, logFile = None )
import Analysis.Tools.logger as logger_an
logger_an = logger_an.get_logger(options.logLevel, logFile = None )


from Analysis.Tools.cardFileWriter.CombineResults import CombineResults
cardFile = "/afs/hephy.at/data/cms05/StopsDileptonLegacy/results/v7/COMBINED/fitAll/cardFiles/T2tt/observed/T2tt_800_100.txt"
plotDirectory = "/afs/hephy.at/user/r/rschoefbeck/www/StopsDilepton/postFit"

Results = CombineResults( cardFile=cardFile, plotDirectory=plotDirectory, year=0, bkgOnly=True, isSearch=True )
yields = Results.getEstimates( postFit=True )

#L1p=Results.getNuisanceYields("L1prefire", postFit=True)
#SFb   =Results.getNuisanceYields("SFb_2016", postFit=True)

systematics = [ 
    ("Jet energy scale", "JEC", ["JEC_2018", "JEC_2016", "JEC_2017"]), 
    ("Integrated luminosity", "Lumi", ["Lumi_2016", "Lumi_2017", "Lumi_2018"]),
    ("Pileup modeling", "PU", ["PU_2016", "PU_2017", "PU_2018"]),
    ("Jet energy resolution", "JER", ["JER_2018", "JER_2017",  "JER_2016"]),
    ("Modeling of unclust. en.", "unclEn", ["unclEn_2016", "unclEn_2017", "unclEn_2018"]), 
    ("Trigger efficiency", "trigger", ["trigger_2018", "trigger_2016", "trigger_2017"]),
    ("b tagging light flavor", "SFl", ["SFl_2018", "SFl_2017", "SFl_2016"]),
    ("b tagging heavy flavor", "SFb", ["SFb_2017", "SFb_2016", "SFb_2018"]),
    ("Lepton scale factors", "leptonSF", ["leptonSF"]),
    ("L1 prefire correction", "L1prefire", ["L1prefire"]),
    ("0 missing hit scale factor", "leptonHit0SF", ["leptonHit0SF"]),
    ("Impact parameter scale factor", "leptonSIP3DSF", ["leptonSIP3DSF"]),
    ("PDF choice", "PDF", ["PDF"]),
    ("$p_{T}(\\textrm{top})$", "topPt", ["topPt"]),
    ("t#bar{t} cross section", "topXSec", ["topXSec"]),
    ("t#bar{t}Z background", "ttZ_SR", ["ttZ_SR"]),
    ("fake/non-prompt leptons", "topFakes", ["topFakes"]),
    ("Drell-Yan background", "DY_SR", ["DY_SR"]),
    ("non-gaussian jet mismeasurements", "topNonGauss", ["topNonGauss"]),
    ("multiboson background", "MB_SR", ["MB_SR"]),
    ("$\mu_{R}$ and $\mu_{F}$ choice $(t\bar{t})$", "scaleTT", ["scaleTT"]),
    ("rare background", "other", ["other"]),
    ("$\mu_{R}$ and $\mu_{F}$ choice $(t\bar{t}Z)$", "scaleTTZ", ["scaleTTZ"]),
    ("Drell-Yan tail", "DY_hMT2ll", ["DY_hMT2ll"]), 
]


#"btagFS"
#"leptonFS"

for postFit in [True, False]:

    sys={}
    variations      = {}
    table_strings   = []
    total_scale     = {'typ':0, 'max':0}
    total_lepton    = {'typ':0, 'max':0}
    total_DY        = {'typ':0, 'max':0}
    for tex_name, name, sys_names in systematics:
        variations[name] = []
        max_var  = 0
        max_year = 0
        max_bin  = ""
        for sys_name in sys_names:
            #sys[sys_name] = Results.getNuisanceYields(sys_name, postFit=postFit, addSignal=False)
            sys[sys_name] = Results.getNuisanceYields(sys_name, postFit=postFit)
            for year in [2016, 2017, 2018]:
                for bin in range(20, 46):
                    binName = "dc_%i_Bin%i"%( year, bin )
                    variation = abs(0.5*(sys[sys_name][binName]['relUp']-sys[sys_name][binName]['relDown']))
                    if variation>0:
                        variations[name].append( (variation, binName))
                    if variation>max_var:
                        max_var  = variation   
                        max_year = year
                        max_bin  = binName
        all_variations = [ variation[0] for variation in variations[name] ]
        logger.info( "Systematic %20s. variation quantiles 50 / 68 / 90 / max %4.3f %4.3f %4.3f %4.3f max: %i %s %4.3f", name, \
            np.percentile(all_variations, 50) if len(all_variations)>0 else 0, 
            np.percentile(all_variations, 68) if len(all_variations)>0 else 0, 
            np.percentile(all_variations, 90) if len(all_variations)>0 else 0, 
            max(all_variations) if len(all_variations)>0 else 0,
            max_year, max_bin, max_var
                )
        if name.count('scale'):
            logger.info("Adding to scale uncertainty")
            total_scale['typ'] += (np.percentile(all_variations, 90)**2 if len(all_variations)>0 else 0)
            total_scale['max'] += max(all_variations)**2
        elif name.count('lepton'):
            logger.info("Adding to lepton uncertainty")
            total_lepton['typ'] += (np.percentile(all_variations, 90)**2 if len(all_variations)>0 else 0)
            total_lepton['max'] += max(all_variations)**2
        elif name.count('DY'):
            logger.info("Adding to DY uncertainty")
            total_DY['typ'] += (np.percentile(all_variations, 90)**2 if len(all_variations)>0 else 0)
            total_DY['max'] += max(all_variations)**2
        else:
            #table_strings.append( "%s & $ %5.1f $\\\\"%( tex_name, 100*(np.percentile(all_variations, 90) if len(all_variations)>0 else 0)) )
            table_strings.append( "{:35} & $ {:<8.1f} $ & $ {:<8.1f} $ ".format( tex_name, 100*(np.percentile(all_variations, 90) if len(all_variations)>0 else 0), 100*max(all_variations)) )
    table_strings.append( "{:35} & $ {:<8.1f} $ & $ {:<8.1f} $ ".format( '#mu_{R} and #mu_{F} choice',  100*math.sqrt(total_scale['typ']),  100*math.sqrt(total_scale['max'])) )
    table_strings.append( "{:35} & $ {:<8.1f} $ & $ {:<8.1f} $ ".format( 'Lepton scale factors',        100*math.sqrt(total_lepton['typ']), 100*math.sqrt(total_lepton['max'])) )
    table_strings.append( "{:35} & $ {:<8.1f} $ & $ {:<8.1f} $ ".format( 'Drell-Yan background',        100*math.sqrt(total_DY['typ']),     100*math.sqrt(total_DY['max'])) )

    #unc    = Results.getUncertainties( postFit=True)
    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in xrange(0, len(lst), n):
            yield lst[i:i + n]

    if len(table_strings)%2:
        table_strings.append("&")

    pair_strings = map( lambda l: " & ".join( l )+"\\\\", chunks( table_strings, 2 ) )

    table_string = \
    """
    \\begin{{table}}
       \caption{{Typical values (90\% quantiles) and maxmimal values of the {post}-fit variations of all signal regions.}} \label{{tab:{post}-fit}}
       \center
          \\begin{{tabular}}{{r|c|c|r|c|c}}
                systematic  & typical & max & systematic & typical & max \\\\\\hline
    {systematics}
          \end{{tabular}}
    \end{{table}}
    """.format(systematics = "\n".join(pair_strings), post="post" if postFit else "pre")

    print table_string

#cardFile  = "/afs/hephy.at/data/cms05/StopsDileptonLegacy/results/v7/2016/fitAll/cardFiles/T2tt/observed/T2tt_800_100.txt"
#cardFile = "/afs/hephy.at/data/llechner01/TTGammaEFT/cache/analysis/2016/limits/cardFiles/defaultSetup/observed/SR4pM3_VG4p_misDY4p_misTT2_incl.txt"


#signal         TTZ            TTJets         TTXNoZ         DY             multiBoson

# Bin20: SF MET_significance>=12&&MET_significance<50&&dl_mt2blbl>=0&&dl_mt2blbl<100&&dl_mt2ll>=100&&dl_mt2ll<140
# Bin21: EMu MET_significance>=12&&MET_significance<50&&dl_mt2blbl>=0&&dl_mt2blbl<100&&dl_mt2ll>=100&&dl_mt2ll<140
# Bin22: SF MET_significance>=50&&dl_mt2blbl>=0&&dl_mt2blbl<100&&dl_mt2ll>=100&&dl_mt2ll<140
# Bin23: EMu MET_significance>=50&&dl_mt2blbl>=0&&dl_mt2blbl<100&&dl_mt2ll>=100&&dl_mt2ll<140
# Bin24: SF MET_significance>=12&&MET_significance<50&&dl_mt2blbl>=100&&dl_mt2blbl<200&&dl_mt2ll>=100&&dl_mt2ll<140
# Bin25: EMu MET_significance>=12&&MET_significance<50&&dl_mt2blbl>=100&&dl_mt2blbl<200&&dl_mt2ll>=100&&dl_mt2ll<140
# Bin26: SF MET_significance>=50&&dl_mt2blbl>=100&&dl_mt2blbl<200&&dl_mt2ll>=100&&dl_mt2ll<140
# Bin27: EMu MET_significance>=50&&dl_mt2blbl>=100&&dl_mt2blbl<200&&dl_mt2ll>=100&&dl_mt2ll<140
# Bin28: SF MET_significance>=12&&MET_significance<50&&dl_mt2blbl>=200&&dl_mt2ll>=100&&dl_mt2ll<140
# Bin29: EMu MET_significance>=12&&MET_significance<50&&dl_mt2blbl>=200&&dl_mt2ll>=100&&dl_mt2ll<140
# Bin30: SF MET_significance>=50&&dl_mt2blbl>=200&&dl_mt2ll>=100&&dl_mt2ll<140
# Bin31: EMu MET_significance>=50&&dl_mt2blbl>=200&&dl_mt2ll>=100&&dl_mt2ll<140
# Bin32: SF MET_significance>=12&&MET_significance<50&&dl_mt2blbl>=0&&dl_mt2blbl<100&&dl_mt2ll>=140&&dl_mt2ll<240
# Bin33: EMu MET_significance>=12&&MET_significance<50&&dl_mt2blbl>=0&&dl_mt2blbl<100&&dl_mt2ll>=140&&dl_mt2ll<240
# Bin34: SF MET_significance>=50&&dl_mt2blbl>=0&&dl_mt2blbl<100&&dl_mt2ll>=140&&dl_mt2ll<240
# Bin35: EMu MET_significance>=50&&dl_mt2blbl>=0&&dl_mt2blbl<100&&dl_mt2ll>=140&&dl_mt2ll<240
# Bin36: SF MET_significance>=12&&MET_significance<50&&dl_mt2blbl>=100&&dl_mt2blbl<200&&dl_mt2ll>=140&&dl_mt2ll<240
# Bin37: EMu MET_significance>=12&&MET_significance<50&&dl_mt2blbl>=100&&dl_mt2blbl<200&&dl_mt2ll>=140&&dl_mt2ll<240
# Bin38: SF MET_significance>=50&&dl_mt2blbl>=100&&dl_mt2blbl<200&&dl_mt2ll>=140&&dl_mt2ll<240
# Bin39: EMu MET_significance>=50&&dl_mt2blbl>=100&&dl_mt2blbl<200&&dl_mt2ll>=140&&dl_mt2ll<240
# Bin40: SF MET_significance>=12&&MET_significance<50&&dl_mt2blbl>=200&&dl_mt2ll>=140&&dl_mt2ll<240
# Bin41: EMu MET_significance>=12&&MET_significance<50&&dl_mt2blbl>=200&&dl_mt2ll>=140&&dl_mt2ll<240
# Bin42: SF MET_significance>=50&&dl_mt2blbl>=200&&dl_mt2ll>=140&&dl_mt2ll<240
# Bin43: EMu MET_significance>=50&&dl_mt2blbl>=200&&dl_mt2ll>=140&&dl_mt2ll<240
# Bin44: SF MET_significance>=12&&dl_mt2blbl>=0&&dl_mt2ll>=240
# Bin45: EMu MET_significance>=12&&dl_mt2blbl>=0&&dl_mt2ll>=240
