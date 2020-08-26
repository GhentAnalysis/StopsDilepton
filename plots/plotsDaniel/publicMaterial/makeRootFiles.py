import ROOT

from StopsDilepton.tools.helpers import getObjFromFile, writeObjToFile

objects = {}

objects['T2tt'] = {
    'observed': '/afs/hephy.at/user/d/dspitzbart//www/stopsDileptonLegacy/limits/T2tt/v9/comb/CWR_smooth_it1_k5a/limitXSEC.root',
    'expected': '/afs/hephy.at/user/d/dspitzbart//www/stopsDileptonLegacy/limits/T2tt/v9/comb/FR_smooth_it1_k5a_expected/limitXSEC.root',
    'objects': ['contour_exp', 'contour_exp_up', 'contour_exp_down', 'contour_obs', 'contour_obs_up', 'contour_obs_down', 'contour_corr_obs', 'contour_corr_obs_up', 'contour_corr_obs_down']
}

objects['T2bW'] = {
    'observed': '/afs/hephy.at/user/d/dspitzbart//www/stopsDileptonLegacy/limits/T2bW/v9/comb/FR/limitXSEC.root',
    'expected': '/afs/hephy.at/user/d/dspitzbart//www/stopsDileptonLegacy/limits/T2bW/v9/comb/FR_expected/limitXSEC.root',
    'objects': ['contour_exp', 'contour_exp_up', 'contour_exp_down', 'contour_obs', 'contour_obs_up', 'contour_obs_down']
}

objects['T8bbllnunu_XCha0p5_XSlep0p95'] = {
    'observed': '/afs/hephy.at/user/d/dspitzbart//www/stopsDileptonLegacy/limits/T8bbllnunu_XCha0p5_XSlep0p95/v9/comb/FR/limitXSEC.root',
    'expected': '/afs/hephy.at/user/d/dspitzbart//www/stopsDileptonLegacy/limits/T8bbllnunu_XCha0p5_XSlep0p95/v9/comb/FR_expected/limitXSEC.root',
    'objects': ['contour_exp', 'contour_exp_up', 'contour_exp_down', 'contour_obs', 'contour_obs_up', 'contour_obs_down']
}

objects['T8bbllnunu_XCha0p5_XSlep0p5'] = {
    'observed': '/afs/hephy.at/user/d/dspitzbart//www/stopsDileptonLegacy/limits/T8bbllnunu_XCha0p5_XSlep0p5/v9/comb/FR_smooth_it1_k5a/limitXSEC.root',
    'expected': '/afs/hephy.at/user/d/dspitzbart//www/stopsDileptonLegacy/limits/T8bbllnunu_XCha0p5_XSlep0p5/v9/comb/FR_smooth_it1_k5a_expected/limitXSEC.root',
    'objects': ['contour_exp', 'contour_exp_up', 'contour_exp_down', 'contour_obs', 'contour_obs_up', 'contour_obs_down']
}

objects['T8bbllnunu_XCha0p5_XSlep0p05'] = {
    'observed': '/afs/hephy.at/user/d/dspitzbart//www/stopsDileptonLegacy/limits/T8bbllnunu_XCha0p5_XSlep0p05/v9/comb/FR_smooth_it1_k5a/limitXSEC.root',
    'expected': '/afs/hephy.at/user/d/dspitzbart//www/stopsDileptonLegacy/limits/T8bbllnunu_XCha0p5_XSlep0p05/v9/comb/FR_smooth_it1_k5a_expected/limitXSEC.root',
    'objects': ['contour_exp', 'contour_exp_up', 'contour_exp_down', 'contour_obs', 'contour_obs_up', 'contour_obs_down']
}

for signal in objects.keys():
    print signal

    objects['objects'] = {}

    # first observed
    c_obs = getObjFromFile(objects[signal]['observed'], 'cCONT_asdf')
    tmp = c_obs.GetPrimitive('temperature')
    tmp.SetName("temp_obs")
    writeObjToFile('%s.root'%signal, tmp)
    
    for obj in objects[signal]['objects']:
        tmp = c_obs.GetPrimitive(obj)
        writeObjToFile('%s.root'%signal, tmp, update=True)

    del c_obs

    c_exp = getObjFromFile(objects[signal]['expected'], 'cCONT_asdf')
    tmp = c_exp.GetPrimitive('temperature')
    tmp.SetName("temp_exp")
    writeObjToFile('%s.root'%signal, tmp, update=True)
    del c_exp

#writeObjToFile( , , update=True)

