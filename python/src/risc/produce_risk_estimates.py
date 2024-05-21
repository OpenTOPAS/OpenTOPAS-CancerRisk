import glob
import os
import pandas as pd

from risc.cancer_risk import TsRiskModelWrapper

def unisex_SEER_to_specific(organ, sex):
    if organ in ['testesuterus', 'prostateovary']:
        lookup = {'testesuterus': {'m': 'testes',
                                   'f': 'uterus'},
                  'prostateovary': {'m': 'prostate',
                                    'f': 'ovary'}
                 }
        return lookup[organ][sex]
    else:
        return organ

def run_risk_model(resultsdir, patient_info, risk_parameters, SEER_dir, dvhFile):
    TsRiskModelWrapper(patient_info['patient'],
                       risk_parameters['cancer_model'],
                       patient_info['sex'],
                       risk_parameters['parameter_set'],
                       risk_parameters['organ'],
                       patient_info['age'],
                       SEER_dir,
                       risk_parameters['attained_age'],
                       risk_parameters['fractions'],
                       dvhFile)

def process_directory_with_parameters(resultsdir, dvh_dir, patient_info, 
                                      risk_parameters, exclude_organs, 
                                      SEER_dir):
    dvh_dir = os.path.join(dvh_dir, '')
    os.system('mkdir -p {:s}'.format(resultsdir))
    os.chdir(resultsdir)

    files = glob.glob(dvh_dir + '*VolHist.csv')
    for file in files:
        organ = os.path.basename(file).split('_')[0]
        organ = unisex_SEER_to_specific(organ, patient_info['sex'])
        avail_organs = get_available_organs(risk_parameters['parameter_set'],
                                            risk_parameters['cancer_model'],
                                            SEER_dir)
        if organ in exclude_organs:
            print('Skipping (excluded): ', organ)
            continue
        if organ not in avail_organs:
            print('Skipping (not in parameter set): ', organ)
            continue
        risk_parameters['organ'] = organ
        run_risk_model(resultsdir, patient_info, risk_parameters, SEER_dir, file)
    os.chdir('..')

def run_different_parameter_sets(resultsdir, patient_info, risk_parameters, 
                                 exclude_organs, dvh_dir, SEER_dir):
    os.system('mkdir -p {:s}'.format(resultsdir))
    cwd = os.getcwd()
    os.chdir(resultsdir)

    param_sets = ['BEIRVII', 'Schneider']
    cancer_types = ['sarcoma', 'carcinoma']

    for param_set in param_sets:
        for cancer_type in cancer_types:
            risk_parameters['parameter_set'] = param_set
            risk_parameters['cancer_model'] = cancer_type
            resultsdir = '{:s}_{:s}/'.format(param_set, cancer_type)
            process_directory_with_parameters(resultsdir, dvh_dir,
                                              patient_info, risk_parameters,
                                              exclude_organs, SEER_dir)
    os.chdir(cwd)

def get_available_organs(parameter_set, cancer_type, data_dir):
    parameter_set_files = {'BEIRVII': {'sarcoma': 'BEIRVII.txt',
                                       'carcinoma': 'BEIRVII_repopulation.txt'},
                           'Schneider': {'sarcoma': 'Schneider.txt',
                                         'carcinoma': 'Schneider.txt'}
                          }
    file = parameter_set_files[parameter_set][cancer_type]
    df = pd.read_csv(data_dir + file, comment='#', delimiter=r'\s+')
    SEER_names = list(set(df['SEER_name'].values))
    print(SEER_names)
    return SEER_names
