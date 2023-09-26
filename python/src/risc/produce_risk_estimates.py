import glob
import os

from risc.cancer_risk import TsRiskModelWrapper

def run_risk_model(resultsdir, patient_info, risk_parameters, SEER_dir, dvhFile):
    print(patient_info, risk_parameters)
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
    os.system('mkdir -p {:s}'.format(resultsdir))
    os.chdir(resultsdir)

    files = glob.glob(dvh_dir + '*VolHist.csv')
    for file in files:
        organ = os.path.basename(file).split('_')[0]
        if organ in exclude_organs:
            print('Skipping: ', organ)
            continue
        print(organ)
        risk_parameters['organ'] = organ
        run_risk_model(resultsdir, patient_info, risk_parameters, SEER_dir, file)
    os.chdir('..')

def run_different_parameter_sets(resultsdir, patient_info, risk_parameters, 
                                 exclude_organs, dvh_dir, SEER_dir):
    os.system('mkdir -p {:s}'.format(resultsdir))
    os.chdir(resultsdir)

    param_sets = ['BEIR', 'Schneider']
    cancer_types = ['sarcoma', 'carcinoma']

    for param_set in param_sets:
        for cancer_type in cancer_types:
            risk_parameters['parameter_set'] = param_set
            risk_parameters['cancer_model'] = cancer_type
            resultsdir = '{:s}_{:s}/'.format(param_set, cancer_type)
            process_directory_with_parameters(resultsdir, dvh_dir,
                                              patient_info, risk_parameters,
                                              exclude_organs, SEER_dir)

