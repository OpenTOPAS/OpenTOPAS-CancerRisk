import risc

patient = 'TestPatient'
cancerModel = 'sarcoma'
sex = 'F'
parameterSet = 'Schneider'
organAtRisk = 'Liver'
ageAtExposure = 35
# Data path must be absolute for now
SEERDirectory = '/Users/isaacmeyer/topas_extensions/TOPAS-CancerRisk/data/'
attainedAge = 70  
numberOfFractions = 28
dvhFile = '../../examples/AF_SOBP_LIVER/DoseLiver_VolHist.csv'

risc.cancer_risk.TsRiskModelWrapper(patient, cancerModel, sex, parameterSet, 
                                    organAtRisk, ageAtExposure, SEERDirectory, 
                                    attainedAge, numberOfFractions, dvhFile)

