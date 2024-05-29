import ctypes
from ctypes import (c_int, c_double, POINTER, c_wchar_p, c_char_p)
import numpy as np
import os
import subprocess
from numpy.ctypeslib import ndpointer

if 'TOPAS_CANCER_RISK_BUILD_DIR' not in os.environ:
    print('Please set the TOPAS_CANCER_RISK_BUILD_DIR environment variable')
    exit(1)
else:
    libpath = os.environ['TOPAS_CANCER_RISK_BUILD_DIR'] + '/extensions/libextensions_shared.dylib'

_topaslib = ctypes.cdll.LoadLibrary(libpath)

_topaslib.TsRiskModelWrapper.argtypes = (c_char_p, c_char_p, c_char_p, 
                                         c_char_p, c_char_p, c_int,
                                         c_char_p, c_int, c_int,
                                         c_char_p)
_topaslib.TsRiskModelWrapper.restype = None

_topaslib.testHistogramSumming.argtypes = (c_char_p, c_char_p, c_char_p)
_topaslib.testHistogramSumming.restype = None

def TsRiskModelWrapper(patient, cancerModel, sex, parameterSet, 
                       organAtRisk, ageAtExposure, SEERDirectory, 
                       attainedAge, numberOfFractions, dvhFile):
    patient = patient.encode('utf-8')
    cancerModel = cancerModel.encode('utf-8')
    sex = sex.encode('utf-8')
    parameterSet = parameterSet.encode('utf-8')
    organAtRisk = organAtRisk.encode('utf-8').lower()
    SEERDirectory = SEERDirectory.encode('utf-8')
    dvhFile = dvhFile.encode('utf-8')

    _topaslib.TsRiskModelWrapper(c_char_p(patient), 
                                 c_char_p(cancerModel), 
                                 c_char_p(sex), 
                                 c_char_p(parameterSet), 
                                 c_char_p(organAtRisk), 
                                 c_int(ageAtExposure), 
                                 c_char_p(SEERDirectory),
                                 c_int(attainedAge), 
                                 c_int(numberOfFractions),
                                 c_char_p(dvhFile)
                                )

def testHistogramSumming(file1, file2, outfile):
    file1 = file1.encode('utf-8')
    file2 = file2.encode('utf-8')
    outfile = outfile.encode('utf-8')

    _topaslib.testHistogramSumming(c_char_p(file1),
                                   c_char_p(file2),
                                   c_char_p(outfile)
                                  )

def run_TsRiskModelWrapper():
    # TODO: move to examples
    patient = 'TestPatient'
    cancerModel = 'sarcoma'
    sex = 'm'
    parameterSet = 'BEIRVII'
    # parameterSet = 'Schneider'
    organAtRisk = 'liver'
    ageAtExposure = 30
    attainedAge = 30
    numberOfFractions = 20
    dvhFile = '/Users/isaacmeyer/research/secondary_risk/topas_use_dvh_risk/DoseAtPhantom_VolHist.csv'

    SEERDirectory = '/Users/isaacmeyer/topas_extensions/TOPAS-CancerRisk/data/'
    TsRiskModelWrapper(patient, cancerModel, sex, parameterSet, 
                       organAtRisk, ageAtExposure, SEERDirectory,
                       attainedAge, numberOfFractions, dvhFile)

if __name__=='__main__':
    run_TsRiskModelWrapper()
