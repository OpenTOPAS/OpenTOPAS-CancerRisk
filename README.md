# CancerRisk extension
Repository for code related to Cancer Risk project

Authors:
 - [Isaac Meyer](imeyer@mgh.harvard.edu)
 - [Alejandro Bertolet](abertoletreina@mgh.harvard.edu)
 - [Harald Paganetti](hpaganetti@mgh.harvard.edu)

# TODO:
    - add examples
    - clean up this README

# SEER data
This extension uses statistics from SEER (Surveillance, Epidemiology and End Results) Program to estimate risks in relation to each organ. This information is included in *cvs* files that are provided with this extension. All these files need to be included in the same directory, and its full path needs to be specified in the parameter file as follows:
- s:Sc/CancerRisk/DataDirectory = "/Users/.../SEERData/"

# How to use it
To include the cancer risk prediction as an outcome of the calculation, please include "CancerRisk" within the options in the parameter sv:Sc/OAR/OutcomeModelName (OAR represents the name of each organ at risk for which outcome is desired). This process needs to be done for each OAR to be considered. Dose-volume histogram needs to be calculated for each OAR of interest. For each OAR, it is possible to specify the parameter set used for outcome (BEIR or Schneider). It is assumed by default the Schneider set. To change to BEIR, it is *optional* to use the following parameter for each OAR:
- s:Sc/OAR/CancerRisk/ParameterSet = "BEIR" (otherwise Schneider set will be used).

The name of the OAR *needs* to be specified:
s:Sc/OAR/OrganName.
This is used to seek for the specific SEER parameters for that organs. The available data include: "brain", "breast", "colonrectum", "esophagus", "kidneys", "laryinx", "liver", "lung", "lungbronchus", "pancreas", "pharynx", "prostateovary", "stomach", "testesuterus", "thyroid", "urinarybladder". Plase make sure that you use one of these exact strings.

The rest of parameters, however, need to be specified only once.
It is *optional* to specify the following parameters:
- s:Sc/CancerRisk/CancerModel = "carcinoma" or "sarcoma". If not specified, it will use models for carcinoma.
- s:Sc/CancerRisk/PatientName. If not specified, it will use "Patient1". This is used to produce output files.

It is *mandatory* to specify the following parameters:
- s:Sc/CancerRisk/Sex = "male" or "female".
- i:Sc/CancerRisk/AttainedAge = attained age of the patient.
- i:Sc/CancerRisk/AgeAtExposure = age at which patient received radiation.
- i:Sc/CancerRisk/NumberOfFractions = number of fractions for the treatment.

# MRCP Phantoms and scoring
If instead of calculating on a patient, a MRCP tetrahedron-based phantom is preferred, this can be loaded using the classes in the module "MRCP".
Phantoms need to be defined as geometry components in the parameter file (*Type="TsMRCP"*), and the directory where the files .node, .ele and .material provided by IRCP are located needs to be specified (see example "15F_SOBP_Liver.txt").
A specific class of scorer is needed for MRCP phantoms by specifying (*Quantity="TsMRCPScorer"*). By default, dose to water is computed, but other material can be used to calculate the dose by using the parameter *Material*.
For each organ involved it is necessary to add a new scorer, specifying the list parameter *Organ=1 "Name"*, where *Name* should be one (or more) of the materials listed in the .material file of the MRCP phantom, e.g., "Liver" or "Stomach_contents". This will restrict the dose considered only to the organ(s) of interest.
