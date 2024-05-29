# CancerRisk extension
Repository for code related to Cancer Risk project

Authors:
 - Isaac Meyer
 - Alejandro Bertolet
 - Harald Paganetti

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

# Installation and Examples

## Installing TOPAS with extension

1. Clone OpenTOPAS repo v4.0.0
   git clone -b v4.0.0 https://github.com/OpenTOPAS/OpenTOPAS
2. Clone the cancerrisk extension
   git clone --recurse-submodules https://github.com/OpenTOPAS/TOPAS-CancerRisk
3. Install Geant4v11.1.2 following https://opentopas.github.io/installation.html
4. Build TOPAS with cmake, possible build script below (change directories as appropriate)

```
#!/bin/bash

FLAVOR="CancerRisk"
BUILDDIR="build_$FLAVOR"
rm -rf $BUILDDIR
mkdir $BUILDDIR
cd $BUILDDIR

export Geant4_DIR=/Applications/geant4-v11.1.2-install
cmake -DTOPAS_TYPE=public \
      -DTOPAS_EXTENSIONS_DIR="../../TOPAS-CancerRisk" \
      -DCMAKE_INSTALL_PREFIX="../topas-install-$FLAVOR" \
      ../
make -j12 install 
```

## Running example `example/AF_SOBP_Liver`
1. Check that the environment variable `$TOPAS_G4_DATA` is appropriately set
2. Set appropriate path to ICRP 145 phantom data for parameter `Ge/MeshPhantom/PhantomDirector`

## Installing python wrapper
The python wrapper requires generating a TOPAS shared library. 
1. Add the following to `OpenTOPAS/extensions/CMakeLists.txt.in`:
    ```
    add_library( extensions_shared SHARED
    TsExtensionManager.cc
    TsExtensionManager.hh
    )
    target_link_libraries(extensions_shared
        main
        parameter
        chemistry
        geometry
        extensions
        graphics
        material
        physics
        variance
        filtering
        scoring
        outcome
        io
        sequence
        primary
        gdcmMSFF 
        ${Geant4_LIBRARIES}
    )
    ```
2. Set the environment variable `$TOPAS_CANCER_RISK_BUILD_DIR`
3. The python module can be installed by running `pip install -e .` in the `TOPAS-CancerRisk/python` folder
4. The DVH processing example in `TOPAS-CancerRisk/python/examples/` can be run after running TOPAS on the `example/AF_SOBP_Liver` example

# Issues
Please use the issues tab on github to report problems or request features

