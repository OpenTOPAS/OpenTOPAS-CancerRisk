// Outcome Model for CancerRisk
#include "TsRiskModel.hh"
#include "TsParameterManager.hh"
#include "MeshGeomTools.hh"

#include "G4Integrator.hh"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

using namespace std;

TsRiskModel::TsRiskModel(TsParameterManager* pm, G4String parmName) : TsVOutcomeModel(pm, parmName), fPm(pm), fModelName("CancerRisk")
{
	// Set carcinoma as default
	fCancerModel = "carcinoma";
	fPatient = "Patient1";
	// Parameters used to fit the models
	fLatency = 5; 	// Years
	fDDREF = 1.0;	// Dose and dose-rate effectiveness ratio
	ResolveParameters();
	ReadLifetimeRiskTable();
	ReadSEEROrganSpecificTable(fOrganAtRisk);
}

TsRiskModel::TsRiskModel(G4String patient, G4String cancerModel, G4String sex, G4String parameterSet,
                         G4String organAtRisk, G4int ageAtExposure, G4String SEERDirectory,
                         G4int attainedAge, G4int numberOfFractions, G4String dvhFile)
           : TsVOutcomeModel(NULL, "None"),
           fPatient(patient), fCancerModel(cancerModel), fSex(sex), fParameterSet(parameterSet),
           fOrganAtRisk(organAtRisk), fAgeAtExposure(ageAtExposure), fSEERDirectory(SEERDirectory),
           fAttainedAge(attainedAge), fNumberOfFractions(numberOfFractions)
{
    G4cout << "Constructing risk model " <<  G4endl;
	fLatency = 5; 	// Years
	fDDREF = 1.0;	// Dose and dose-rate effectiveness ratio
    G4cout << "Patient: " << fPatient << G4endl;
    G4cout << "Age: " << fAgeAtExposure << G4endl;
    G4cout << "OrganAtRisk: " << fOrganAtRisk << G4endl;

    // Warnings for incorrect input
	if (fParameterSet != "BEIR" && fParameterSet != "Schneider")
	{
		G4cout << "Parameter set not specified or not recognized. Schneider parameter set will be used." << G4endl;
	}

	if (fSex == "m" || fSex == "Male" || fSex == "M" || fSex == "man" || fSex == "Man") fSex = "male";
	if (fSex == "f" || fSex == "Female" || fSex == "F" || fSex == "woman" || fSex == "Woman" || fSex == "w" || fSex == "W" ) fSex = "female";
	if (fSex != "male" && fSex != "female")
	{
		G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << "Sc/CancerRisk/Sex needs to be set to one of the following options: " << G4endl;
		G4cerr << "male, female, m, f, Male, Female, M, F, man, woman, Man, Woman, w or W." << G4endl;
		exit(1);
	}
	ReadLifetimeRiskTable();
	ReadSEEROrganSpecificTable(fOrganAtRisk);
	ReadDVHcsv(dvhFile);
    // TODO: check expected format of dose and volume inputs
    G4double value = Initialize(DVHDose, DVHVolume);
    G4cout << "    Probability: " << value << " %" << G4endl;
}

TsRiskModel::~TsRiskModel() { }

void TsRiskModel::ResolveParameters()
{
	fSEERDirectory = fPm->GetStringParameter("Sc/CancerRisk/DataDirectory");
	if (*fSEERDirectory.rbegin() != '/') fSEERDirectory += "/";
	fParameterSet = fPm->GetStringParameter(GetFullParmName(fModelName, "ParameterSet"));
	if (fParameterSet != "BEIR" && fParameterSet != "Schneider")
	{
		G4cout << "Parameter set not specified or not recognized. Schneider parameter set will be used." << G4endl;
	}
	fOrganAtRisk = fPm->GetStringParameter(GetFullParmName(fModelName, "OrganName"));

	fCancerModel = fPm->GetStringParameter("Sc/CancerRisk/CancerModel");
	if (fCancerModel != "sarcoma" && fCancerModel != "carcinoma") fCancerModel = "carcinoma";
	fPatient = fPm->GetStringParameter("Sc/CancerRisk/PatientName");
	fSex = fPm->GetStringParameter("Sc/CancerRisk/Sex");
	if (fSex == "m" || fSex == "Male" || fSex == "M" || fSex == "man" || fSex == "Man") fSex = "male";
	if (fSex == "f" || fSex == "Female" || fSex == "F" || fSex == "woman" || fSex == "Woman" || fSex == "w" || fSex == "W" ) fSex = "female";
	if (fSex != "male" && fSex != "female")
	{
		G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << "Sc/CancerRisk/Sex needs to be set to one of the following options: " << G4endl;
		G4cerr << "male, female, m, f, Male, Female, M, F, man, woman, Man, Woman, w or W." << G4endl;
		exit(1);
	}
	fAttainedAge = fPm->GetIntegerParameter("Sc/CancerRisk/AttainedAge");
	fAgeAtExposure = fPm->GetIntegerParameter("Sc/CancerRisk/AgeAtExposure");
	fNumberOfFractions = fPm->GetIntegerParameter("Sc/CancerRisk/NumberOfFractions");
}

G4double TsRiskModel::Initialize(std::vector<G4double> dose, std::vector<G4double> volume)
{
	SetOrganSpecificParameters();

	G4String fileName = fOrganAtRisk + "_diffDVH.csv";
	fstream fout;
	// Creates a new file, overwritting if exists
	fout.open(fileName, ios::out | ios::trunc);
	// Insert headers
	for (G4int i = 0; i < dose.size(); i++)
	{
		fout << dose[i] << "," << volume[i] << "\n";
	}
	fout.close();

	JobResults.push_back(fPatient);

	G4double OED;
	if (fCancerModel == "sarcoma") OED = GetOEDSarcomaModel(dose, volume);
	if (fCancerModel == "carcinoma") OED = GetOEDCarcinomaModel(dose, volume);
	vector<G4double> LARERRCum = LARBasedOnERR(OED);
	vector<G4double> LAREARCum = LARBasedOnEAR(OED);
	G4double avgLAR = AverageLAReAndLARa(LARERRCum, LAREARCum);

	// Write summary output
	JobHeaders.push_back("DVHFile"); JobHeaders.push_back("MeanDose");
	JobHeaders.push_back("ERR5"); JobHeaders.push_back("ERR10"); JobHeaders.push_back("ERR15"); JobHeaders.push_back("ERR50"); JobHeaders.push_back("ERR70");
	JobHeaders.push_back("EAR5"); JobHeaders.push_back("EAR10"); JobHeaders.push_back("EAR15"); JobHeaders.push_back("EAR50"); JobHeaders.push_back("EAR70");
	JobHeaders.push_back("LAR5"); JobHeaders.push_back("LAR10"); JobHeaders.push_back("LAR15"); JobHeaders.push_back("LAR50"); JobHeaders.push_back("LAR70");

	G4String outputFileName = fPatient + "_Table_" + fOrganAtRisk + ".csv";
	WriteOutputFile(outputFileName, JobResults, JobHeaders);
	G4cout << "Results for ERR, EAR and LAR 5, 10, 15, 50 and 70 were written into the file " + outputFileName << G4endl;
	G4cout << "NOTE: Probability below shows the Lifetime Atributtable Risk (LAR) at the attained age averaged using ERR and EAR formalisms, expressed in risk per 10,000 person-years." << G4endl;
	return avgLAR / 100;
}

std::vector<G4double> TsRiskModel::LARBasedOnERR(G4double OED)
{
	vector<G4double> vLARerr, vLARerrCum, vERR, vAge, Sa, SaOrgan;
	G4double JobERR5, JobERR10, JobERR15, JobERR50, JobERR70;
	G4double S;
	G4double LARExcess = 0.0;
	if (fSex == "female")
	{
		Sa = SaFemale;
		SaOrgan = SaOrganFemale;
		S = 0.17;
	}
	if (fSex == "male")
	{
		Sa = SaMale;
		SaOrgan = SaOrganMale;
		S = -0.17;
	}
	G4double eStar = 0;
	if (fAgeAtExposure < 30) eStar = fAgeAtExposure - 30;
	if (fOrganAtRisk == "thyroid" || fOrganAtRisk == "breast") eStar = fAgeAtExposure;
	for (G4int age = fAgeAtExposure + fLatency; age <= fAttainedAge + fLatency; age++)
	{
		G4double ERR = OED * fBetaERR * exp(fGamma_e_ERR * eStar + fGamma_a_ERR * log((G4double)age/(G4double)fParAttainedAge)) * (1 + S);
		G4double LARe = ERR * Sa[age] / Sa[fAgeAtExposure] * SaOrgan[age];
		vAge.push_back(age);
		vLARerr.push_back(LARe / 1000);
		vERR.push_back(ERR);
		if (age == fAgeAtExposure + 5) JobERR5 = ERR;
		if (age == fAgeAtExposure + 10) JobERR10 = ERR;
		if (age == fAgeAtExposure + 15) JobERR15 = ERR;
		if (age == 50) JobERR50 = ERR;
		if (age == 70) JobERR70 = ERR;
		LARExcess += LARe;
		vLARerrCum.push_back(LARExcess / 1000);
	}
	LARExcess /= 1000;
	JobResults.push_back(to_string(JobERR5)); JobResults.push_back(to_string(JobERR10)); JobResults.push_back(to_string(JobERR15));
	JobResults.push_back(to_string(JobERR50)); JobResults.push_back(to_string(JobERR70));
	// Write outputs
	G4String nameOutput = fPatient + "_ERR_and_LAR_ERR_" + fOrganAtRisk + ".csv";
	vector<G4String> headers;
	vector<vector<G4double>> quantities;
	quantities.push_back(vAge); headers.push_back("Age");
	quantities.push_back(vERR); headers.push_back("ERR");
	quantities.push_back(vLARerrCum); headers.push_back("LAR_err_Cum");
	WriteOutputFile(nameOutput, quantities, headers);
	G4cout << "Results for LAR based on ERR were written into the file " + nameOutput << G4endl;

	return vLARerrCum;
}

std::vector<G4double> TsRiskModel::LARBasedOnEAR(G4double OED)
{
	vector<G4double> vEAR, vLARear, vLARearCum, vAge, Sa, SaOrgan;
	G4double JobEAR5, JobEAR10, JobEAR15, JobEAR50, JobEAR70;
	G4double LARabs = 0.0;
	if (fSex == "female")
	{
		Sa = SaFemale;
		SaOrgan = SaOrganFemale;
	}
	if (fSex == "male")
	{
		Sa = SaMale;
		SaOrgan = SaOrganMale;
	}
	G4double eStar = 0;
	if (fAgeAtExposure < 30) eStar = fAgeAtExposure - 30;
	if (fOrganAtRisk == "thyroid" || fOrganAtRisk == "breast") eStar = fAgeAtExposure;
	for (G4int age = fAgeAtExposure + fLatency; age <= fAttainedAge + fLatency; age++)
	{
		G4double EAR = OED * fBetaEAR * exp(fGamma_e_EAR * eStar + fGamma_a_EAR * log((G4double)age/(G4double)fParAttainedAge));
		G4double LARe = EAR * Sa[age] / Sa[fAgeAtExposure];
		LARabs += LARe;
		vAge.push_back(age);
		vLARear.push_back(LARe / 100);
		vEAR.push_back(EAR / 100);
		if (age == fAgeAtExposure + 5) JobEAR5 = EAR / 100;
		if (age == fAgeAtExposure + 10) JobEAR10 = EAR / 100;
		if (age == fAgeAtExposure + 15) JobEAR15 = EAR / 100;
		if (age == 50) JobEAR50 = EAR / 100;
		if (age == 70) JobEAR70 = EAR / 100;
		vLARearCum.push_back(LARabs / 100);
	}
	JobResults.push_back(to_string(JobEAR5)); JobResults.push_back(to_string(JobEAR10)); JobResults.push_back(to_string(JobEAR15));
	JobResults.push_back(to_string(JobEAR50)); JobResults.push_back(to_string(JobEAR70));
	LARabs /= 100;
	// Write outputs
	G4String nameOutput = fPatient + "_EAR_and_LAR_EAR_" + fOrganAtRisk + ".csv";
	vector<G4String> headers;
	vector<vector<G4double>> quantities;
	quantities.push_back(vAge); headers.push_back("Age");
	quantities.push_back(vEAR); headers.push_back("EAR");
	quantities.push_back(vLARearCum); headers.push_back("LAR_ear_Cum");
	WriteOutputFile(nameOutput, quantities, headers);
	G4cout << "Results for LAR based on EAR were written into the file " + nameOutput << G4endl;

	return vLARearCum;
}

G4double TsRiskModel::AverageLAReAndLARa(std::vector<G4double> LARERRCum, std::vector<G4double> LAREARCum)
{
	G4double wERR = 0.0;
	G4double wEAR = 1.0;
	G4double wAvg = 0.0;
	vector<G4double> vLARavg, vAge;
	G4double JobLAR5, JobLAR10, JobLAR15, JobLAR50, JobLAR70, LARAttainedAge;
	for (G4int age = fAgeAtExposure + fLatency; age <= fAttainedAge + fLatency; age++)	{ vAge.push_back(age); }
	if (fParameterSet == "BEIR")
	{
		if (fOrganAtRisk == "breast" || fOrganAtRisk == "bone_marrow") wERR = 0.0;
		else if (fOrganAtRisk == "thyroid" || fOrganAtRisk == "skin") wERR = 1.0;
		else if (fOrganAtRisk == "lung") wERR = 0.3;
		else wERR = 0.5;
		wEAR = 1 - wERR;
	}
	for (G4int i = 0; i <= fAttainedAge - fAgeAtExposure; i++)
	{
		wAvg = pow(LARERRCum[i], wERR) * pow(LAREARCum[i], wEAR) / fDDREF;
		vLARavg.push_back(wAvg);
		if (i == fAgeAtExposure + 5 - (fAgeAtExposure + fLatency)) JobLAR5 = wAvg;
		if (i == fAgeAtExposure + 10 - (fAgeAtExposure + fLatency)) JobLAR10 = wAvg;
		if (i == fAgeAtExposure + 15 - (fAgeAtExposure + fLatency)) JobLAR15 = wAvg;
		if (i == 50 - (fAgeAtExposure + fLatency)) JobLAR50 = wAvg;
		if (i == 70 - (fAgeAtExposure + fLatency)) JobLAR70 = wAvg;
		if (i == fAttainedAge - fAgeAtExposure - 1) LARAttainedAge = wAvg;
	}
	JobResults.push_back(to_string(JobLAR5)); JobResults.push_back(to_string(JobLAR10)); JobResults.push_back(to_string(JobLAR15));
	JobResults.push_back(to_string(JobLAR50)); JobResults.push_back(to_string(JobLAR70));
	// Write output
	G4String nameOutput = fPatient + "_LAR_" + fOrganAtRisk + ".csv";
	vector<G4String> headers;
	vector<vector<G4double>> quantities;
	quantities.push_back(vAge); headers.push_back("Age");
	quantities.push_back(vLARavg); headers.push_back("LAR averaged");
	WriteOutputFile(nameOutput, quantities, headers);
	G4cout << "Results for Averaged LAR were written into the file " + nameOutput << G4endl;
	return LARAttainedAge;
}

G4double TsRiskModel::GetOEDSarcomaModel(std::vector<G4double> dose, std::vector<G4double> volume)
{
	G4double beta = fAlpha / fAlphaBeta;
	G4double sumVolume = 0;
	G4double OEDs = 0;
	G4double meanDose = 0;
	for (G4int i = 0; i < dose.size(); i++)
	{
		G4double dvhD = dose[i];
		G4double dvhV = volume[i];
		G4double alphaPrime = fAlpha + beta * dvhD / fNumberOfFractions;
		sumVolume += dvhV;
		meanDose += dvhD * dvhV;
		if (fRepopulationFactor != 0)
		{
			OEDs += dvhV * exp(-alphaPrime * dvhD) / (alphaPrime * fRepopulationFactor) * (1 - 2 * fRepopulationFactor +
					exp(alphaPrime * dvhD) * pow(fRepopulationFactor, 2) - alphaPrime * fRepopulationFactor * dvhD -
					exp(-alphaPrime*fRepopulationFactor/(1-fRepopulationFactor)*dvhD) * pow(1-fRepopulationFactor, 2));
		}
		else
		{
			OEDs += dvhV * dvhD * exp(-alphaPrime * dvhD);
		}
	}
	meanDose /= sumVolume;
	JobResults.push_back(to_string(meanDose));
	return OEDs / sumVolume;
}

G4double TsRiskModel::GetOEDCarcinomaModel(std::vector<G4double> dose, std::vector<G4double> volume)
{
	G4double beta = fAlpha / fAlphaBeta;
	G4double sumVolume = 0;
	G4double OEDc = 0;
	G4double meanDose = 0;
	for (G4int i = 0; i < dose.size(); i++)
	{
		G4double dvhD = dose[i];
		G4double dvhV = volume[i];
		G4double alphaPrime = fAlpha + beta * dvhD / (G4double)fNumberOfFractions;
		sumVolume += dvhV;
		meanDose += dvhD * dvhV;
		if (fRepopulationFactor != 0)
		{
			OEDc += dvhV * exp(-alphaPrime * dvhD) / (alphaPrime * fRepopulationFactor) * (1 - 2 * fRepopulationFactor +
					exp(alphaPrime * dvhD) * pow(fRepopulationFactor, 2) - exp(-alphaPrime*fRepopulationFactor/(1-fRepopulationFactor)*dvhD) *
					pow(1-fRepopulationFactor, 2));
		}
		else
		{
			OEDc += dvhV * dvhD * exp(-alphaPrime * dvhD);
		}
	}
	meanDose /= sumVolume;
	JobResults.push_back(to_string(meanDose));
	return OEDc / sumVolume;
}

void TsRiskModel::ReadLifetimeRiskTable()
{
	G4String fileName = fSEERDirectory + "LifetimeRisks_Sa.csv";
	// File pointer
	fstream f;
	// Open csv file
	f.open(fileName, ios::in);

	if (f.fail())
	{
		G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << "Lifetime risks file not found in provided directory. Please check the parameter Sc/CancerRisk/DataDirectory." << G4endl;
		exit(1);
	}
	G4int count = 0;
	G4String line, value;
	vector<G4String> row;

	// Read file
	while (f >> line)
	{
		row.clear();
		stringstream s(line);
		while (getline(s, value, ',')) { row.push_back(value); }
		vAttainedAge.push_back(count);
		SaMale.push_back(stod(row[0]));
		SaFemale.push_back(stod(row[1]));
		count++;
	}
	// Close file
	f.close();
}

void TsRiskModel::ReadDVHcsv(G4String dvhFile)
{
    MeshGeomTools::readCSV(dvhFile, DVHDose, DVHVolume);
}

void TsRiskModel::ReadSEEROrganSpecificTable(G4String organ)
{
	// List of available organs
	vector<G4String> availableOrgans;
	availableOrgans.push_back("brain");
	availableOrgans.push_back("breast");
	availableOrgans.push_back("colonrectum");
	availableOrgans.push_back("esophagus");
	availableOrgans.push_back("kidneys");
	availableOrgans.push_back("larynx");
	availableOrgans.push_back("liver");
	availableOrgans.push_back("lung");
	availableOrgans.push_back("lungbronchus");
	availableOrgans.push_back("pancreas");
	availableOrgans.push_back("pharynx");
	availableOrgans.push_back("prostateovary");
	availableOrgans.push_back("stomach");
	availableOrgans.push_back("testesuterus");
	availableOrgans.push_back("thyroid");
	availableOrgans.push_back("urinarybladder");

	// Check if specified organ is in the list
	if (std::find(std::begin(availableOrgans), std::end(availableOrgans), organ) == std::end(availableOrgans))
	{
		G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << "There is not specific data for the specified organ (" + organ + ")." << G4endl;
		G4cerr << "Please use one of the following names for your OAR: " << G4endl;
		for (G4int i = 0; i < availableOrgans.size() - 1; i++)
		{
			G4cerr << availableOrgans[i] << ",";
		}
		G4cerr << availableOrgans[availableOrgans.size()-1] << "." << G4endl;
		exit(1);
	}

	SaOrganMale.clear();
	SaOrganFemale.clear();

	G4String fileName = fSEERDirectory + "SEER_" + organ + ".csv";

	fstream f;
	f.open(fileName, ios::in);
	if (f.fail())
	{
		G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << "SEER data file not found in provided directory. Please check the parameter Sc/CancerRisk/DataDirectory." << G4endl;
		exit(1);
	}
	G4String line, value;
	vector<G4String> row;
	// Read file
	while (f >> line)
	{
		row.clear();
		stringstream s(line);
		while (getline(s, value, ',')) { row.push_back(value); }
		SaOrganMale.push_back(stod(row[0]));
		SaOrganFemale.push_back(stod(row[1]));
	}
	// Close file
	f.close();
}

void TsRiskModel::SetOrganSpecificParameters()
{
    // Sets fAlpha, fAlphaBeta, fBetaERR, fBetaEAR, fGamma_e_ERR, fGamma_e_EAR, fGamma_a_ERR, fGamma_e_EAR,
    // fRepopulationFactor
    // Based on fSex, fAge, fOrganAtRisk

    // Open the input file
	G4String parameterFile;
    parameterFile = fSEERDirector + fParameterSet + ".txt";
    std::ifstream inputFile(parameterFile);
    
    if (!inputFile.is_open()) {
		G4cerr << "Could not open parameter table, no file: " << parameterFile << G4endl;
        exit(1);
    }

    // skip headers
    std::string line;
    while (std::getline(inputFile, line)) {
        if (line[0] != '#'){break;}
    }
    
    // Get column names
    std::string token;
    std::vector<std::string> column_names;
    std::stringstream linestream(line);
    while (linestream >> token) {
        column_names.push_back(token);
    }

    // Get rows
    int col_n = 0;
    int row_n = 0;
    std::map<G4String, G4String> row;
    while (std::getline(inputFile, line)) {
        std::stringstream ss(line);
        col_n = 0;
        while (ss >> token) {
            row[column_names[col_n]] = token;
            col_n++;
        }
        if (row["SEER_name"] == fOrganAtRisk) {break;}
        row_n++;
        if (row_n > (column_names.size() - 1)){
            G4cerr << "Organ: " << fOrganAtRisk << " not available in parameter file: " << parameterFile << G4endl;
            exit(1);
        }
    }
    // Close the input file
    inputFile.close();

    fAlpha              = GetMapDouble(row, "Alpha");
    fAlphaBeta          = GetMapDouble(row, "AlphaBeta");
    fRepopulationFactor = GetMapDouble(row, "RepopulationFactor");
    fRepopulationFactor = GetMapDouble(row, "RepopulationFactor");
    fParAttainedAge = GetMapDouble(row, "Attained_Age_Ref");
    if (fSex == "male"){
        fBetaERR            = GetMapDouble(row, "Beta_ERR_M");
        fBetaEAR            = GetMapDouble(row, "Beta_EAR_M");
    }
    else if (fSex == "female"){
        fBetaERR            = GetMapDouble(row, "Beta_ERR_F");
        fBetaEAR            = GetMapDouble(row, "Beta_EAR_F");
    }
    fGamma_e_ERR        = GetMapDouble(row, "Gamma_e_ERR");
    fGamma_e_EAR        = GetMapDouble(row, "Gamma_e_EAR");
    fGamma_a_ERR        = GetMapDouble(row, "Gamma_a_ERR");
    fGamma_a_EAR        = GetMapDouble(row, "Gamma_a_EAR");
}

G4double TsRiskModel::GetMapDouble(std::map<G4String, G4String> &table, G4String key){
    return (G4double) std::stod(table[key]);
}

void TsRiskModel::WriteOutputFile(G4String fileName, vector<vector<G4double>> quantitiesToWrite, vector<G4String> headers)
{
	fstream fout;
	// Creates a new file, overwritting if exists
	fout.open(fileName, ios::out | ios::trunc);
	// Insert headers
	for (vector<G4String>::iterator it = headers.begin(); it != headers.end(); ++it)
	{
		fout << *it;
		fout << ",";
	}
	fout << "\n";

	// Insert data
	G4int lengthQ = quantitiesToWrite[0].size();
	for (G4int iLine = 0; iLine < lengthQ; iLine++)
	{
		for (G4int iQ = 0; iQ < quantitiesToWrite.size(); iQ++)
		{
			fout << quantitiesToWrite[iQ][iLine];
			fout << ",";
		}
		fout << "\n";
	}
	fout.close();
}

void TsRiskModel::WriteOutputFile(G4String fileName, vector<G4String> quantitiesToWrite, vector<G4String> headers)
{
	fstream fout;
	// Creates a new file, overwritting if exists
	fout.open(fileName, ios::out | ios::trunc);
	// Insert headers
	for (vector<G4String>::iterator it = headers.begin(); it != headers.end(); ++it)
	{
		fout << *it;
		fout << ",";
	}
	fout << "\n";

	// Insert data
	for (G4int iQ = 0; iQ < quantitiesToWrite.size(); iQ++)
	{
		fout << quantitiesToWrite[iQ];
		fout << ",";
	}
	fout << "\n";

	fout.close();
}

void TsRiskModelWrapper(char* patient, char* cancerModel, char* sex, char* parameterSet,
                        char* organAtRisk, G4int ageAtExposure, char* SEERDirectory,
                        G4int attainedAge, G4int numberOfFractions, char* dvhFile)
{
    TsRiskModel riskModel(patient, cancerModel, sex, parameterSet,
                          organAtRisk, ageAtExposure, SEERDirectory, 
                          attainedAge, numberOfFractions, dvhFile);
}
