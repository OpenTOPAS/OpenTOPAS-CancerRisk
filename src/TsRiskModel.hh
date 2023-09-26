#ifndef TsRiskModel_hh
#define TsRiskModel_hh

#include "globals.hh"

#include "TsVOutcomeModel.hh"
#include <vector>

// To implement as a TOPAS class
class TsRiskModel : public TsVOutcomeModel
{
public:
	TsRiskModel(TsParameterManager* pm, G4String parmName);
	TsRiskModel(G4String patient, G4String cancerModel, G4String fSex,
                G4String parameterSet, G4String organAtRisk, G4int ageAtExposure, G4String SEERDirectory,
                G4int attainedAge, G4int numberOfFractions, G4String dvhFile);
	virtual ~TsRiskModel();

	// Mandatory method from parent class
	G4double Initialize(std::vector<G4double> dose, std::vector<G4double> volume);
	void ResolveParameters();

	// Specific methods for this class
	std::vector<G4double> LARBasedOnERR(G4double OED);
	std::vector<G4double> LARBasedOnEAR(G4double OED);
	G4double AverageLAReAndLARa(std::vector<G4double> LARERRCum, std::vector<G4double> LAREARCum);
	G4double GetOEDSarcomaModel(std::vector<G4double> dose, std::vector<G4double> volume);
	G4double GetOEDCarcinomaModel(std::vector<G4double> dose, std::vector<G4double> volume);
	void ReadLifetimeRiskTable();
	void ReadDVHcsv(G4String dvhFile);
	void ReadSEEROrganSpecificTable(G4String organ);
	void SetOrganSpecificParameters();
	void WriteOutputFile(G4String fileName, std::vector<std::vector<G4double>> quantitiesToWrite, std::vector<G4String> headers);
	void WriteOutputFile(G4String fileName, std::vector<G4String> quantitiesToWrite, std::vector<G4String> headers);

	// Not necessary for TOPAS class
	void inline SetOrganAtRisk(G4String organ)	{ fOrganAtRisk = organ; }
	void inline SetParameterSet(G4String pSet)	{ fParameterSet = pSet; }
	void inline SetSex(G4String sex)			{ fSex = sex; }
	void inline SetAttainedAge(G4int age)		{ fAttainedAge = age; }
	void inline SetAgeAtExposure(G4int age)		{ fAgeAtExposure = age; }
	void inline SetPatientName(G4String pat)	{ fPatient = pat; }
	void inline SetNumberOfFractions(G4int n)	{ fNumberOfFractions = n; }
	void inline SetCancerModel(G4String mod)	{ fCancerModel = mod; }


private:
	TsParameterManager* fPm;
	G4String fModelName;
	G4String fSEERDirectory;

	G4String fPatient, fSex, fParameterSet, fCancerModel;
	G4String fOrganAtRisk;
	G4int fParAttainedAge, fLatency;
	G4int fAttainedAge, fAgeAtExposure;
	G4int fNumberOfFractions;
	G4double fDDREF;
	G4double fAlpha, fAlphaBeta, fRepopulationFactor;
	G4double fBetaERR, fBetaEAR;
	G4double fGamma_e_ERR, fGamma_a_ERR, fGamma_e_EAR, fGamma_a_EAR;

	std::vector<G4double> SaMale;
	std::vector<G4double> SaFemale;
	std::vector<G4int> vAttainedAge;

	std::vector<G4double> SaOrganMale;
	std::vector<G4double> SaOrganFemale;

	std::vector<G4String> JobHeaders;
	std::vector<G4String> JobResults;

	std::vector<G4double> DVHDose;
	std::vector<G4double> DVHVolume;
};

extern "C" void TsRiskModelWrapper(char* patient, char* cancerModel, char* fSex,
                                   char* parameterSet, char* organAtRisk, G4int ageAtExposure,
                                   char* SEERDirectory, G4int attainedAge, G4int numberOfFractions,
                                   char* dvhFile);
#endif
