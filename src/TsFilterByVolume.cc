// Filter for OnlyIncludeIfParticleInVolume,OnlyIncludeIfParticleNotInVolume

#include "TsFilterByVolume.hh"

#include "G4VSolid.hh"
#include <iostream>

TsFilterByVolume::TsFilterByVolume(G4String name, TsParameterManager* pM,
                                   TsMaterialManager* mM, TsGeometryManager* gM,
                                   TsFilterManager* fM, TsVGenerator* generator,
                                   TsVScorer* scorer, TsVFilter* parentFilter)
:TsVFilter(name, pM, mM, gM, fM, generator, scorer, parentFilter)
{
    ResolveParameters();
    if (name=="onlyincludeifparticlenotinvolume") fInvert = true;
}


TsFilterByVolume::~TsFilterByVolume()
{;}


void TsFilterByVolume::ResolveParameters()
{
    fNames = fPm->GetStringVector(GetFullParmName(GetName()));
    fNamesLength = fPm->GetVectorLength(GetFullParmName(GetName()));
    CacheGeometryPointers();
}


void TsFilterByVolume::CacheGeometryPointers() {
	fSolids.clear();

	for (G4int i = 0; i < fNamesLength; i++) {
        G4VPhysicalVolume* volume = GetPhysicalVolume(fNames[i]);
        if (volume) {
            G4VSolid* solid = volume->GetLogicalVolume()->GetSolid();
            fSolids.push_back(solid);
        } 
        else if (volume->GetParameterisation()) {
            G4cerr << "Volume filter not yet implemented for parameterized volumes." << G4endl;
            fPm->AbortSession(1);
        }
        else  {
            G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
            G4cerr << GetName() << " = " << fNames[i] << " refers to an unknown Volume." << G4endl;
            fPm->AbortSession(1);
        }
	}
}

G4bool TsFilterByVolume::Accept(const G4Step* aStep) const {
	if (fParentFilter && !fParentFilter->Accept(aStep)) return false;

    // Loop through solids and return true if particle point is inside any
	for (G4int i = 0; i < (int)fSolids.size(); i++) {
        // Position is a const G4ThreeVector
        if (fSolids[i]->Inside(aStep->GetTrack()->GetStep()->GetPreStepPoint()->GetPosition())) {
            if (fInvert) return false;
            else return true;
        }
    }

	if (fInvert) return true;
	else return false;
}


G4bool TsFilterByVolume::AcceptTrack(const G4Track*) const {
	G4cerr << "Topas is exiting due to a serious error in source setup." << G4endl;
	G4cerr << "Sources cannot be filtered by " << GetName() << G4endl;
	fPm->AbortSession(1);
	return false;
}

