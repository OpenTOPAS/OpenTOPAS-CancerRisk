// Filter by the current volume of the particle
// May be roughly equivalent to "OnlyIncludeIfParticleLastInteractedInVolume" but should allow for 
// a list of possible volumes and avoid iterating over previous interacted volumes.
// Currently no methods for including ancestors or children, volumes must be listed explicitly
// 
// Author: Isaac Meyer imeyer@mgh.harvard.edu

#ifndef TsFilterByVolume_hh
#define TsFilterByVolume_hh

#include "TsVFilter.hh"

class TsFilterByVolume : public TsVFilter
{
    public:
        TsFilterByVolume(G4String name, TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM,
                TsFilterManager* fM, TsVGenerator* generator, TsVScorer* scorer, TsVFilter* parentFilter);
        virtual ~TsFilterByVolume();

        void ResolveParameters();
        void CacheGeometryPointers();

        virtual G4bool Accept(const G4Step*) const;
        virtual G4bool AcceptTrack(const G4Track*) const;

    private:
        G4String* fNames;
        G4int fNamesLength;
        std::vector<G4VSolid*> fSolids;
};

#endif
