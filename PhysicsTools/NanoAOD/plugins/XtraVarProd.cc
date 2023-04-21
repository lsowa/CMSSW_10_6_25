// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      XtraVarProd
// 
/**\class XtraVarProd XtraVarProd.cc PhysicsTools/NanoAOD/plugins/XtraVarProd.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Peruzzi
//         Created:  Tue, 05 Sep 2017 12:24:38 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"


#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "PhysicsTools/NanoAOD/interface/MatchingUtils.h"


//
// class declaration
//

template <typename T>
class XtraVarProd : public edm::global::EDProducer<> {
  public:
    explicit XtraVarProd(const edm::ParameterSet &iConfig):
      src_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("src"))),
      pvsrc_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvsrc"))),
      svsrc_(consumes<edm::View<reco::VertexCompositePtrCandidate>>(iConfig.getParameter<edm::InputTag>("svsrc"))),
      rhosrc_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhosrc")))
    {
      //un prodotto da copiare
      produces<edm::ValueMap<int>>("numDaughtersPt03vec");
      // add
    }
    ~XtraVarProd() override {};
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

    // ----------member data ---------------------------

    edm::EDGetTokenT<edm::View<pat::Jet>> src_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> pvsrc_;
    edm::EDGetTokenT<edm::View<reco::VertexCompositePtrCandidate>> svsrc_;
    edm::EDGetTokenT<double> rhosrc_;

    edm::Handle<std::vector<reco::Vertex>> pvs_;
    edm::Handle<edm::View<reco::VertexCompositePtrCandidate>> svs_;
    edm::Handle<double> rhos_;
 
   void readAdditionalCollections(edm::Event& iEvent, const edm::EventSetup&) {
    iEvent.getByToken(pvsrc_, pvs_);
    iEvent.getByToken(svsrc_, svs_);
    iEvent.getByToken(rhosrc_, rhos_);
  }

};


//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// member functions
//


// ------------ method called to produce the data  ------------
template <typename T>
void
XtraVarProd<T>::produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const{

  edm::Handle<edm::View<T>> src;
  iEvent.getByToken(src_, src);
  //readAdditionalCollections(iEvent, iSetup);
  //readAdditionalCollections(iEvent, iSetup);

  //iEvent.getByToken(pvsrc_, pvsrc_);
  //iEvent.getByToken(svsrc_, svsrc_);
  //iEvent.getByToken(rhosrc_, rhosrc_);
  const auto& srcJet = iEvent.getHandle(src_);
  //const auto& vtxProd = iEvent.get(srcVtx_);
  //const auto& svProd = iEvent.get(srcSV_);
  auto nJet = srcJet->size();
  std::vector<int> numDaughtersPt03vec(nJet, 0);

  int counter = -1;
  for (auto const& j : *src) { 
    counter += 1;

    float cone_boundaries[] = {0.05, 0.1, 0.2, 0.3, 0.4};
    size_t ncone_boundaries = sizeof(cone_boundaries) / sizeof(float);
    std::vector<float> emFractionEnergyRings(ncone_boundaries + 1);
    std::vector<float> chFractionEnergyRings(ncone_boundaries + 1);
    std::vector<float> neFractionEnergyRings(ncone_boundaries + 1);
    std::vector<float> muFractionEnergyRings(ncone_boundaries + 1);
    float jetRawEnergy = j.p4().E() * j.jecFactor("Uncorrected");
    int numDaughtersPt03 = 0;
    for (unsigned int ijcone = 0; ijcone < ncone_boundaries; ijcone++) {
      emFractionEnergyRings[ijcone] = 0;
      muFractionEnergyRings[ijcone] = 0;
      chFractionEnergyRings[ijcone] = 0;
      neFractionEnergyRings[ijcone] = 0;
    }
    for (const auto& d : j.daughterPtrVector()) {
      float candDr = Geom::deltaR(d->p4(), j.p4());
      size_t icone =
          std::lower_bound(&cone_boundaries[0], &cone_boundaries[ncone_boundaries], candDr) - &cone_boundaries[0];
      float candEnergy = d->energy() / jetRawEnergy;
      int pdgid = abs(d->pdgId());
      if (pdgid == 22 || pdgid == 11) {
        emFractionEnergyRings[icone] += candEnergy;
      } else if (pdgid == 13) {
        muFractionEnergyRings[icone] += candEnergy;
      } else if (d->charge() != 0) {
        chFractionEnergyRings[icone] += candEnergy;
      } else {
        neFractionEnergyRings[icone] += candEnergy;
      }
      if (d->pt() > 0.3)
        numDaughtersPt03 += 1;
    }  // end of jet daughters loop

    numDaughtersPt03vec[counter] = numDaughtersPt03;

  }

  std::unique_ptr<edm::ValueMap<int>> numDaughtersPt03_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler fillerRel(*numDaughtersPt03_VM);
  fillerRel.insert(srcJet, numDaughtersPt03vec.begin(), numDaughtersPt03vec.end());
  fillerRel.fill();
  iEvent.put(std::move(numDaughtersPt03_VM), "numDaughtersPt03vec");

}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template <typename T>
void
XtraVarProd<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src")->setComment("jet input collection");
  desc.add<edm::InputTag>("rhosrc")->setComment("rho source");
  desc.add<edm::InputTag>("pvsrc")->setComment("primary vertex input collection");
  desc.add<edm::InputTag>("svsrc")->setComment("secondary vertex input collection");
  
  std::string modname;
  if (typeid(T) == typeid(pat::Jet)) modname+="Jet";
  modname+="XtraVarProd_type";
  descriptions.add(modname,desc);
}

typedef XtraVarProd<pat::Jet> XtraVarProd_type;

//define this as a plug-in
DEFINE_FWK_MODULE(XtraVarProd_type);
