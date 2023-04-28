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
class XtraVarProd : public edm::stream::EDProducer<> {
  public:
    explicit XtraVarProd(const edm::ParameterSet &iConfig):
      src_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("src"))),
      pvsrc_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvsrc"))),
      svsrc_(consumes<edm::View<reco::VertexCompositePtrCandidate>>(iConfig.getParameter<edm::InputTag>("svsrc"))),
      rhosrc_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhosrc")))
    {
      //un prodotto da copiare
      produces<edm::ValueMap<int>>("numDaughtersPt03vec");

      produces<edm::ValueMap<float>>("emFractionEnergyRingsvec0");
      produces<edm::ValueMap<float>>("emFractionEnergyRingsvec1");
      produces<edm::ValueMap<float>>("emFractionEnergyRingsvec2");
      produces<edm::ValueMap<float>>("emFractionEnergyRingsvec3");
      produces<edm::ValueMap<float>>("emFractionEnergyRingsvec4");

      produces<edm::ValueMap<float>>("chFractionEnergyRingsvec0");
      produces<edm::ValueMap<float>>("chFractionEnergyRingsvec1");
      produces<edm::ValueMap<float>>("chFractionEnergyRingsvec2");
      produces<edm::ValueMap<float>>("chFractionEnergyRingsvec3");
      produces<edm::ValueMap<float>>("chFractionEnergyRingsvec4");

      produces<edm::ValueMap<float>>("muFractionEnergyRingsvec0");
      produces<edm::ValueMap<float>>("muFractionEnergyRingsvec1");
      produces<edm::ValueMap<float>>("muFractionEnergyRingsvec2");
      produces<edm::ValueMap<float>>("muFractionEnergyRingsvec3");
      produces<edm::ValueMap<float>>("muFractionEnergyRingsvec4");

      produces<edm::ValueMap<float>>("neFractionEnergyRingsvec0");
      produces<edm::ValueMap<float>>("neFractionEnergyRingsvec1");
      produces<edm::ValueMap<float>>("neFractionEnergyRingsvec2");
      produces<edm::ValueMap<float>>("neFractionEnergyRingsvec3");
      produces<edm::ValueMap<float>>("neFractionEnergyRingsvec4");

      // add

    }
    ~XtraVarProd() override {};
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    //void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const;
    void produce(edm::Event&, const edm::EventSetup&) override;
    void beginStream(edm::StreamID) override;
    void endStream() override;

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
XtraVarProd<T>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<edm::View<T>> src;
  iEvent.getByToken(src_, src);

  const auto& srcJet = iEvent.getHandle(src_);


  auto nJet = srcJet->size();

  std::vector<int> numDaughtersPt03vec(nJet, 0);

  std::vector<float> emFractionEnergyRingsvec0(nJet, 0);
  std::vector<float> emFractionEnergyRingsvec1(nJet, 0);
  std::vector<float> emFractionEnergyRingsvec2(nJet, 0);
  std::vector<float> emFractionEnergyRingsvec3(nJet, 0);
  std::vector<float> emFractionEnergyRingsvec4(nJet, 0);

  std::vector<float> chFractionEnergyRingsvec0(nJet, 0);
  std::vector<float> chFractionEnergyRingsvec1(nJet, 0);
  std::vector<float> chFractionEnergyRingsvec2(nJet, 0);
  std::vector<float> chFractionEnergyRingsvec3(nJet, 0);
  std::vector<float> chFractionEnergyRingsvec4(nJet, 0);

  std::vector<float> neFractionEnergyRingsvec0(nJet, 0);
  std::vector<float> neFractionEnergyRingsvec1(nJet, 0);
  std::vector<float> neFractionEnergyRingsvec2(nJet, 0);
  std::vector<float> neFractionEnergyRingsvec3(nJet, 0);
  std::vector<float> neFractionEnergyRingsvec4(nJet, 0);

  std::vector<float> muFractionEnergyRingsvec0(nJet, 0);
  std::vector<float> muFractionEnergyRingsvec1(nJet, 0);
  std::vector<float> muFractionEnergyRingsvec2(nJet, 0);
  std::vector<float> muFractionEnergyRingsvec3(nJet, 0);
  std::vector<float> muFractionEnergyRingsvec4(nJet, 0);

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

    emFractionEnergyRingsvec0[counter] = emFractionEnergyRings[0];
    emFractionEnergyRingsvec1[counter] = emFractionEnergyRings[1];
    emFractionEnergyRingsvec2[counter] = emFractionEnergyRings[2];
    emFractionEnergyRingsvec3[counter] = emFractionEnergyRings[3];
    emFractionEnergyRingsvec4[counter] = emFractionEnergyRings[4];

    muFractionEnergyRingsvec0[counter] = muFractionEnergyRings[0];
    muFractionEnergyRingsvec1[counter] = muFractionEnergyRings[1];
    muFractionEnergyRingsvec2[counter] = muFractionEnergyRings[2];
    muFractionEnergyRingsvec3[counter] = muFractionEnergyRings[3];
    muFractionEnergyRingsvec4[counter] = muFractionEnergyRings[4];

    chFractionEnergyRingsvec0[counter] = chFractionEnergyRings[0];
    chFractionEnergyRingsvec1[counter] = chFractionEnergyRings[1];
    chFractionEnergyRingsvec2[counter] = chFractionEnergyRings[2];
    chFractionEnergyRingsvec3[counter] = chFractionEnergyRings[3];
    chFractionEnergyRingsvec4[counter] = chFractionEnergyRings[4];

    neFractionEnergyRingsvec0[counter] = neFractionEnergyRings[0];
    neFractionEnergyRingsvec1[counter] = neFractionEnergyRings[1];
    neFractionEnergyRingsvec2[counter] = neFractionEnergyRings[2];
    neFractionEnergyRingsvec3[counter] = neFractionEnergyRings[3];
    neFractionEnergyRingsvec4[counter] = neFractionEnergyRings[4];


  }

  std::unique_ptr<edm::ValueMap<int>> numDaughtersPt03_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler_duaghters(*numDaughtersPt03_VM);
  filler_duaghters.insert(srcJet, numDaughtersPt03vec.begin(), numDaughtersPt03vec.end());
  filler_duaghters.fill();
  iEvent.put(std::move(numDaughtersPt03_VM), "numDaughtersPt03vec");


  std::unique_ptr<edm::ValueMap<float>> emFractionEnergyRingsvec0_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_emFractionEnergyRings0(*emFractionEnergyRingsvec0_VM);
  filler_emFractionEnergyRings0.insert(srcJet, emFractionEnergyRingsvec0.begin(), emFractionEnergyRingsvec0.end());
  filler_emFractionEnergyRings0.fill();
  iEvent.put(std::move(emFractionEnergyRingsvec0_VM), "emFractionEnergyRingsvec0");
  
  std::unique_ptr<edm::ValueMap<float>> emFractionEnergyRingsvec1_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_emFractionEnergyRings1(*emFractionEnergyRingsvec1_VM);
  filler_emFractionEnergyRings1.insert(srcJet, emFractionEnergyRingsvec1.begin(), emFractionEnergyRingsvec1.end());
  filler_emFractionEnergyRings1.fill();
  iEvent.put(std::move(emFractionEnergyRingsvec1_VM), "emFractionEnergyRingsvec1");

  std::unique_ptr<edm::ValueMap<float>> emFractionEnergyRingsvec2_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_emFractionEnergyRings2(*emFractionEnergyRingsvec2_VM);
  filler_emFractionEnergyRings2.insert(srcJet, emFractionEnergyRingsvec2.begin(), emFractionEnergyRingsvec2.end());
  filler_emFractionEnergyRings2.fill();
  iEvent.put(std::move(emFractionEnergyRingsvec2_VM), "emFractionEnergyRingsvec2");

  std::unique_ptr<edm::ValueMap<float>> emFractionEnergyRingsvec3_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_emFractionEnergyRings3(*emFractionEnergyRingsvec3_VM);
  filler_emFractionEnergyRings3.insert(srcJet, emFractionEnergyRingsvec3.begin(), emFractionEnergyRingsvec3.end());
  filler_emFractionEnergyRings3.fill();
  iEvent.put(std::move(emFractionEnergyRingsvec3_VM), "emFractionEnergyRingsvec3");

  std::unique_ptr<edm::ValueMap<float>> emFractionEnergyRingsvec4_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_emFractionEnergyRings4(*emFractionEnergyRingsvec4_VM);
  filler_emFractionEnergyRings4.insert(srcJet, emFractionEnergyRingsvec4.begin(), emFractionEnergyRingsvec4.end());
  filler_emFractionEnergyRings4.fill();
  iEvent.put(std::move(emFractionEnergyRingsvec4_VM), "emFractionEnergyRingsvec4");


  std::unique_ptr<edm::ValueMap<float>> muFractionEnergyRingsvec0_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_muFractionEnergyRings0(*muFractionEnergyRingsvec0_VM);
  filler_muFractionEnergyRings0.insert(srcJet, muFractionEnergyRingsvec0.begin(), muFractionEnergyRingsvec0.end());
  filler_muFractionEnergyRings0.fill();
  iEvent.put(std::move(muFractionEnergyRingsvec0_VM), "muFractionEnergyRingsvec0");
  
  std::unique_ptr<edm::ValueMap<float>> muFractionEnergyRingsvec1_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_muFractionEnergyRings1(*muFractionEnergyRingsvec1_VM);
  filler_muFractionEnergyRings1.insert(srcJet, muFractionEnergyRingsvec1.begin(), muFractionEnergyRingsvec1.end());
  filler_muFractionEnergyRings1.fill();
  iEvent.put(std::move(muFractionEnergyRingsvec1_VM), "muFractionEnergyRingsvec1");

  std::unique_ptr<edm::ValueMap<float>> muFractionEnergyRingsvec2_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_muFractionEnergyRings2(*muFractionEnergyRingsvec2_VM);
  filler_muFractionEnergyRings2.insert(srcJet, muFractionEnergyRingsvec2.begin(), muFractionEnergyRingsvec2.end());
  filler_muFractionEnergyRings2.fill();
  iEvent.put(std::move(muFractionEnergyRingsvec2_VM), "muFractionEnergyRingsvec2");

  std::unique_ptr<edm::ValueMap<float>> muFractionEnergyRingsvec3_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_muFractionEnergyRings3(*muFractionEnergyRingsvec3_VM);
  filler_muFractionEnergyRings3.insert(srcJet, muFractionEnergyRingsvec3.begin(), muFractionEnergyRingsvec3.end());
  filler_muFractionEnergyRings3.fill();
  iEvent.put(std::move(muFractionEnergyRingsvec3_VM), "muFractionEnergyRingsvec3");

  std::unique_ptr<edm::ValueMap<float>> muFractionEnergyRingsvec4_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_muFractionEnergyRings4(*muFractionEnergyRingsvec4_VM);
  filler_muFractionEnergyRings4.insert(srcJet, muFractionEnergyRingsvec4.begin(), muFractionEnergyRingsvec4.end());
  filler_muFractionEnergyRings4.fill();
  iEvent.put(std::move(muFractionEnergyRingsvec4_VM), "muFractionEnergyRingsvec4");


  std::unique_ptr<edm::ValueMap<float>> chFractionEnergyRingsvec0_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_chFractionEnergyRings0(*chFractionEnergyRingsvec0_VM);
  filler_chFractionEnergyRings0.insert(srcJet, chFractionEnergyRingsvec0.begin(), chFractionEnergyRingsvec0.end());
  filler_chFractionEnergyRings0.fill();
  iEvent.put(std::move(chFractionEnergyRingsvec0_VM), "chFractionEnergyRingsvec0");
  
  std::unique_ptr<edm::ValueMap<float>> chFractionEnergyRingsvec1_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_chFractionEnergyRings1(*chFractionEnergyRingsvec1_VM);
  filler_chFractionEnergyRings1.insert(srcJet, chFractionEnergyRingsvec1.begin(), chFractionEnergyRingsvec1.end());
  filler_chFractionEnergyRings1.fill();
  iEvent.put(std::move(chFractionEnergyRingsvec1_VM), "chFractionEnergyRingsvec1");

  std::unique_ptr<edm::ValueMap<float>> chFractionEnergyRingsvec2_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_chFractionEnergyRings2(*chFractionEnergyRingsvec2_VM);
  filler_chFractionEnergyRings2.insert(srcJet, chFractionEnergyRingsvec2.begin(), chFractionEnergyRingsvec2.end());
  filler_chFractionEnergyRings2.fill();
  iEvent.put(std::move(chFractionEnergyRingsvec2_VM), "chFractionEnergyRingsvec2");

  std::unique_ptr<edm::ValueMap<float>> chFractionEnergyRingsvec3_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_chFractionEnergyRings3(*chFractionEnergyRingsvec3_VM);
  filler_chFractionEnergyRings3.insert(srcJet, chFractionEnergyRingsvec3.begin(), chFractionEnergyRingsvec3.end());
  filler_chFractionEnergyRings3.fill();
  iEvent.put(std::move(chFractionEnergyRingsvec3_VM), "chFractionEnergyRingsvec3");

  std::unique_ptr<edm::ValueMap<float>> chFractionEnergyRingsvec4_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_chFractionEnergyRings4(*chFractionEnergyRingsvec4_VM);
  filler_chFractionEnergyRings4.insert(srcJet, chFractionEnergyRingsvec4.begin(), chFractionEnergyRingsvec4.end());
  filler_chFractionEnergyRings4.fill();
  iEvent.put(std::move(chFractionEnergyRingsvec4_VM), "chFractionEnergyRingsvec4");


  std::unique_ptr<edm::ValueMap<float>> neFractionEnergyRingsvec0_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_neFractionEnergyRings0(*neFractionEnergyRingsvec0_VM);
  filler_neFractionEnergyRings0.insert(srcJet, neFractionEnergyRingsvec0.begin(), neFractionEnergyRingsvec0.end());
  filler_neFractionEnergyRings0.fill();
  iEvent.put(std::move(neFractionEnergyRingsvec0_VM), "neFractionEnergyRingsvec0");
  
  std::unique_ptr<edm::ValueMap<float>> neFractionEnergyRingsvec1_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_neFractionEnergyRings1(*neFractionEnergyRingsvec1_VM);
  filler_neFractionEnergyRings1.insert(srcJet, neFractionEnergyRingsvec1.begin(), neFractionEnergyRingsvec1.end());
  filler_neFractionEnergyRings1.fill();
  iEvent.put(std::move(neFractionEnergyRingsvec1_VM), "neFractionEnergyRingsvec1");

  std::unique_ptr<edm::ValueMap<float>> neFractionEnergyRingsvec2_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_neFractionEnergyRings2(*neFractionEnergyRingsvec2_VM);
  filler_neFractionEnergyRings2.insert(srcJet, neFractionEnergyRingsvec2.begin(), neFractionEnergyRingsvec2.end());
  filler_neFractionEnergyRings2.fill();
  iEvent.put(std::move(neFractionEnergyRingsvec2_VM), "neFractionEnergyRingsvec2");

  std::unique_ptr<edm::ValueMap<float>> neFractionEnergyRingsvec3_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_neFractionEnergyRings3(*neFractionEnergyRingsvec3_VM);
  filler_neFractionEnergyRings3.insert(srcJet, neFractionEnergyRingsvec3.begin(), neFractionEnergyRingsvec3.end());
  filler_neFractionEnergyRings3.fill();
  iEvent.put(std::move(neFractionEnergyRingsvec3_VM), "neFractionEnergyRingsvec3");

  std::unique_ptr<edm::ValueMap<float>> neFractionEnergyRingsvec4_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_neFractionEnergyRings4(*neFractionEnergyRingsvec4_VM);
  filler_neFractionEnergyRings4.insert(srcJet, neFractionEnergyRingsvec4.begin(), neFractionEnergyRingsvec4.end());
  filler_neFractionEnergyRings4.fill();
  iEvent.put(std::move(neFractionEnergyRingsvec4_VM), "neFractionEnergyRingsvec4");


  // edm::LogWarning("DEBUGGER") << "Produce Done";


}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
template <typename T>
void
XtraVarProd<T>::beginStream(edm::StreamID) 
{
  //edm::LogWarning("DEBUGGER") << "beginStream";
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
template <typename T>
void
XtraVarProd<T>::endStream() {
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
