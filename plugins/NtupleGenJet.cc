// -*- C++ -*-
//
// Package:    validation/NtupleGenJet
// Class:      NtupleGenJet
// 
/**\class NtupleGenJet NtupleGenJet.cc validation/NtupleGenJet/plugins/NtupleGenJet.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lata Panwar
//         Created:  Wed, 14 Mar 2018 12:41:16 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "TTree.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TMath.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class NtupleGenJet : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit NtupleGenJet(const edm::ParameterSet&);
  ~NtupleGenJet();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  
  // ----------member data ---------------------------
  edm::EDGetTokenT <reco::GenParticleCollection> genparticlesToken;
	
  int genHiggs_n_=0;
  int genY_n_=0;
  TH1D *nHiggs_histo;
  TH1F *hpt_b;
  TH1F *heta_b;
  TH1F *hphi_b;
  TH1F *hpt_photon;
  TH1F *heta_photon;
  TH1F *hphi_photon;
  TH1F *hpt_y;
  TH1F *heta_y;
  TH1F *hphi_y;
  TH1F *hm_y;
  TH1F *hn_y;
  TH1F *hpt_higgs;
  TH1F *heta_higgs;
  TH1F *hphi_higgs;
  TH1F *hm_higgs;
  TH1F *hn_higgs;
  TH1F *hcostheta_hh_cs;
  TH1F *hm_bbaa;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
NtupleGenJet::NtupleGenJet(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   genparticlesToken = consumes <reco::GenParticleCollection> (std::string("genParticles"));
   edm::Service<TFileService> fs; 
   hn_higgs = fs->make<TH1F>("N_higgs" , ";N_{H};Events;;" , 50 , 0 , 50 );
   hpt_b = fs->make<TH1F>("pt_b" , ";p_{T} of b[GeV];Events;;" , 100 , 0 , 500 );
   heta_b=fs->make<TH1F>("eta_b" , ";#eta of b;Events;;" , 20 , -5 , 5 );
   hphi_b=fs->make<TH1F>("phi_b" , ";#phi of b;Events;;" , 20 , -5 , 5 );
   hpt_higgs= fs->make<TH1F>("pt_higgs" , ";p_{T} of H[GeV];Events;;" , 100 , 0 , 500 );
   heta_higgs=fs->make<TH1F>("eta_higgs" , ";#eta of H;Events;;" , 20 , -5 , 5 );
   hphi_higgs=fs->make<TH1F>("phi_higgs" , ";#phi of H;Events;;" , 20 , -5 , 5 );
   hm_higgs=fs->make<TH1F>("m_higgs" , ";mass of Higgs[GeV];Events;;" , 100,0,500);
   hpt_photon = fs->make<TH1F>("pt_photon" , ";p_{T} of photon[GeV];Events;;" , 100 , 0 , 500 );
   heta_photon=fs->make<TH1F>("eta_photon" , ";#eta of photon;Events;;" , 20 , -5 , 5 );
   hphi_photon=fs->make<TH1F>("phi_photon" , ";#phi of photon;Events;;" , 20 , -5 , 5 );
   hpt_y= fs->make<TH1F>("pt_y" , ";p_{T} of Y[GeV];Events;;" , 100 , 0 , 500 );
   heta_y=fs->make<TH1F>("eta_y" , ";#eta of Y;Events;;" , 20 , -5 , 5 );
   hphi_y=fs->make<TH1F>("phi_y" , ";#phi of Y;Events;;" , 20 , -5 , 5 );
   hm_y=fs->make<TH1F>("m_y" , ";mass of Y[GeV];Events;;" , 100,0,500);
   hn_y = fs->make<TH1F>("N_y" , ";N_{Y};Events;;" , 20 , 0 , 50 );
   hcostheta_hh_cs = fs->make<TH1F>("costheta",";|cos#theta_{HY/HH}^{CS}|;Events;;", 10, 0, 1);
   hm_bbaa=fs->make<TH1F>("m_bbaa","", 1000, 299.5, 300.5);
}



NtupleGenJet::~NtupleGenJet()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
NtupleGenJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genparticlesToken, genParticles);
  //for(reco::GenParticle jet : *(gen_h.product())){
  // for(const auto& jet : genparticles){
  TLorentzVector h1, h2;
  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const reco::GenParticle & p = (*genParticles)[i];
    int id = p.pdgId();
    //    int st = p.status();
    //     double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
    if (id == 25){
      int n = p.numberOfDaughters();
      if(n < 2 ) continue;
      const reco::Candidate * d1 = p.daughter(0);
      const reco::Candidate * d2 = p.daughter(1);
      if (std::abs(d1->pdgId())==22 && std::abs(d2->pdgId())==22){
	TLorentzVector v1, v2;
	v1.SetPtEtaPhiM(d1->pt(),d1->eta(),d1->phi(),d1->mass());
	v2.SetPtEtaPhiM(d2->pt(),d2->eta(),d2->phi(),d2->mass());
	h1 = v1+v2;
	++genHiggs_n_;
	hn_higgs->Fill(genHiggs_n_);
	hpt_higgs->Fill(p.pt());
	heta_higgs->Fill(p.eta());
	hphi_higgs ->Fill(p.phi());
	hm_higgs->Fill(p.mass());
	//// filling for photon
	hpt_photon->Fill(d1->pt());
	heta_photon->Fill(d1->eta());
	hphi_photon->Fill(d1->phi());
	
	hpt_photon->Fill(d2->pt());
	heta_photon->Fill(d2->eta());
	hphi_photon->Fill(d2->phi());
      }
    }
    if (id == 35){
      int n = p.numberOfDaughters();
      if(n < 2 ) continue;
      const reco::Candidate * d1 = p.daughter(0);
      const reco::Candidate * d2 = p.daughter(1);
      if (std::abs(d1->pdgId())==5 && std::abs(d2->pdgId())==5){
	TLorentzVector v1, v2;
	v1.SetPtEtaPhiM(d1->pt(),d1->eta(),d1->phi(),d1->mass());
	v2.SetPtEtaPhiM(d2->pt(),d2->eta(),d2->phi(),d2->mass());
	h2 = v1+v2;
	++genY_n_;
	hn_y->Fill(genY_n_);
	hpt_y->Fill(p.pt());
	heta_y->Fill(p.eta());
	hphi_y ->Fill(p.phi());
	hm_y->Fill(p.mass());
	//// filling for b's
	hpt_b->Fill(d1->pt());
	heta_b->Fill(d1->eta());
	hphi_b->Fill(d1->phi());
	
	hpt_b->Fill(d2->pt());
	heta_b->Fill(d2->eta());
	hphi_b->Fill(d2->phi());
      }
    }
  }
  TLorentzVector hh= h1+h2;
  h1.Boost(-hh.BoostVector());
  hcostheta_hh_cs->Fill(TMath::Abs(h1.CosTheta()));
  hm_bbaa->Fill(hh.M());
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NtupleGenJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(NtupleGenJet);
