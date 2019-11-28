#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

constexpr Double_t MFTLayerZ[] = {-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76.1, -77.5};

constexpr Double_t pMax = 4;
constexpr Double_t pMin = 0.1;
constexpr Double_t etaMin = -3.9;
constexpr Double_t etaMax = -2.0;

void MFT_Tracker_inspector(const Char_t *SimFile = "o2sim.root", const Char_t *trkFile = "mfttracks.root") {

  using o2::itsmft::Hit;
  using o2::MCTrackT;
  using o2::mft::TrackCA;
  using o2::mft::TrackLTF;
  using eventFoundTracks = std::vector<bool>;

  using trackHasHitsinDisks = std::array<bool,5>; // Disks with hits from a MFT track

  std::unique_ptr<TH1F> MCTrackspT = std::make_unique<TH1F> ("MC Tracks pT", "MC Tracks pT", 100, 0, pMax);
  std::unique_ptr<TH1F> MCTracksp = std::make_unique<TH1F> ("MC Tracks p", "MC Tracks p", 100, 0, pMax);
  std::unique_ptr<TH1F> MCTrackRap = std::make_unique<TH1F> ("MC Tracks eta", "MC Tracks Rapidity", 100, etaMin, etaMax);

  std::unique_ptr<TH1F> MFTTrackspT = std::make_unique<TH1F> ("MFT Tracks pT", "MFT Tracks pT", 100, 0, pMax);
  std::unique_ptr<TH1F> MFTTracksp = std::make_unique<TH1F> ("MFT Tracks p", "MFT Tracks p", 100, 0, pMax);
  std::unique_ptr<TH1F> MFTTrackRap = std::make_unique<TH1F> ("MFT Tracks eta", "MFT Tracks Rapidity", 100, etaMin, etaMax);

  std::unique_ptr<TH1F> LTFTrackspT = std::make_unique<TH1F> ("LTF Tracks pT", "LTF Tracks pT", 100, 0, 5);
  std::unique_ptr<TH1F> LTFTracksp = std::make_unique<TH1F> ("LTF Tracks p", "LTF Tracks p", 100, 0, 5);
  std::unique_ptr<TH1F> LTFTrackRap = std::make_unique<TH1F> ("LTF Tracks eta", "LTF Tracks Rapidity", 100, -4.0, -2);

  std::unique_ptr<TH1F> CATrackspT = std::make_unique<TH1F> ("CA Tracks pT", "CA Tracks pT", 100, 0, 5);
  std::unique_ptr<TH1F> CATracksp = std::make_unique<TH1F> ("CA Tracks p", "CA Tracks p", 100, 0, 5);
  std::unique_ptr<TH1F> CATrackRap = std::make_unique<TH1F> ("CA Tracks eta", "LTF Tracks Rapidity", 100, -4.0, -2);

  std::unique_ptr<TH1F> MCTracksEta5 = std::make_unique<TH1F> ("MC Tracks 5 eta", "-5 cm < zVertex < 5 cm", 100, etaMin, etaMax);
  std::unique_ptr<TH1F> MCTracksEta5_10pos = std::make_unique<TH1F> ("MC Tracks -5 -10 eta", "-10 cm < zVertex < -5 cm", 100, etaMin, etaMax);
  std::unique_ptr<TH1F> MCTracksEta5_10neg = std::make_unique<TH1F> ("MC Tracks 5 10 eta", "5 cm < zVertex < 10 cm", 100, etaMin, etaMax);

  std::unique_ptr<TH1F> MFTTracksEta5 = std::make_unique<TH1F> ("MFT Tracks 5 eta", "-5 cm < zVertex < 5 cm", 100, etaMin, etaMax);
  std::unique_ptr<TH1F> MFTTracksEta5_10pos = std::make_unique<TH1F> ("MFT Tracks -5 -10 eta", "-10 cm < zVertex < -5 cm", 100, etaMin, etaMax);
  std::unique_ptr<TH1F> MFTTracksEta5_10neg = std::make_unique<TH1F> ("MFT Tracks 5 10 eta", "5 cm < zVertex < 10 cm", 100, etaMin, etaMax);


  std::unique_ptr<TH1F> MCTracksp5 = std::make_unique<TH1F> ("MC Tracks 5 p", "-5 cm < zVertex < 5 cm", 100, 0, pMax);
  std::unique_ptr<TH1F> MCTracksp5_10pos = std::make_unique<TH1F> ("MC Tracks -5 -10 p", "-10 cm < zVertex < -5 cm", 100, 0, pMax);
  std::unique_ptr<TH1F> MCTracksp5_10neg = std::make_unique<TH1F> ("MC Tracks 5 10 p", "5 cm < zVertex < 10 cm", 100, 0, pMax);

  std::unique_ptr<TH1F> MFTTracksp5 = std::make_unique<TH1F> ("MFT Tracks 5 p", "-5 cm < zVertex < 5 cm", 100, 0, pMax);
  std::unique_ptr<TH1F> MFTTracksp5_10pos = std::make_unique<TH1F> ("MFT Tracks -5 -10 p", "-10 cm < zVertex < -5 cm", 100, 0, pMax);
  std::unique_ptr<TH1F> MFTTracksp5_10neg = std::make_unique<TH1F> ("MFT Tracks 5 10 p", "5 cm < zVertex < 10 cm", 100, 0, pMax);

  std::unique_ptr<TH1I> Trackablility = std::make_unique<TH1I> ("Trackablility", "In how many disks the tracks has hits", 6, 0, 6);

  //Histos for Missed (missed tracks that could be tracked)
  std::unique_ptr<TH1F> MissedlepT = std::make_unique<TH1F> ("Missed Tracks pT", "Missed Tracks pT", 100, 0, pMax);
  std::unique_ptr<TH1F> Missedp = std::make_unique<TH1F> ("Missed Tracks p", "Missed Tracks p", 100, 0, pMax);
  std::unique_ptr<TH1F> MissedRap = std::make_unique<TH1F> ("Missed Tracks eta", "Missed Rapidity", 100, etaMin, etaMax);




  //Histos for Trackables
  std::unique_ptr<TH1F> TrackablepT = std::make_unique<TH1F> ("Trackables Tracks pT", "Trackables Tracks pT", 100, 0, pMax);
  std::unique_ptr<TH1F> Trackablep = std::make_unique<TH1F> ("Trackables Tracks p", "Trackables Tracks p", 100, 0, pMax);
  std::unique_ptr<TH1F> TrackableRap = std::make_unique<TH1F> ("Trackables Tracks eta", "Trackables Rapidity", 100, etaMin, etaMax);


  //2D Histos
  std::unique_ptr<TH2F> MFTTrackedEtaZ = std::make_unique<TH2F> ("MFT_Tracked_eta_z", "Reconstructed Tracks: Rapidity vs zVertex", 31, -15, 16, 25, etaMin, etaMax);
  MFTTrackedEtaZ->GetXaxis()->SetTitle("Vertex PosZ [cm]");
  std::unique_ptr<TH2F> MFTAccepEtaZ = std::make_unique<TH2F> ("MFT_Acceptance_eta_z", "MFT Acceptance (Trackables): Rapidity vs zVertex", 31, -15, 16, 25, etaMin, etaMax);
  MFTAccepEtaZ->GetXaxis()->SetTitle("Vertex PosZ [cm]");
  std::unique_ptr<TH2F> MCTracksEtaZ = std::make_unique<TH2F> ("MCTracks_eta_z", "MC Tracks: Rapidity vs zVertex", 31, -15, 16, 25, etaMin, etaMax);
  MCTracksEtaZ->GetXaxis()->SetTitle("Vertex PosZ [cm]");



  TFile *simFileIn = new TFile(SimFile);
  TFile *trkFileIn = new TFile(trkFile);
  TFile outFile("MFT_eff_z.root","RECREATE");


  TTree *o2SimTree = (TTree*) simFileIn -> Get("o2sim");
  TTree *mftTrackTree = (TTree*) trkFileIn -> Get("o2sim");

  Int_t numberOfEvents = o2SimTree -> GetEntries();
  std::cout << "numberOfEvents = " << numberOfEvents << std::endl;

  Int_t nCleanTracksLTF = 0, nCleanTracksCA = 0, nBadTracksLTF =0, nBadTracksCA = 0;

  vector<Hit>* hit = nullptr;
  o2SimTree -> SetBranchAddress("MFTHit",&hit);
  vector<MCTrackT<float>>* mcTr = nullptr;
  o2SimTree -> SetBranchAddress("MCTrack",&mcTr);
  std::vector<o2::mft::TrackCA> trackCAVec, *trackCAVecP = &trackCAVec;
  mftTrackTree->SetBranchAddress("MFTTrackCA", &trackCAVecP);
  std::vector<o2::mft::TrackLTF> trackLTFVec, *trackLTFVecP = &trackLTFVec;
  mftTrackTree->SetBranchAddress("MFTTrackLTF", &trackLTFVecP);

  mftTrackTree->GetEntry(0);
  o2SimTree -> GetEntry(0);
  Int_t numberOfTracksPerEvent = mcTr->size(); // Number of tracks in first MCEvent (assumed to be the same for all events)
  vector<eventFoundTracks> allFoundTracksLTF(numberOfEvents+10), allFoundTracksCA(numberOfEvents+10); // True for reconstructed tracks - one vector of bool per event
  for (auto event = 0 ; event < numberOfEvents ; event++) { // Resize vector to accomodate found status of all tracks in all events
    //std::cout << "Resizing allFoundTracks for event " << event <<  " with ntracks = " << numberOfTracksPerEvent << std::endl;
    allFoundTracksLTF[event].resize(numberOfTracksPerEvent+10,false);
    allFoundTracksCA[event].resize(numberOfTracksPerEvent+10,false);
  }

  // 1. Loop over all reconstructed tracks to identify clean and mixed/noise tracks
  //   1.1 Clean tracks have at least 80% of its clusters from the same track
  //   1.2 If track is not clear it is a mixed/noise track
  // 2. Loop over all MC events to
  //   2.1 Identify trackable tracks (clusters in at least 4 disks)
  //   2.2 Identify successfully reconstructed tracks
  //   2.3 Fill Histograms


  // TracksLTF - Discover reconstructed tracks
  for (const auto &trackLTF : trackLTFVec) {
  auto thisTrackMCCompLabels = trackLTF.getMCCompLabels();
  std::map<Int_t, Int_t> trkIDs;
  for (auto iLabel = 0; iLabel < trackLTF.getNPoints(); iLabel++) { // Count trackIDs of this track
    auto id = thisTrackMCCompLabels[iLabel].getTrackID();
    trkIDs[id]=trkIDs[id]+1;
    //std::cout << thisTrackMCCompLabels[iLabel].getTrackID() << " ";
  }
  //std::cout << std::endl;

  Int_t thisTrkID=-1, thisTEventIDLabel=-1;
  for (auto iLabel = 0; iLabel < trackLTF.getNPoints(); iLabel++) { //Decide if track was successfully tracked
    auto id = thisTrackMCCompLabels[iLabel].getTrackID();
    //std::cout << "ID ; count ; NPoints ; count/NPoints : "<< id << " " << trkIDs[id] << " " << trackLTF.getNPoints()  << " " << 1.0*trkIDs[id]/trackLTF.getNPoints() << std::endl;
    if(1.0*trkIDs[id]/trackLTF.getNPoints() >= 0.8) {
      thisTrkID = id;
      thisTEventIDLabel = iLabel;
    }
  }
  auto eventID =  thisTrackMCCompLabels[thisTEventIDLabel].getEventID();
  //std::cout << "This TrackLTF ID = " << thisTrkID << " in eventID = " << eventID << std::endl;

  if( (thisTrkID > 0 & thisTrkID != 0x7FFFFFF) & (eventID <= numberOfEvents) ) {
   allFoundTracksLTF[eventID][thisTrkID]=true;
   nCleanTracksLTF++;
   }
   else {
    //std::cout << "Noise or Mixed TrackLTF!" << std::endl;
    nBadTracksLTF++;
  }

}


// TracksCA - Discover reconstructed tracks
for (const auto &trackCA : trackCAVec) {
auto thisTrackMCCompLabels = trackCA.getMCCompLabels();
std::map<Int_t, Int_t> trkIDs;
for (auto iLabel = 0; iLabel < trackCA.getNPoints(); iLabel++) { // Count trackIDs of this track
  auto id = thisTrackMCCompLabels[iLabel].getTrackID();
  trkIDs[id]=trkIDs[id]+1;
  //std::cout << thisTrackMCCompLabels[iLabel].getTrackID() << " ";
}
//std::cout << std::endl;

Int_t thisTrkID=-1, thisTEventIDLabel=-1;
for (auto iLabel = 0; iLabel < trackCA.getNPoints(); iLabel++) { //Decide if track was successfully tracked
  auto id = thisTrackMCCompLabels[iLabel].getTrackID();
  //std::cout << "ID ; count ; NPoints ; count/NPoints : "<< id << " " << trkIDs[id] << " " << trackCA.getNPoints()  << " " << 1.0*trkIDs[id]/trackCA.getNPoints() << std::endl;
  if(1.0*trkIDs[id]/trackCA.getNPoints() >= 0.8) {
    thisTrkID = id;
    thisTEventIDLabel = iLabel;
  }
}

auto eventID =  thisTrackMCCompLabels[thisTEventIDLabel].getEventID();
//std::cout << "This TrackCA ID = " << thisTrkID << " in eventID = " << eventID << std::endl;
if( (thisTrkID > 0 & thisTrkID != 0x7FFFFFF) & (eventID <= numberOfEvents) ) {
  allFoundTracksCA[eventID][thisTrkID]=true;
  nCleanTracksCA++;
}
else {
  //std::cout << "Noise or Mixed TrackCA!" << std::endl;
  nBadTracksCA++;
}

}



  for (Int_t event=0; event<numberOfEvents ; event++) { // Loop over events in o2sim
    //std::cout << "Loop over events in o2sim. Event = " << event << std::endl;
    o2SimTree -> GetEntry(event);
    Int_t nbH = hit->size(); // Number of hits in this event
    //std::cout << "Event " << event << " has " << numberOfTracksPerEvent << " tracks and " << nbH << " hits\n";

    std::vector<trackHasHitsinDisks> trackLayersHits(numberOfTracksPerEvent+10,{0,0,0,0,0}); //

    //std::cout << "Loop over " << nbH << " hits to discover trackable tracks in event " <<  event << std::endl;
    for (Int_t n_hit=0 ; n_hit < nbH; n_hit++) { // Loop over hits to discover trackable tracks
      Hit* hitp = &(*hit).at(n_hit);
      Int_t trID = hitp->GetTrackID(); // ID of the tracks having given the hit
      //std::cout << "n_hit = " << n_hit << " ** trID = " << trID << std::endl;

      Float_t z = hitp->GetZ(); // Z position of the hit = discover disk.
      for(auto disk: {0,1,2,3,4}) if( z < MFTLayerZ[disk*2] + .3  & z > MFTLayerZ[disk*2+1] -.3 ) trackLayersHits[trID][disk] = true;
      }

    for (Int_t trID=0 ; trID < numberOfTracksPerEvent; trID++) { // Loop on tracks
      //std::cout << "Loop on tracks to build histos. Track " << trID << " at event " << event << " -> " ;

      //fill MC histograms
      MCTrackT<float>* thisTrack =  &(*mcTr)[trID];
      auto z = thisTrack->GetStartVertexCoordinatesZ();
      auto eta = thisTrack->GetRapidity();
      auto p = thisTrack->GetP();
      MCTrackspT->Fill(thisTrack->GetPt());
      MCTracksp->Fill(p);
      MCTrackRap->Fill(eta);
      MCTracksEtaZ->Fill(z,eta);
      if( (z>-5) & (z<5) ) {
        MCTracksEta5->Fill(eta);
        MCTracksp5->Fill(p);
      }
      if( (z>5) & (z<10) ) {
        MCTracksEta5_10pos->Fill(eta);
        MCTracksp5_10pos->Fill(p);
      }
      if( (z>-10) & (z<-5) ) {
        MCTracksEta5_10neg->Fill(eta);
        MCTracksp5_10neg->Fill(p);
      }





      // Count disks "touched" by the track
      int nDisksHasHits = 0;
      for(auto disk: {0,1,2,3,4}) nDisksHasHits+= int(trackLayersHits[trID][disk]);
      Trackablility->Fill(nDisksHasHits);
      //std::cout << "nDisksHasHits = " << nDisksHasHits; // << std::endl;

      if(nDisksHasHits>=4) {   //Track is trackable if has left hits on at least 4 disks

        MFTAccepEtaZ->Fill(z,eta);


        bool wasFound = allFoundTracksLTF[event][trID] | allFoundTracksCA[event][trID];

        if(wasFound) {
          MFTTrackspT->Fill(thisTrack->GetPt());
          MFTTracksp->Fill(p);
          MFTTrackRap->Fill(eta);


          if( (z>-5) & (z<5) ) {
            MFTTracksEta5->Fill(eta);
            MFTTracksp5->Fill(p);
          }
          if( (z>5) & (z<10) ) {
            MFTTracksEta5_10pos->Fill(eta);
            MFTTracksp5_10pos->Fill(p);
          }
          if( (z>-10) & (z<-5) ) {
            MFTTracksEta5_10neg->Fill(eta);
            MFTTracksp5_10neg->Fill(p);
          }


          if(allFoundTracksLTF[event][trID]) {
            LTFTrackspT->Fill(thisTrack->GetPt());
            LTFTracksp->Fill(thisTrack->GetP());
            LTFTrackRap->Fill(thisTrack->GetRapidity());
            MFTTrackedEtaZ->Fill(thisTrack->GetStartVertexCoordinatesZ(),thisTrack->GetRapidity());


          }
          if(allFoundTracksCA[event][trID]) {
            CATrackspT->Fill(thisTrack->GetPt());
            CATracksp->Fill(thisTrack->GetP());
            CATrackRap->Fill(thisTrack->GetRapidity());
          }
        }
      } else {  // Fill histograms for Missed Tracks
        MissedlepT->Fill(thisTrack->GetPt());
        Missedp->Fill(thisTrack->GetP());
        MissedRap->Fill(thisTrack->GetRapidity());
      }
    //std::cout << " Finished Track " << trID << std::endl;

    } // end Loop on tracks
    //std::cout << "Finished event " << event << std::endl;
  } // end loop over events

//std::cout << "Building efficiencies histos..." << std::endl;
TH1F MFTEfficiencypT = (*MFTTrackspT)/ (*MCTrackspT);
TH1F MFTTEfficiencyp = (*MFTTracksp) / (*MCTracksp);
TH1F MFTEfficiencyRap = (*MFTTrackRap) / (*MCTrackRap);
TH2F MFTTrackerEfficiency = (*MFTTrackedEtaZ) / (*MFTAccepEtaZ);
TH2F MFTEfficiency2D = (*MFTTrackedEtaZ) / (*MCTracksEtaZ);

TH1F MFTEffsEta5 = (*MFTTracksEta5)/(*MCTracksEta5);
TH1F MFTEffEta5_10pos = (*MFTTracksEta5_10pos)/(*MCTracksEta5_10pos);
TH1F MFTEffEta5_10neg = (*MFTTracksEta5_10neg)/(*MCTracksEta5_10neg);


TH1F MFTEffsp5 = (*MFTTracksp5)/(*MCTracksp5);
TH1F MFTEffp5_10pos = (*MFTTracksp5_10pos)/(*MCTracksp5_10pos);
TH1F MFTEffp5_10neg = (*MFTTracksp5_10neg)/(*MCTracksp5_10neg);






MFTEfficiencypT.SetNameTitle("MFT Efficiency pT", "MFT Efficiency pT");
MFTTEfficiencyp.SetNameTitle("MFT Efficiency p", "MFT Efficiency p");
MFTEfficiencyRap.SetNameTitle("MFT Efficiency eta", "MFT Efficiency Rapidity");
MFTTrackerEfficiency.SetNameTitle("MFT Tracker Efficiency", "MFT Tracker Efficiency");
MFTTrackerEfficiency.GetXaxis()->SetTitle("Vertex PosZ [cm]");

MFTEfficiency2D.SetNameTitle("MFT Efficiency", "MFT Efficiency");
MFTEfficiency2D.GetXaxis()->SetTitle("Vertex PosZ [cm]");


MFTEffsEta5.SetNameTitle("MFT Eta Efficiency5_5", "-5 cm < z < 5 cm");
MFTEffsEta5.GetXaxis()->SetTitle("Rapidity");
MFTEffEta5_10pos.SetNameTitle("MFT Eta Efficiency_5_10pos", "5 cm < z < 10 cm");
MFTEffEta5_10pos.GetXaxis()->SetTitle("Rapidity");
MFTEffEta5_10neg.SetNameTitle("MFT Eta Efficiency_10_5neg", "-10 cm < z < -5 cm");
MFTEffEta5_10neg.GetXaxis()->SetTitle("Rapidity");

MFTEffsp5.SetNameTitle("MFT P Efficiency5_5", "-5 cm < z < 5 cm");
MFTEffsp5.GetXaxis()->SetTitle("P (GeV)");
MFTEffp5_10pos.SetNameTitle("MFT P Efficiency_5_10pos", "5 cm < z < 10 cm");
MFTEffp5_10pos.GetXaxis()->SetTitle("P (GeV)");
MFTEffp5_10neg.SetNameTitle("MFT P Efficiency_10_5neg", "-10 cm < z < -5 cm");
MFTEffp5_10neg.GetXaxis()->SetTitle("P (GeV)");

// Write histograms to file
//std::cout << "Writting histograms to file..." << std::endl;

MCTrackspT->Write();
MCTracksp->Write();
MCTrackRap->Write();

MFTTrackspT->Write();
MFTTracksp->Write();
MFTTrackRap->Write();

LTFTrackspT->Write();
LTFTracksp->Write();
LTFTrackRap->Write();

CATrackspT->Write();
CATracksp->Write();
CATrackRap->Write();

MFTEfficiencypT.Write();
MFTTEfficiencyp.Write();
MFTEfficiencyRap.Write();

MCTracksEta5->Write();
MCTracksEta5_10pos->Write();
MCTracksEta5_10neg->Write();

MFTTracksEta5->Write();
MFTTracksEta5_10pos->Write();
MFTTracksEta5_10neg->Write();

MCTracksp5->Write();
MCTracksp5_10pos->Write();
MCTracksp5_10neg->Write();


MFTTracksp5->Write();
MFTTracksp5_10pos->Write();
MFTTracksp5_10neg->Write();

MFTEffsEta5.Write();
MFTEffEta5_10pos.Write();
MFTEffEta5_10neg.Write();

MFTEffsp5.Write();
MFTEffp5_10pos.Write();
MFTEffp5_10neg.Write();


MissedlepT->Write();
Missedp->Write();
MissedRap->Write();

Trackablility->Write();


MCTracksEtaZ->SetOption("CONT4");
MCTracksEtaZ->Write();

MFTAccepEtaZ->SetOption("CONT4");
MFTAccepEtaZ->Write();

MFTTrackedEtaZ->SetOption("CONT4");
MFTTrackedEtaZ->Write();


MFTTrackerEfficiency.SetOption("CONT4");
MFTTrackerEfficiency.Write();

MFTEfficiency2D.SetOption("CONT4");
MFTEfficiency2D.Write();

outFile.Close();

Int_t totalRecoMFTTracks = nCleanTracksLTF + nCleanTracksCA + nBadTracksLTF + nBadTracksCA;
std::cout << "Total Reconstructed MFT Tracks = " << totalRecoMFTTracks << std::endl;
std::cout << "nCleanTracksLTF = " << nCleanTracksLTF << std::endl;
std::cout << "nCleanTracksCA = " << nCleanTracksCA << std::endl;
std::cout << "nBadTracksLTF = " << nBadTracksLTF  << " (" << 100.f*nBadTracksLTF/(nCleanTracksLTF+nBadTracksLTF) << " %)" << std::endl;
std::cout << "nBadTracksCA = " << nBadTracksCA << " (" << 100.f*nBadTracksCA/(nCleanTracksCA+nBadTracksCA) << " %)" << std::endl;
}
