#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

constexpr Double_t LayerZ[] = {-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76.1, -77.5};

constexpr Double_t pMax = 1.5;
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

  std::unique_ptr<TH1I> Trackablility = std::make_unique<TH1I> ("Trackablility", "In how many disks the tracks has hits", 6, 0, 6);

  //Histos for untrackables
  std::unique_ptr<TH1F> UntrackablepT = std::make_unique<TH1F> ("Untrackables Tracks pT", "Untrackables Tracks pT", 100, 0, pMax);
  std::unique_ptr<TH1F> Untrackablep = std::make_unique<TH1F> ("Untrackables Tracks p", "Untrackables Tracks p", 100, 0, pMax);
  std::unique_ptr<TH1F> UntrackableRap = std::make_unique<TH1F> ("Untrackables Tracks eta", "Untrackables Rapidity", 100, etaMin, etaMax);


  TFile *simFileIn = new TFile(SimFile);
  TFile *trkFileIn = new TFile(trkFile);
  TFile outFile("MFT_eff_z.root","RECREATE");


  TTree *o2SimTree = (TTree*) simFileIn -> Get("o2sim");
  TTree *mftTrackTree = (TTree*) trkFileIn -> Get("o2sim");

  Int_t numberOfEvents = o2SimTree -> GetEntries();
  std::cout << "numberOfEvents = " << numberOfEvents << std::endl;

  Int_t NInvalid_tracks = 0;
  std::vector<int> eventsWithInvalidIDs;

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
  vector<eventFoundTracks> allFoundTracksLTF(numberOfEvents), allFoundTracksCA(numberOfEvents); // True for reconstructed tracks - one vector of bool per event
  for (auto event = 0 ; event < numberOfEvents ; event++) { // Resize vector to accomodate found status of all tracks in all events
    //std::cout << "Resizing allFoundTracks for event " << event <<  " with ntracks = " << numberOfTracksPerEvent << std::endl;
    allFoundTracksLTF[event].resize(numberOfTracksPerEvent,false);
    allFoundTracksCA[event].resize(numberOfTracksPerEvent,false);

  }

  // 1. Loop over all reconstructed tracks to build found tracks identities
  //   1.1 Clean tracks have at least 80% of its clusters from the same tracks
  //   1.2 If track is not clear it is a mixed and noise track
  // 2. Loop over all MC events to
  //   2.1 Discover trackable tracks
  //   2.2 Discover which tracks have been successfully reconstructed
  //   2.3 Fill Histograms


  // TracksLTF
  for (const auto &trackLTF : trackLTFVec) {
  auto thisTrackMCCompLabels = trackLTF.getMCCompLabels();
  auto firstTrackID = thisTrackMCCompLabels[0].getTrackID();
  auto eventID =  thisTrackMCCompLabels[0].getEventID();
  if(!allFoundTracksLTF[eventID][firstTrackID]) allFoundTracksLTF[eventID][firstTrackID]=true; else std::cout << "TRACK Found Elsewere!!!\n";
  //std::cout << "Found TrackLTF " << firstTrackID << " in event " << eventID << std::endl;
  //for (auto iLabel = 0; iLabel < trackLTF.getNPoints(); iLabel++) {
  //  std::cout << thisTrackMCCompLabels[iLabel].getTrackID() << " ";
  //}
  //  std::cout << std::endl;
  }


  // TracksCA
  for (const auto &trackCA : trackCAVec) {
  auto thisTrackMCCompLabels = trackCA.getMCCompLabels();
  auto firstTrackID = thisTrackMCCompLabels[0].getTrackID();
  auto eventID =  thisTrackMCCompLabels[0].getEventID();
  //std::cout << "Found TrackCA " << firstTrackID << " in event " << eventID << std::endl;
  //for (auto iLabel = 0; iLabel < trackCA.getNPoints(); iLabel++) {
  //  std::cout << thisTrackMCCompLabels[iLabel].getTrackID() << " ";
  //}
    if (eventID < numberOfEvents & firstTrackID < numberOfTracksPerEvent ) {
      allFoundTracksCA.at(eventID).at(firstTrackID)=true;
    } else {
      std::cout << "TrackCA InvalidID: " << firstTrackID << " in event " << eventID << std::endl;
      NInvalid_tracks++;
      eventsWithInvalidIDs.push_back(eventID);
      for (auto iLabel = 0; iLabel < trackCA.getNPoints(); iLabel++) {
        std::cout << thisTrackMCCompLabels[iLabel].getTrackID() << " ";
      }
      std::cout << std::endl;

    }
    }

  for (Int_t event=0; event<numberOfEvents ; event++) { // Loop over events in o2sim
    //std::cout << "Loop over events in o2sim. Event = " << event << std::endl;
    o2SimTree -> GetEntry(event);
    Int_t nbH = hit->size(); // Number of hits in this event
    //std::cout << "Event " << event << " has " << numberOfTracksPerEvent << " tracks and " << nbH << " hits\n";
    bool invalid=false;
    for(auto Id: eventsWithInvalidIDs) if (Id==event) {
      invalid=true;
      std::cout << " ***** Event " << event << " had a track with an invalid ID. Skipping...  *****" << std::endl;
    }
    if(invalid) continue;
    std::vector<trackHasHitsinDisks> trackLayersHits(numberOfTracksPerEvent,{0,0,0,0,0}); //

    //std::cout << "Loop over " << nbH << " hits to discover trackable tracks in event " <<  event << std::endl;
    for (Int_t n_hit=0 ; n_hit < nbH; n_hit++) { // Loop over hits to discover trackable tracks
      Hit* hitp = &(*hit).at(n_hit);
      Int_t trID = hitp->GetTrackID(); // ID of the tracks having given the hit
      //std::cout << "n_hit = " << n_hit << " ** trID = " << trID << std::endl;

      Float_t z = hitp->GetZ(); // Z position of the hit = discover disk.
      for(auto disk: {0,1,2,3,4}) if( z < LayerZ[disk*2] + .3  & z > LayerZ[disk*2+1] -.3 ) trackLayersHits[trID][disk] = true;
      }

    for (Int_t trID=0 ; trID < numberOfTracksPerEvent; trID++) { // Loop on tracks
      //std::cout << "Loop on tracks to build histos. Track " << trID << " at event " << event << std::endl;

      //fill MC histograms
      MCTrackT<float>* thisTrack =  &(*mcTr)[trID];
      MCTrackspT->Fill(thisTrack->GetPt());
      MCTracksp->Fill(thisTrack->GetP());
      MCTrackRap->Fill(thisTrack->GetRapidity());

      // Count disks "touched" by the track
      int nDisksHasHits = 0;
      for(auto disk: {0,1,2,3,4}) nDisksHasHits+= int(trackLayersHits[trID][disk]);
      Trackablility->Fill(nDisksHasHits);
      //std::cout << "nDisksHasHits = " << nDisksHasHits << std::endl;

      if(nDisksHasHits>=4) {   //Track is trackable if has left hits on at least 4 disks
        bool wasFound = allFoundTracksLTF[event][trID] | allFoundTracksCA[event][trID];

        if(wasFound) {
          MFTTrackspT->Fill(thisTrack->GetPt());
          MFTTracksp->Fill(thisTrack->GetP());
          MFTTrackRap->Fill(thisTrack->GetRapidity());
          if(allFoundTracksLTF[event][trID]) {
            LTFTrackspT->Fill(thisTrack->GetPt());
            LTFTracksp->Fill(thisTrack->GetP());
            LTFTrackRap->Fill(thisTrack->GetRapidity());
          }
          if(allFoundTracksCA[event][trID]) {
            CATrackspT->Fill(thisTrack->GetPt());
            CATracksp->Fill(thisTrack->GetP());
            CATrackRap->Fill(thisTrack->GetRapidity());
          }
        }
      } else {  // Fill histograms for untrackables
        UntrackablepT->Fill(thisTrack->GetPt());
        Untrackablep->Fill(thisTrack->GetP());
        UntrackableRap->Fill(thisTrack->GetRapidity());

      }

    } // end Loop on tracks
    //std::cout << "Finished event" << event << std::endl;
  } // end loop over events

std::cout << "Building efficiencies histos..." << std::endl;
TH1F MFTEfficiencypT = (*MFTTrackspT)/ (*MCTrackspT);
TH1F MFTTEfficiencyp = (*MFTTracksp) / (*MCTracksp);
TH1F MFTEfficiencyRap = (*MFTTrackRap) / (*MCTrackRap);

MFTEfficiencypT.SetNameTitle("MFT Efficiency pT", "MFT Efficiency pT");
MFTTEfficiencyp.SetNameTitle("MFT Efficiency p", "MFT Efficiency p");
MFTEfficiencyRap.SetNameTitle("MFT Efficiency eta", "MFT Efficiency Rapidity");

// Write histograms to file
std::cout << "Writting histograms to file..." << std::endl;

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

UntrackablepT->Write();
Untrackablep->Write();
UntrackableRap->Write();

Trackablility->Write();

outFile.Close();

std::cout << "NInvalid_tracks = " << NInvalid_tracks << std::endl;
}
