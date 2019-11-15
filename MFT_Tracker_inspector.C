#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

constexpr Double_t LayerZ[] = {-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76.1, -77.5};


void MFT_Tracker_inspector(const Char_t *SimFile = "o2sim.root", const Char_t *trkFile = "mfttracks.root") {

  using o2::itsmft::Hit;
  using o2::MCTrackT;
  using o2::mft::TrackCA;
  using o2::mft::TrackLTF;
  using eventFoundTracks = std::vector<bool>;

  using trackHasHitsinDisks = std::array<bool,5>; // Disks with hits from a MFT track

  std::unique_ptr<TH1F> MCTrackspT = std::make_unique<TH1F> ("MC Tracks pT", "MC Tracks pT", 100, 0, 5);
  std::unique_ptr<TH1F> MCTracksp = std::make_unique<TH1F> ("MC Tracks p", "MC Tracks p", 100, 0, 5);
  std::unique_ptr<TH1F> MCTrackRap = std::make_unique<TH1F> ("MC Tracks eta", "MC Tracks Rapidity", 100, -4.0, -2);

  std::unique_ptr<TH1F> MFTTrackspT = std::make_unique<TH1F> ("MFT Tracks pT", "MFT Tracks pT", 100, 0, 5);
  std::unique_ptr<TH1F> MFTTracksp = std::make_unique<TH1F> ("MFT Tracks p", "MFT Tracks p", 100, 0, 5);
  std::unique_ptr<TH1F> MFTTrackRap = std::make_unique<TH1F> ("MFT Tracks eta", "MFT Tracks Rapidity", 100, -4.0, -2);

  //std::unique_ptr<TH1F> MFTEfficiencypT = std::make_unique<TH1F> ("MFT Efficiency pT", "MFT Efficiency pT", 100, 0, 5);
  //std::unique_ptr<TH1F> MFTTEfficiencyp = std::make_unique<TH1F> ("MFT Efficiency p", "MFT Efficiency p", 100, 0, 5);
  //std::unique_ptr<TH1F> MFTEfficiencyRap = std::make_unique<TH1F> ("MFT Efficiency eta", "MFT Efficiency Rapidity", 100, -4.0, -2);

  std::unique_ptr<TH1F> LTFTrackspT = std::make_unique<TH1F> ("LTF Tracks pT", "LTF Tracks pT", 100, 0, 5);
  std::unique_ptr<TH1F> LTFTracksp = std::make_unique<TH1F> ("LTF Tracks p", "LTF Tracks p", 100, 0, 5);
  std::unique_ptr<TH1F> LTFTrackRap = std::make_unique<TH1F> ("LTF Tracks eta", "LTF Tracks Rapidity", 100, -4.0, -2);

  std::unique_ptr<TH1F> CATrackspT = std::make_unique<TH1F> ("CA Tracks pT", "CA Tracks pT", 100, 0, 5);
  std::unique_ptr<TH1F> CATracksp = std::make_unique<TH1F> ("CA Tracks p", "CA Tracks p", 100, 0, 5);
  std::unique_ptr<TH1F> CATrackRap = std::make_unique<TH1F> ("CA Tracks eta", "LTF Tracks Rapidity", 100, -4.0, -2);

  std::unique_ptr<TH1I> Trackablility = std::make_unique<TH1I> ("Trackablility", "In how many disks the tracks has hits", 6, 0, 6);

  //Histos for untrackables
  std::unique_ptr<TH1F> UntrackablepT = std::make_unique<TH1F> ("Untrackables Tracks pT", "Untrackables Tracks pT", 100, 0, 5);
  std::unique_ptr<TH1F> Untrackablep = std::make_unique<TH1F> ("Untrackables Tracks p", "Untrackables Tracks p", 100, 0, 5);
  std::unique_ptr<TH1F> UntrackableRap = std::make_unique<TH1F> ("Untrackables Tracks eta", "Untrackables Rapidity", 100, -4.0, -2);


  TFile *simFileIn = new TFile(SimFile);
  TFile *trkFileIn = new TFile(trkFile);
  TFile outFile("MFT_eff_z.root","RECREATE");


  TTree *o2SimTree = (TTree*) simFileIn -> Get("o2sim");
  TTree *mftTrackTree = (TTree*) trkFileIn -> Get("o2sim");

  Int_t numberOfEvents = o2SimTree -> GetEntries();
  std::cout << "numberOfEvents = " << numberOfEvents << std::endl;



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
    allFoundTracksLTF[event].resize(numberOfTracksPerEvent,false);
    allFoundTracksCA[event].resize(numberOfTracksPerEvent,false);
  }

  // TracksCA
  for (const auto &trackCA : trackCAVec) {
  auto thisTrackMCCompLabels = trackCA.getMCCompLabels();
  auto firstTrackID = thisTrackMCCompLabels[0].getTrackID();
  auto eventID =  thisTrackMCCompLabels[0].getEventID();
  allFoundTracksCA[eventID][firstTrackID]=true;
  //std::cout << "Found TrackCA " << firstTrackID << " in event " << eventID << std::endl;
  //for (auto iLabel = 0; iLabel < trackCA.getNPoints(); iLabel++) {
  //  std::cout << thisTrackMCCompLabels[iLabel].getTrackID() << " ";
  //  Count how many times the clusters from initial MC tracks are present
  //  }
  //  std::cout << std::endl;
  }

  // TracksLTF
  for (const auto &trackLTF : trackLTFVec) {
  auto thisTrackMCCompLabels = trackLTF.getMCCompLabels();
  auto firstTrackID = thisTrackMCCompLabels[0].getTrackID();
  auto eventID =  thisTrackMCCompLabels[0].getEventID();
  allFoundTracksLTF[eventID][firstTrackID]=true;
  //std::cout << "Found TrackLTF " << firstTrackID << " in event " << eventID << std::endl;
  //for (auto iLabel = 0; iLabel < trackLTF.getNPoints(); iLabel++) {
  //  std::cout << thisTrackMCCompLabels[iLabel].getTrackID() << " ";
  //  Count how many times the clusters from initial MC tracks are present
  //  }
  //  std::cout << std::endl;
  }




  for (Int_t event=0; event<numberOfEvents ; event++) { // Loop over events in o2sim
    o2SimTree -> GetEntry(event);
    //Int_t numberOfTracksPerEvent = mcTr->size(); // Number of tracks in this event
    Int_t nbH = hit->size(); // Number of hits in this event
    //std:cout << "Event " << event << " has " << numberOfTracksPerEvent << " tracks and " << nbH << " hits\n";


    std::vector<trackHasHitsinDisks> trackLayersHits(numberOfTracksPerEvent,{0,0,0,0,0}); //

    for (Int_t n_hit=0 ; n_hit < nbH; n_hit++) { // Loop over hits to discover trackable tracks
      Hit* hitp = &(*hit)[n_hit];
      Int_t trID = hitp->GetTrackID(); // ID of the tracks having given the hit
      Float_t z = hitp->GetZ(); // Z position of the hit = discover disk.
      for(auto disk: {0,1,2,3,4}) if( z < LayerZ[disk*2] + .3  & z > LayerZ[disk*2+1] -.3 ) trackLayersHits[trID][disk] = true;
       }

    for (Int_t trID=0 ; trID < numberOfTracksPerEvent; trID++) { // Loop on tracks

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
        bool reconstructed = allFoundTracksLTF[event][trID] | allFoundTracksCA[event][trID];

        if(reconstructed) {
          MFTTrackspT->Fill(thisTrack->GetPt());
          MFTTracksp->Fill(thisTrack->GetP());
          MFTTrackRap->Fill(thisTrack->GetRapidity());
        }
      } else {  // Fill historgrams for untrackables
        UntrackablepT->Fill(thisTrack->GetPt());
        Untrackablep->Fill(thisTrack->GetP());
        UntrackableRap->Fill(thisTrack->GetRapidity());

      }

    } // end Loop on tracks
  } // end loop over events

TH1F MFTEfficiencypT = (*MFTTrackspT)/ (*MCTrackspT);
TH1F MFTTEfficiencyp = (*MFTTracksp) / (*MCTracksp);
TH1F MFTEfficiencyRap = (*MFTTrackRap) / (*MCTrackRap);

MFTEfficiencypT.SetNameTitle("MFT Efficiency pT", "MFT Efficiency pT");
MFTTEfficiencyp.SetNameTitle("MFT Efficiency p", "MFT Efficiency p");
MFTEfficiencyRap.SetNameTitle("MFT Efficiency eta", "MFT Efficiency Rapidity");

// Write Histograms to file
MCTrackspT->Write();
MCTracksp->Write();
MCTrackRap->Write();

MFTTrackspT->Write();
MFTTracksp->Write();
MFTTrackRap->Write();

MFTEfficiencypT.Write();
MFTTEfficiencyp.Write();
MFTEfficiencyRap.Write();

UntrackablepT->Write();
Untrackablep->Write();
UntrackableRap->Write();

Trackablility->Write();

outFile.Close();

}
