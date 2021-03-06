#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

//constexpr Double_t MFTLayerZ[] = {-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76.1, -77.5};

constexpr Double_t pMax = 5;
constexpr Double_t pMin = 0.1;
constexpr Double_t etaMin = -4.0;
constexpr Double_t etaMax = -2.0;

void MultiplicityEstimator_noReco(const Char_t *SimFile = "o2sim.root") {

  using o2::itsmft::Hit;
  using o2::MCTrackT;
  o2::itsmft::ChipMappingMFT mftChipMapper;
  std::vector<Int_t> trackables_per_event;
  Int_t nMFTTrackables=0;

  using trackHasHitsinMFTDisks = std::array<bool,5>; // Disks with hits from a MFT track

  std::unique_ptr<TH1F> MCTrackspT = std::make_unique<TH1F> ("MC Tracks pT", "MC Tracks pT", 100, 0, pMax);
  MCTrackspT->GetXaxis()->SetTitle("Transverse p");

  std::unique_ptr<TH1F> MCTracksp = std::make_unique<TH1F> ("MC Tracks p", "MC Tracks p", 100, 0, pMax);
  MCTracksp->GetXaxis()->SetTitle("Total p");

  std::unique_ptr<TH1F> MCTracksEta = std::make_unique<TH1F> ("MC Tracks eta", "MC Tracks Pseudorapidity", 100, etaMin, etaMax);
  MCTracksEta->GetXaxis()->SetTitle("Pseudorapidity");

  std::unique_ptr<TH2F> hitVertex = std::make_unique<TH2F> ("HitVertexes", "Vertexes of Tracks with hits in the MFT", 1000, -200, 200, 1000, -200, 200);
  hitVertex->GetXaxis()->SetTitle("Z");
  hitVertex->GetYaxis()->SetTitle("X");

  std::unique_ptr<TH1I> Trackablility = std::make_unique<TH1I> ("Trackablility", "In how many disks the tracks has hits", 6, 0, 6);
  Trackablility->GetXaxis()->SetTitle("Number of disks");


  std::unique_ptr<TH1I> MultiplicityDistrib = std::make_unique<TH1I> ("Multiplicity", "MFT Trackables Distribution", 10000, 0, 10000);
  MultiplicityDistrib->GetXaxis()->SetTitle("Trackable Multiplicity");




  //Histos for Trackables
  std::unique_ptr<TH1F> TrackablepT = std::make_unique<TH1F> ("Trackables Tracks pT", "Trackables Tracks pT", 100, 0, pMax);
  TrackablepT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> Trackablep = std::make_unique<TH1F> ("Trackables Tracks p", "Trackables Tracks p", 100, 0, pMax);
  Trackablep->GetXaxis()->SetTitle("Total p");

  std::unique_ptr<TH1F> MFTTrackablesEta = std::make_unique<TH1F> ("Trackables Tracks eta", "Trackables Pseudorapidity", 100, etaMin, etaMax);
  MFTTrackablesEta->GetXaxis()->SetTitle("Pseudorapidity");





  TFile *simFileIn = new TFile(SimFile);
  //TFile *trkFileIn = new TFile(trkFile);
  TFile outFile("MultiplicityEstimator_noReco.root","RECREATE");


  TTree *o2SimTree = (TTree*) simFileIn -> Get("o2sim");
  //TTree *mftTrackTree = (TTree*) trkFileIn -> Get("o2sim");

  Int_t numberOfEvents = o2SimTree -> GetEntries();
  std::cout << "numberOfEvents = " << numberOfEvents << std::endl;

  //Int_t nCleanTracksLTF = 0, nCleanTracksCA = 0, nInvalidTracksLTF =0, nInvalidTracksCA = 0;

  vector<Hit>* mfthit = nullptr;
  o2SimTree -> SetBranchAddress("MFTHit",&mfthit);
  vector<MCTrackT<float>>* mcTr = nullptr;
  o2SimTree -> SetBranchAddress("MCTrack",&mcTr);


  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  o2SimTree -> SetBranchAddress("MCEventHeader.",&eventHeader);

//  mftTrackTree->GetEntry(0);
  o2SimTree -> GetEntry(0);
  //Int_t numberOfTracksPerEvent = mcTr->size(); // Number of tracks in first MCEvent (assumed to be the same for all events)
//  vector<eventFoundTracks> allFoundTracksLTF(numberOfEvents), allFoundTracksCA(numberOfEvents); // True for reconstructed tracks - one vector of bool per event


/*  for (auto event = 0 ; event < numberOfEvents ; event++) { // Resize vector to accomodate found status of all tracks in all events
    o2SimTree -> GetEntry(event);
    //o2::dataformats::MCEventStats& eventstats = eventHeader->getMCEventStats();
    auto numberOfTracksThisEvent = eventHeader->getMCEventStats().getNKeptTracks();
    //std::cout << "Resizing allFoundTracks for event " << event <<  " with ntracks = " << numberOfTracksThisEvent << std::endl;
    allFoundTracksLTF[event].resize(numberOfTracksThisEvent,false);
    allFoundTracksCA[event].resize(numberOfTracksThisEvent,false);
  }

  // Part 1: Quality of reconstructed MFT tracks
  //   - Loop over reconstructed tracks to identify clean and mixed/noise tracks
  //   - Clean tracks have at least 80% of its clusters from the same track
  //   - If track is not clear it is a mixed/noise track
  // Part 2: MC hits and tracks
  //   - Identify trackable tracks (clusters in at least 4 disks)
  //   - Identify successfully reconstructed tracks
  // Part 3: Calculate Efficiencies


  // Part 1: Reconstructed MFT Tracks
  //std::cout << "Starting Part 1: Reconstructed MFT Tracks!" << std::endl;

  // TracksLTF - Identify reconstructed tracks
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
  for (auto iLabel = 0; iLabel < trackLTF.getNPoints(); iLabel++) { // Decide if track was successfully tracked
    auto id = thisTrackMCCompLabels[iLabel].getTrackID();
    //std::cout << "ID ; count ; NPoints ; count/NPoints : "<< id << " " << trkIDs[id] << " " << trackLTF.getNPoints()  << " " << 1.0*trkIDs[id]/trackLTF.getNPoints() << std::endl;
    if(1.0*trkIDs[id]/trackLTF.getNPoints() >= 0.8) { // Must have at least 80% of its clusters from the same MC Track
      thisTrkID = id;
      thisTEventIDLabel = iLabel;
    }
  }
  auto eventID =  thisTrackMCCompLabels[thisTEventIDLabel].getEventID();
  //std::cout << "This TrackLTF ID = " << thisTrkID << " in eventID = " << eventID << std::endl;

  if( (thisTrkID > 0 & thisTrkID != 0x7FFFFFF) & (eventID <= numberOfEvents) ) { // If is good match and not noise ...
   allFoundTracksLTF[eventID][thisTrkID]=true;
   nCleanTracksLTF++;
   }
   else {
    //std::cout << "Noise or Mixed TrackLTF!" << std::endl;
    nInvalidTracksLTF++;
  }

  } // Loop on TracksLTF


  // TracksCA - Identify reconstructed tracks
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
  if( (thisTrkID > 0 & thisTrkID != 0x7FFFFFF) & (eventID <= numberOfEvents) ) { // If is good match and not noise ...
    allFoundTracksCA[eventID][thisTrkID]=true;
    nCleanTracksCA++;
  }
  else {
    //std::cout << "Noise or Mixed TrackCA!" << std::endl;
    nInvalidTracksCA++;
  }

  } // Loop on TracksCA
  */


  // Part 2: MC hits and tracks
  // std::cout << "Starting Part 2: MC hits and tracks!" << std::endl;
  for (Int_t event=0; event<numberOfEvents ; event++) { // Loop over events in o2sim
    // std::cout << "Loop over events in o2sim. Event = " << event << std::endl;
    o2SimTree -> GetEntry(event);
    MCTrackT<float>* thisTrack;// =  &(*mcTr)[trID];
    Int_t nMFTHits = mfthit->size(); // Number of mft hits in this event
    Int_t trackables_in_this_event=0;

    // std::cout << "Event " << event << " has " << eventHeader->getMCEventStats().getNKeptTracks() << " tracks and " << nMFTHits << " hits\n";

    std::vector<trackHasHitsinMFTDisks> mcTrackHasHitsInMFTDisks(eventHeader->getMCEventStats().getNKeptTracks(),{0,0,0,0,0}); //

    // std::cout << "Loop over " << nMFTHits << " mfthits to identify trackable MFT tracks in event " <<  event << std::endl;
    for (Int_t n_hit=0 ; n_hit < nMFTHits; n_hit++) { // Loop over mfthits to identify trackable tracks
      Hit* hitp = &(*mfthit).at(n_hit);
      Int_t trID = hitp->GetTrackID(); // ID of the tracks having given the hit
      //std::cout << "n_hit = " << n_hit << " ** trID = " << trID << std::endl;

      Float_t z = hitp->GetZ(); // Z position of the hit => Identify MFT disk
      Int_t chipID = hitp->GetDetectorID();
      thisTrack = &(*mcTr)[trID];
      auto vz = thisTrack->GetStartVertexCoordinatesZ();
      auto vx = thisTrack->GetStartVertexCoordinatesX();
      hitVertex->Fill(vz,vx);
      mcTrackHasHitsInMFTDisks[trID][mftChipMapper.chip2Layer(hitp->GetDetectorID())/2] = true;
      }

    for (Int_t trID=0 ; trID < eventHeader->getMCEventStats().getNKeptTracks(); trID++) { // Loop on MC tracks
      //std::cout << "Loop on tracks to build histos. Track " << trID << " at event " << event << " -> " ;

      //fill MC histograms
      thisTrack =  &(*mcTr)[trID];
      auto z = thisTrack->GetStartVertexCoordinatesZ();
      auto p = thisTrack->GetP();
      auto eta = atanh (thisTrack->GetStartVertexMomentumZ()/p); // eta;
      MCTrackspT->Fill(thisTrack->GetPt());
      MCTracksp->Fill(p);
      MCTracksEta->Fill(eta);


      // Count disks "touched" by the track
      int nMFTDisksHasHits = 0;
      for(auto disk: {0,1,2,3,4}) nMFTDisksHasHits+= int(mcTrackHasHitsInMFTDisks[trID][disk]);
      Trackablility->Fill(nMFTDisksHasHits);
      //std::cout << "nMFTDisksHasHits = " << nMFTDisksHasHits; // << std::endl;

      if(nMFTDisksHasHits>=4) {   //Track is trackable if has left hits on at least 4 disks

        MFTTrackablesEta->Fill(eta);
        trackables_in_this_event++;
}

    //std::cout << " Finished Track " << trID << std::endl;

    } // end Loop on tracks
    trackables_per_event.push_back(trackables_in_this_event);
    nMFTTrackables+=trackables_in_this_event;
    MultiplicityDistrib->Fill(trackables_in_this_event);

    //std::cout << "Finished event " << event << std::endl;
    } // end loop over events


// Write histograms to file
//std::cout << "Writting histograms to file..." << std::endl;

MCTrackspT->Write();
MCTracksp->Write();
MCTracksEta->Write();


Trackablility->Write();
MultiplicityDistrib->Write();

hitVertex->Write();


MFTTrackablesEta->Write();



outFile.Close();

//Int_t totalRecoMFTTracks = nCleanTracksLTF + nCleanTracksCA + nInvalidTracksLTF + nInvalidTracksCA;
std::cout << "Number of MFT trackables = " << nMFTTrackables << std::endl;
//new TBrowser;
}
