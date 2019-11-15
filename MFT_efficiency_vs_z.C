#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

constexpr Double_t LayerZ[] = {-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76.1, -77.5};


void MFT_efficiency_vs_z(const Char_t *nameFile = "o2sim.root") {

  using o2::itsmft::Hit;
  using o2::MCTrackT;
  using trackHasHitsinDisks = std::array<bool,5>; // Disks with hits from a MFT track

  std::unique_ptr<TH1F> MCTrackspT = std::make_unique<TH1F> ("MC Tracks pT", "MC Tracks pT", 100, 0, 5);
  std::unique_ptr<TH1F> MCTracksp = std::make_unique<TH1F> ("MC Tracks p", "MC Tracks p", 100, 0, 5);
  std::unique_ptr<TH1F> MCTrackRap = std::make_unique<TH1F> ("MC Tracks eta", "MC Tracks Rapidity", 100, -4.0, -2);

  std::unique_ptr<TH1F> MFTTrackspT = std::make_unique<TH1F> ("MFT Tracks pT", "MFT Tracks pT", 100, 0, 5);
  std::unique_ptr<TH1F> MFTTracksp = std::make_unique<TH1F> ("MFT Tracks p", "MFT Tracks p", 100, 0, 5);
  std::unique_ptr<TH1F> MFTTrackRap = std::make_unique<TH1F> ("MFT Tracks eta", "MFT Tracks Rapidity", 100, -4.0, -2);

  std::unique_ptr<TH1I> Trackablility = std::make_unique<TH1I> ("Trackablility", "In how many disks the tracks has hits", 6, 0, 6);

  //Histos for untrackables
  std::unique_ptr<TH1F> UntrackablepT = std::make_unique<TH1F> ("Untrackables Tracks pT", "Untrackables Tracks pT", 100, 0, 5);
  std::unique_ptr<TH1F> Untrackablep = std::make_unique<TH1F> ("Untrackables Tracks p", "Untrackables Tracks p", 100, 0, 5);
  std::unique_ptr<TH1F> UntrackableRap = std::make_unique<TH1F> ("Untrackables Tracks eta", "Untrackables Rapidity", 100, -4.0, -2);





  TFile *fileIn = new TFile(nameFile);
  TTree *tr = (TTree*) fileIn -> Get("o2sim");

  TFile outFile("MFT_eff_z.root","RECREATE");
  vector<Hit>* hit = nullptr;
  tr -> SetBranchAddress("MFTHit",&hit);
  vector<MCTrackT<float>>* mcTr = nullptr;
  tr -> SetBranchAddress("MCTrack",&mcTr);

  Int_t nbE = tr -> GetEntries();

  std::cout << "nbE = " << nbE << std::endl;
  for (Int_t event=0; event<nbE ; event++) { // Loop over events in o2sim
    tr -> GetEntry(event);
    Int_t nbTr = mcTr->size(); // Number of tracks in this event
    Int_t nbH = hit->size(); // Number of hits in this event
    std:cout << "Event " << event << " has " << nbTr << " tracks and " << nbH << " hits\n";
    //std::cout << "nbTr = " << nbTr << std::endl;
    //s/td::cout << "nbH = " << nbH << std::endl;

    std::vector<trackHasHitsinDisks> trackLayersHits(nbTr,{0,0,0,0,0}); //

    for (Int_t n_hit=0 ; n_hit < nbH; n_hit++) { // Loop over hits to discover trackable tracks
      Hit* hitp = &(*hit)[n_hit];
      Int_t trID = hitp->GetTrackID(); // ID of the tracks having given the hit
      Float_t z = hitp->GetZ(); // Z position of the hit = discover disk.
      for(auto disk: {0,1,2,3,4}) if( z < LayerZ[disk*2] + .3  & z > LayerZ[disk*2+1] -.3 ) trackLayersHits[trID][disk] = true;
       }

    for (Int_t trID=0 ; trID < nbTr; trID++) { // Loop on tracks

      //fill MC histograms
      MCTrackT<float>* thisTrack =  &(*mcTr)[trID];
      MCTrackspT->Fill(thisTrack->GetPt());
      MCTracksp->Fill(thisTrack->GetP());
      MCTrackRap->Fill(thisTrack->GetRapidity());

      // Count disks "touched" by the track
      int nDisksHasHits = 0;
      for(auto disk: {0,1,2,3,4}) nDisksHasHits+= int(trackLayersHits[trID][disk]);
      Trackablility->Fill(nDisksHasHits);
      std::cout << "nDisksHasHits = " << nDisksHasHits << std::endl;

      if(nDisksHasHits>=4) {   //Track is trackable if has left hits on at least 4 disks
        bool reconstructed = true; // TODO: discover if track has been found

        // Load mfttracks.root
        // Loop on Events
        // Loop on TracksLFT and TracksCA and collect trackIDd -> Build an vector/array of tracked Tracks IDs



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

// Write Histograms to file
MCTrackspT->Write();
MCTracksp->Write();
MCTrackRap->Write();

MFTTrackspT->Write();
MFTTracksp->Write();
MFTTrackRap->Write();

UntrackablepT->Write();
Untrackablep->Write();
UntrackableRap->Write();

Trackablility->Write();

outFile.Close();

}
