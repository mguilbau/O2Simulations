#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "DataFormatsITSMFT/ROFRecord.h"
#include "MFTTracking/TrackCA.h"
#include <FairLogger.h>
#include <TChain.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <vector>
#endif

void MFTTrackInspector(
    std::string outName =
        "MFTTrackInspector.root",           // name of the output binary file
    std::string inpName = "mfttracks.root", // name of the input MFT tracks
    std::string trackTreeName = "o2sim",    // name of the tracks tree
    std::string trackBranchName = "MFTTrackLTF", // name of the tracks branch
    std::string rofRecName =
        "MFTTracksROF" // name of the ROF records tree and its branch
) {
  TStopwatch swTot;
  swTot.Start();
  using ROFR = o2::itsmft::ROFRecord;
  using ROFRVEC = std::vector<o2::itsmft::ROFRecord>;

  ///-------> input
  TChain trackTree(trackTreeName.c_str());
  TChain rofTree(rofRecName.c_str());

  trackTree.AddFile(inpName.c_str());
  rofTree.AddFile(inpName.c_str());

  // TrackLTF record entries in the track tree
  std::vector<o2::mft::TrackLTF> trackLTFVec, *trackLTFVecP = &trackLTFVec;
  if (!trackTree.GetBranch(trackBranchName.c_str())) {
    LOG(FATAL) << "Failed to find the branch " << trackBranchName
               << " in the tree " << trackTreeName;
  }
  trackTree.SetBranchAddress(trackBranchName.c_str(), &trackLTFVecP);

  // ROF record entries in the track tree
  ROFRVEC rofRecVec, *rofRecVecP = &rofRecVec;
  if (!rofTree.GetBranch(rofRecName.c_str())) {
    LOG(FATAL) << "Failed to find the branch " << rofRecName << " in the tree "
               << rofRecName;
  }
  rofTree.SetBranchAddress(rofRecName.c_str(), &rofRecVecP);

  ///-------< input

  // Loop on
  // rofTree-------------------------------------------------------------------------------<<<<
  int lastTreeID = -1;
  long offs = 0, nEntProc = 0;
  LOG(INFO) << "rofTree.GetEntries(): " << rofTree.GetEntries() << " \n";
  for (int i = 0; i < rofTree.GetEntries(); i++) {
    rofTree.GetEntry(i);
    if (rofTree.GetTreeNumber() >
        lastTreeID) {       // this part is needed for chained input
      if (lastTreeID > 0) { // new chunk, increase the offset
        offs += trackTree.GetTree()->GetEntries();
      }
      lastTreeID = rofTree.GetTreeNumber();
    }

    for (const auto &rofRec : rofRecVec) {
      auto rofEntry = rofRec.getROFEntry();
      int nTrackROF = rofRec.getNROFEntries();
      LOG(INFO) << "Processing ROF:" << rofRec.getROFrame() << " with "
                << nTrackROF << " entries";
      if (!nTrackROF) {
        LOG(INFO) << "Frame is empty"; // ??
        continue;
      }
      if (rofEntry.getEvent() != trackTree.GetReadEntry() + offs || !nEntProc) {
        trackTree.GetEntry(rofEntry.getEvent() +
                           offs); // read tree entry containing needed ROF data
        nEntProc++;
      }
      int trackIndex =
          rofEntry.getIndex(); // needed ROF tracks start from this one
      int maxTrackIndex = trackIndex + nTrackROF;
      LOG(INFO) << "BV===== trackIndex " << trackIndex << " maxTrackIndex "
                << maxTrackIndex << "\n";
    }
  } // loop over multiple ROFvectors (in case of chaining)

  // Loop on
  // trackTree-------------------------------------------------------------------------------<<<<
  int lastTracktreeID = -1;
  long trackoffs = 0, nEntProcTrack = 0;
  LOG(INFO) << "trackTree.GetEntries(): " << trackTree.GetEntries();
  for (int i = 0; i < trackTree.GetEntries(); i++) {
    trackTree.GetEntry(i);
    if (trackTree.GetTreeNumber() >
        lastTracktreeID) {       // this part is needed for chained input
      if (lastTracktreeID > 0) { // new chunk, increase the offset
        trackoffs += trackTree.GetTree()->GetEntries();
      }
      lastTracktreeID = trackTree.GetTreeNumber();
    }
    int trackNumber = 0;
    for (const auto &trackLTF : trackLTFVec) {
      LOG(INFO) << "===== TrackLTF # " << trackNumber << " nPoints = " <<  trackLTF.getNPoints() << " =====";
      auto thisTrackMCCompLabels = trackLTF.getMCCompLabels();
      for (auto iLabel = 0; iLabel < trackLTF.getNPoints(); iLabel++) {
        LOG(INFO) << "MCCompLabels.TrackID = "
                  << thisTrackMCCompLabels[iLabel].getTrackID();
      }

      auto thisTrackLayers = trackLTF.getLayers();
      for (auto iLayer = 0; iLayer < trackLTF.getNPoints(); iLayer++) {
        LOG(INFO) << "trackLTF Layer = " << thisTrackLayers[iLayer];
      }

      auto thisXCoordinates = trackLTF.getXCoordinates();
      for (auto iCoordinate = 0; iCoordinate < trackLTF.getNPoints(); iCoordinate++) {
        LOG(INFO) << "trackLTF x = " << thisXCoordinates[iCoordinate];
      }

      trackNumber++;
    }
  } // loop over multiple ROFvectors (in case of chaining)

  //
  swTot.Stop();
  swTot.Print();
}
