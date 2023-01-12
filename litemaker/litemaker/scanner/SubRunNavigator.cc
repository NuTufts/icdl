#include "SubRunNavigator.h"

#include "canvas/Persistency/Provenance/BranchType.h"
#include "canvas/Persistency/Provenance/rootNames.h"
#include "canvas/Utilities/Exception.h"
#include "gallery/throwFunctions.h"

#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>

namespace gallery {

  SubRunNavigator::SubRunNavigator(std::vector<std::string> const& iFileNames)
    : fileNames_{iFileNames}, numberOfFiles_(fileNames_.size())
  {
    if (fileNames_.empty()) {
      fileEntry_ = 0;
    } else {
      while (subrunEntry_ == entriesInCurrentFile_ &&
             fileEntry_ != numberOfFiles_) {
        nextFile();
      }
      firstFileWithSubRun_ = fileEntry_;
    }
    if (isValid()) {
      subrunsTree_->LoadTree(subrunEntry_);
      std::cout << "SubRunNavigator: valid tree. nentries=" << subrunsTree_->GetEntries() << std::endl;
      subrunsTree_->Print();
    }
  }

  void
  SubRunNavigator::toBegin()
  {
    if (fileEntry_ == firstFileWithSubRun_) {
      subrunEntry_ = 0;
      return;
    }
    fileEntry_ = firstFileWithSubRun_ - 1;
    nextFile();
    if (isValid()) {
      subrunsTree_->LoadTree(subrunEntry_);
    }
  }

  void
  SubRunNavigator::next()
  {
    if (subrunEntry_ != entriesInCurrentFile_) {
      // normal case, go to next subrun in same file
      ++subrunEntry_;
    } else {
      // handle odd cases, calling next when already at the end
      // or after a call to nextFile() left us on a file with
      // no subruns. (Neither of these cases should occur in a
      // normal loop over subruns)
      nextFile();
    }
    // if we hit the end of the file go to the next file
    // also skip empty files
    while (subrunEntry_ == entriesInCurrentFile_ &&
           fileEntry_ != numberOfFiles_) {
      nextFile();
    }
    if (isValid()) {
      subrunsTree_->LoadTree(subrunEntry_);
    }
  }

  void
  SubRunNavigator::previous()
  {
    // We will never get here is we're dealing with more than one file.
    --subrunEntry_;
    if (isValid()) {
      subrunsTree_->LoadTree(subrunEntry_);
    }
  }

  void
  SubRunNavigator::goToEntry(long long entry)
  {
    subrunEntry_ = entry;
    if (isValid()) {
      subrunsTree_->LoadTree(subrunEntry_);
    }
  }

  void
  SubRunNavigator::nextFile()
  {
    // Be careful with this function. If the next file is empty this
    // will leave it not pointing at a valid subrun.
    if (atEnd()) {
      throw art::Exception{art::errors::LogicError}
        << "Illegal call to SubRunNavigator::nextFile() when atEnd() is true";
    }

    if (file_) {
      std::cout << "Closing file, read " << file_->GetBytesRead()
                << " bytes in " << file_->GetReadCalls() << " transactions\n";
    }

    ++fileEntry_;
    if (atEnd()) {
      entriesInCurrentFile_ = 0;
      subrunEntry_ = 0;
      file_ = nullptr;
      subrunsTree_ = nullptr;
      subrunAuxiliaryBranch_ = nullptr;
      subrunHistoryTree_ = nullptr;
      subrunHistoryBranch_ = nullptr;
      previousSubRunAuxiliaryEntry_ = -1;
      previousSubRunHistoryEntry_ = -1;
      historyMap_.clear();
      return;
    }
    file_.reset(TFile::Open(fileNames_[fileEntry_].c_str()));
    if (!file_ || file_->IsZombie()) {
      throw art::Exception(art::errors::FileOpenError)
        << "Failed opening file \'" << fileNames_[fileEntry_] << "\'";
    }
    std::cout << "Successfully opened file " << fileNames_[fileEntry_]
              << std::endl;
    initializeTTreePointers();
    initializeTBranchPointers();
    entriesInCurrentFile_ = subrunsTree_->GetEntries();
    subrunEntry_ = 0;
    previousSubRunAuxiliaryEntry_ = -1;
    previousSubRunHistoryEntry_ = -1;
    historyMap_.clear();
    return;
  }

  art::SubRunAuxiliary const&
  SubRunNavigator::subrunAuxiliary() const
  {
    if (previousSubRunAuxiliaryEntry_ != subrunEntry_) {
      subrunAuxiliaryBranch_->GetEntry(subrunEntry_);
      previousSubRunAuxiliaryEntry_ = subrunEntry_;
    }
    return subrunAuxiliary_;
  }

  art::History const&
  SubRunNavigator::history() const
  {
    if (previousSubRunHistoryEntry_ != subrunEntry_) {
      subrunHistoryBranch_->GetEntry(subrunEntry_);
      previousSubRunHistoryEntry_ = subrunEntry_;
    }
    return subrunHistory_;
  }

  art::ProcessHistoryID const&
  SubRunNavigator::processHistoryID() const
  {
    if (previousSubRunHistoryEntry_ != subrunEntry_) {
      subrunHistoryBranch_->GetEntry(subrunEntry_);
      previousSubRunHistoryEntry_ = subrunEntry_;
    }
    return subrunHistory_.processHistoryID();
  }

  art::ProcessHistory const&
  SubRunNavigator::processHistory() const
  {

    if (historyMap_.empty()) {

      std::unique_ptr<TTree> metaDataTree{
        file_->Get<TTree>(art::rootNames::metaDataTreeName().c_str())};

      if (!metaDataTree) {
        throwTreeNotFound(art::rootNames::metaDataTreeName());
      }

      auto pHistMapPtr = &historyMap_;
      TBranch* processHistoryBranch = metaDataTree->GetBranch(
        art::rootNames::metaBranchRootName<art::ProcessHistoryMap>());

      if (!processHistoryBranch) {
        throwBranchNotFound(
          art::rootNames::metaBranchRootName<art::ProcessHistoryMap>());
      }
      processHistoryBranch->SetAddress(&pHistMapPtr);
      processHistoryBranch->GetEntry(0);
    }
    return historyMap_[processHistoryID()];
  }

  TFile*
  SubRunNavigator::getTFile() const
  {
    return file_.get();
  }

  TTree*
  SubRunNavigator::getTTree() const
  {
    return subrunsTree_;
  }

  void
  SubRunNavigator::initializeTTreePointers()
  {
    subrunsTree_ = file_->Get<TTree>(art::rootNames::dataTreeName(art::InSubRun).c_str());
    if (subrunsTree_ == nullptr) {
      throwTreeNotFound(art::rootNames::dataTreeName(art::InSubRun));
    }
    if (subrunsTree_->GetEntries() < 0) {
      throw art::Exception(art::errors::LogicError)
        << "Unable to get the number of entries in subruns TTree.\n"
           "This might be a corrupted file.\n";
    }
    subrunHistoryTree_ =
      file_->Get<TTree>(art::rootNames::eventHistoryTreeName().c_str());
    if (subrunHistoryTree_ == nullptr) {
      throwTreeNotFound(art::rootNames::eventHistoryTreeName());
    }
  }

  void
  SubRunNavigator::initializeTBranchPointers()
  {
    subrunAuxiliaryBranch_ = subrunsTree_->GetBranch(
      art::BranchTypeToAuxiliaryBranchName(art::InSubRun).c_str());

    if (subrunAuxiliaryBranch_ == nullptr) {
      throwBranchNotFound(art::BranchTypeToAuxiliaryBranchName(art::InSubRun));
    }
    subrunAuxiliaryBranch_->SetAddress(&pSubRunAuxiliary_);

    subrunHistoryBranch_ = subrunHistoryTree_->GetBranch(
      art::rootNames::eventHistoryBranchName().c_str());

    if (subrunHistoryBranch_ == nullptr) {
      throwBranchNotFound(art::rootNames::eventHistoryBranchName());
    }
    subrunHistoryBranch_->SetAddress(&pSubRunHistory_);
  }
} // namespace gallery
