#ifndef gallery_SubRunNavigator_h
#define gallery_SubRunNavigator_h

// Manages iteration over the subruns in a vector of
// input files.  It handles the iteration, opening
// the files, opening TTrees, reading SubRunAuxilliary
// objects and reading History objects from the input
// files.

#include "canvas/Persistency/Provenance/SubRunAuxiliary.h"
#include "canvas/Persistency/Provenance/History.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"
#include "canvas/Persistency/Provenance/ProcessHistoryID.h"

#include "TFile.h"

#include <memory>
#include <string>
#include <vector>

class TBranch;
class TTree;

namespace gallery {

  class SubRunNavigator {
  public:
    explicit SubRunNavigator(std::vector<std::string> const& iFileNames);

    // In a normal iteration using the next function, isValid and
    // atEnd will always return opposite values. If nextFile was
    // called directly and there was an empty file, they might
    // return different values.

    bool
    isValid() const
    {
      return fileEntry_ != numberOfFiles_ &&
             subrunEntry_ != entriesInCurrentFile_;
    }

    bool
    atEnd() const
    {
      return fileEntry_ == numberOfFiles_;
    }

    void toBegin();
    void next();
    void previous();
    void goToEntry(long long entry);
    void nextFile();

    art::SubRunAuxiliary const& subrunAuxiliary() const;
    art::History const& history() const;
    art::ProcessHistoryID const& processHistoryID() const;
    art::ProcessHistory const& processHistory() const;

    TFile* getTFile() const;
    TTree* getTTree() const;
    TBranch*
    subrunAuxiliaryBranch() const
    {
      return subrunAuxiliaryBranch_;
    }

    long long
    fileEntry() const
    {
      return fileEntry_;
    }
    long long
    entriesInCurrentFile() const
    {
      return entriesInCurrentFile_;
    }
    long long
    subrunEntry() const
    {
      return subrunEntry_;
    }

  private:
    void initializeTTreePointers();
    void initializeTBranchPointers();

    std::vector<std::string> fileNames_;
    long long numberOfFiles_;
    long long fileEntry_{-1};
    long long firstFileWithSubRun_{};

    long long entriesInCurrentFile_{};
    long long subrunEntry_{};

    std::unique_ptr<TFile> file_{nullptr};

    TTree* subrunsTree_{nullptr};
    TBranch* subrunAuxiliaryBranch_{nullptr};
    mutable art::SubRunAuxiliary subrunAuxiliary_{};
    art::SubRunAuxiliary* pSubRunAuxiliary_{&subrunAuxiliary_};
    mutable long long previousSubRunAuxiliaryEntry_{-1};

    TTree* subrunHistoryTree_{nullptr};
    TBranch* subrunHistoryBranch_{nullptr};
    mutable art::History subrunHistory_{};
    art::History* pSubRunHistory_{&subrunHistory_};
    mutable long long previousSubRunHistoryEntry_{-1};

    mutable art::ProcessHistoryMap historyMap_{};
  };
} // namespace gallery

#endif /* gallery_SubRunNavigator_h */

// Local Variables:
// mode: c++
// End:
