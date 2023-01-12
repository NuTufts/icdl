#include "SubRun.h"

#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/TypeID.h"
#include "gallery/throwFunctions.h"

#include "SubRunHistoryGetter.h"
#include "SubRunNavigator.h"
#include "SubRunDataGetterHelper.h"

namespace gallery {

  SubRun::SubRun(std::vector<std::string> const& fileNames,
		 bool const useTTreeCache,
		 unsigned int const subrunsToLearnUsedBranches)
    : randomAccessOK_{fileNames.size() == 1}
    , subrunNavigator_{std::make_unique<SubRunNavigator>(fileNames)}
    , dataGetterHelper_{std::make_unique<SubRunDataGetterHelper>(
        subrunNavigator_.get(),
        std::make_shared<SubRunHistoryGetter>(subrunNavigator_.get()))}
    , useTTreeCache_{useTTreeCache}
    , subrunsToLearnUsedBranches_{subrunsToLearnUsedBranches}
  {

    if (subrunsToLearnUsedBranches_ < 1)
      subrunsToLearnUsedBranches_ = 1;
    if (!atEnd()) {
      bool constexpr initializeTheCache{false};
      dataGetterHelper_->updateFile(subrunNavigator_->getTFile(),
                                    subrunNavigator_->getTTree(),
                                    initializeTheCache);
    }
  }

  art::SubRunAuxiliary const&
  SubRun::subrunAuxiliary() const
  {
    return subrunNavigator_->subrunAuxiliary();
  }

  art::BranchDescription const&
  SubRun::getProductDescription(art::ProductID const pid) const
  {
    return dataGetterHelper_->getProductDescription(pid);
  }

  long long
  SubRun::numberOfSubRunsInFile() const
  {
    return subrunNavigator_->entriesInCurrentFile();
  }

  long long
  SubRun::subrunEntry() const
  {
    return subrunNavigator_->subrunEntry();
  }

  long long
  SubRun::fileEntry() const
  {
    return subrunNavigator_->fileEntry();
  }

  void
  SubRun::goToEntry(long long const entry)
  {
    if (!randomAccessOK_)
      throwIllegalRandomAccess();
    if (entry < 0)
      throwIllegalRandomAccess();
    if (entry >= numberOfSubRunsInFile())
      throwIllegalRandomAccess();
    subrunNavigator_->goToEntry(entry);
  }

  bool
  SubRun::isValid() const
  {
    return subrunNavigator_->isValid();
  }

  bool
  SubRun::atEnd() const
  {
    return subrunNavigator_->atEnd();
  }

  void
  SubRun::toBegin()
  {
    long long const oldSubRunEntry = subrunEntry();
    long long const oldFileEntry = fileEntry();
    subrunNavigator_->toBegin();
    if (oldSubRunEntry == subrunEntry() && oldFileEntry == fileEntry()) {
      return;
    }
    updateAfterSubRunChange(oldFileEntry);
  }

  SubRun&
  SubRun::operator++()
  {
    auto const oldFileEntry = fileEntry();
    subrunNavigator_->next();
    updateAfterSubRunChange(oldFileEntry);
    return *this;
  }

  SubRun&
  SubRun::operator--()
  {
    if (!randomAccessOK_)
      throwIllegalRandomAccess();
    if (atEnd())
      throwIllegalDecrement();
    auto const oldFileEntry = fileEntry();
    subrunNavigator_->previous();
    updateAfterSubRunChange(oldFileEntry);
    return *this;
  }

  void
  SubRun::updateAfterSubRunChange(long long const oldFileEntry)
  {
    ++subrunsProcessed_;
    if (atEnd())
      return;
    if (oldFileEntry != fileEntry()) {
      bool initializeTheCache =
        useTTreeCache_ && subrunsProcessed_ >= subrunsToLearnUsedBranches_;
      dataGetterHelper_->updateFile(subrunNavigator_->getTFile(),
                                    subrunNavigator_->getTTree(),
                                    initializeTheCache);
    } else {
      dataGetterHelper_->updateEvent();
      if (useTTreeCache_ && subrunsProcessed_ == subrunsToLearnUsedBranches_) {
        dataGetterHelper_->initializeTTreeCache();
      }
    }
  }

  void
  SubRun::next()
  {
    operator++();
  }

  void
  SubRun::previous()
  {
    operator--();
  }

  TFile*
  SubRun::getTFile() const
  {
    return subrunNavigator_->getTFile();
  }

  TTree*
  SubRun::getTTree() const
  {
    return subrunNavigator_->getTTree();
  }

  ProductWithID
  SubRun::getByLabel(std::type_info const& typeInfoOfWrapper,
                    art::InputTag const& inputTag) const
  {
    return dataGetterHelper_->getByLabel(typeInfoOfWrapper, inputTag);
  }

  std::vector<ProductWithID>
  SubRun::getManyByType(std::type_info const& typeInfoOfWrapper) const
  {
    return dataGetterHelper_->getManyByType(typeInfoOfWrapper);
  }

  void
  SubRun::throwProductNotFoundException(std::type_info const& typeInfo,
                                       art::InputTag const& tag) const
  {
    auto e = makeProductNotFoundException(typeInfo, tag);
    throw *e;
  }

  std::shared_ptr<art::Exception const>
  SubRun::makeProductNotFoundException(std::type_info const& typeInfo,
                                      art::InputTag const& tag) const
  {
    auto e = std::make_shared<art::Exception>(art::errors::ProductNotFound);
    *e << "Failed to find product for \n  type = '"
       << art::TypeID{typeInfo}.className() << "'\n  module = '" << tag.label()
       << "'\n  productInstance = '"
       << ((!tag.instance().empty()) ? tag.instance().c_str() : "")
       << "'\n  process='"
       << ((!tag.process().empty()) ? tag.process().c_str() : "") << "'\n";
    return e;
  }

  void
  SubRun::checkForEnd() const
  {
    if (atEnd()) {
      throw art::Exception(art::errors::LogicError)
        << "You have requested data past the last subrun\n";
    }
  }
} // namespace gallery
