#include "SubRunDataGetterHelper.h"
// vim: set sw=2 expandtab :

#include "canvas/Persistency/Provenance/BranchDescription.h"
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Persistency/Provenance/canonicalProductName.h"

#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/uniform_type_name.h"
#include "gallery/AssnsBranchData.h"

#include "canvas_root_io/Streamers/AssnsStreamer.h"
#include "canvas_root_io/Streamers/BranchDescriptionStreamer.h"
#include "canvas_root_io/Streamers/CacheStreamers.h"
#include "canvas_root_io/Streamers/ProductIDStreamer.h"
#include "canvas_root_io/Streamers/RefCoreStreamer.h"
#include "canvas_root_io/Streamers/TransientStreamer.h"
#include "canvas_root_io/Streamers/setPtrVectorBaseStreamer.h"

#include "SubRunHistoryGetter.h"
#include "SubRunNavigator.h"

#include "TClass.h"
#include "TTree.h"

#include <iostream>
#include <algorithm>
#include <cassert>

namespace {

  constexpr char underscore{'_'};
  constexpr char period{'.'};

  void
  initializeStreamers()
  {
    art::setCacheStreamers();
    art::setProvenanceTransientStreamers();
    art::detail::setBranchDescriptionStreamer();
    art::detail::setPtrVectorBaseStreamer();
    art::configureProductIDStreamer();
    art::configureRefCoreStreamer();
  }

  std::string
  buildBranchName(gallery::InfoForTypeLabelInstance const& info,
                  std::string const& processName)
  {
    std::string branchName{info.type().friendlyClassName()};
    unsigned int const branchNameSize =
      branchName.size() + info.label().size() + info.instance().size() +
      processName.size() + 4;
    branchName.reserve(branchNameSize);
    branchName += underscore;
    branchName += info.label();
    branchName += underscore;
    branchName += info.instance();
    branchName += underscore;
    branchName += processName;
    branchName += period;
    return branchName;
  }

} // unnamed namespace

namespace gallery {

  SubRunDataGetterHelper::SubRunDataGetterHelper(
    SubRunNavigator const* subrunNavigator,
    std::shared_ptr<SubRunHistoryGetter> historyGetter)
    : subrunNavigator_{subrunNavigator}, historyGetter_{historyGetter}
  {
    initializeStreamers();
  }

  ProductWithID
  SubRunDataGetterHelper::getByLabel(std::type_info const& typeInfoOfWrapper,
				     art::InputTag const& inputTag) const
  {
    if (!initializedForProcessHistory_) {
      initializeForProcessHistory();
    }
    art::TypeID const type{typeInfoOfWrapper};
    auto const& info =
      getInfoForTypeLabelInstance(type, inputTag.label(), inputTag.instance());

    if (inputTag.process().empty()) {
      // search in reverse order of the ProcessHistory
      for (auto reverseIter = crbegin(info.productIDsOrderedByHistory()),
                iEnd = crend(info.productIDsOrderedByHistory());
           reverseIter != iEnd;
           ++reverseIter) {
        if (auto edProduct = readProduct(*reverseIter, type)) {
          // If the product was present in the input file and we
          // successfully read it then we are done
          return std::make_pair(edProduct, *reverseIter);
        }
      }
    } else { // process is not empty
      auto itProcess = processNameToProcessIndex_.find(inputTag.process());
      if (itProcess != cend(processNameToProcessIndex_)) {
        unsigned int const processIndex = itProcess->second;
        auto branchData = getBranchData(info, processIndex);
        if (branchData) {
          auto pd = branchMapReader_.productDescription(info, itProcess->first);
          assert(pd);
          auto product = branchData->uniqueProduct_(type);
          return std::make_pair(product, pd->productID());
        }
      }
    }
    return std::make_pair(nullptr, art::ProductID::invalid());
  }

  std::vector<art::BranchDescription const*>
  SubRunDataGetterHelper::getProductDescriptions(
    art::TypeID const& typeIDOfWrapper) const
  {
    std::vector<art::BranchDescription const*> result;
    if (!initializedForProcessHistory_) {
      initializeForProcessHistory();
    }
    auto const fcn = typeIDOfWrapper.friendlyClassName();
    for (auto const& pr : branchMapReader_.productDescriptions()) {
      auto const& pd = pr.second;
      if (pd.friendlyClassName() != fcn) {
        continue;
      }

      result.push_back(&pd);
    }
    return result;
  }

  std::vector<art::InputTag>
  SubRunDataGetterHelper::getInputTags(std::type_info const& typeInfoOfWrapper) const
  {
    std::vector<art::InputTag> result;
    art::TypeID const type{typeInfoOfWrapper};
    for (auto const* pd : getProductDescriptions(type)) {
      result.push_back(pd->inputTag());
    }
    return result;
  }

  std::vector<ProductWithID>
  SubRunDataGetterHelper::getManyByType(std::type_info const& typeInfoOfWrapper) const
  {
    std::vector<ProductWithID> result;
    art::TypeID const type{typeInfoOfWrapper};
    for (auto const* pd : getProductDescriptions(type)) {
      auto itProcess = processNameToProcessIndex_.find(pd->processName());
      if (itProcess == cend(processNameToProcessIndex_)) {
        continue;
      }

      auto const& info = getInfoForTypeLabelInstance(
        type, pd->moduleLabel(), pd->productInstanceName());

      unsigned int const processIndex = itProcess->second;
      auto branchData = getBranchData(info, processIndex);
      if (!branchData) {
        continue;
      }

      if (auto product = branchData->uniqueProduct_(type)) {
        result.emplace_back(product, pd->productID());
      }
    }
    return result;
  }

  void
  SubRunDataGetterHelper::updateFile(TFile* iFile,
				     TTree* iTree,
				     bool const initializeTheCache)
  {
    std::cout << "SubRunDataGetterHelper::updateFile" << std::endl;
    tree_ = iTree;
    art::configureProductIDStreamer(); // Null out the ProductID streamer
    branchMapReader_.updateFile(iFile);
    art::configureProductIDStreamer(branchMapReader_.branchIDLists());
    if (initializeTheCache) {
      tree_->SetCacheSize();
      tree_->AddBranchToCache(subrunNavigator_->subrunAuxiliaryBranch(), kTRUE);
    }
    std::cout << "SubRunDataGetterHelper::updateFile - infoVector.size()=" << infoVector_.size() << std::endl;
    for (auto& info : infoVector_) {
      std::vector<IndexProductIDPair> old;
      old.swap(info.processIndexToProductID());
      info.processIndexToProductID().reserve(processNames_.size());
      for (unsigned int processIndex{}; processIndex < processNames_.size();
           ++processIndex) {
        auto const& processName = processNames_[processIndex];
        auto bd = branchMapReader_.productDescription(info, processName);
        if (bd == nullptr) {
          // Product not available.
          continue;
        }

        auto const productID = getMaybeValidProductID(old, processIndex);
        if (productID.isValid()) {
          auto& branchData = *branchDataMap_.at(productID);
          auto branch = tree_->GetBranch(branchData.branchName().c_str());
          // This will update the pointer to the TBranch in
          // BranchData.  This loop is sufficient to update all
          // BranchData objects in the vector.
          branchData.updateFile(branch);

          // A little paranoia here. What if in one input file
          // Assns<A,B,C> is written into the file and Assns<B,A,C> is
          // written into another? The following deals properly with
          // that case by constructing a new AssnsBranchData object,
          // although I imagine it might never happen in practice.
          if (info.isAssns() && branch != nullptr) {
            TClass* tClass = getTClass(info, processName);

            art::TypeID const typeIDInDescription{tClass->GetTypeInfo()};
            art::TypeID const typeIDInBranchData{
              branchData.tClass()->GetTypeInfo()};
            if (typeIDInDescription != typeIDInBranchData) {
              std::string branchName = branchData.branchName();
	      std::cout << "branchName=" << branchName << " : removed AssnsBranchData" << std::endl;
              // branchDataMap_[productID] =
              //   std::make_unique<AssnsBranchData>(typeIDInDescription,
              //                                     tClass,
              //                                     branch,
              //                                     subrunNavigator_,
              //                                     this,
              //                                     move(branchName),
              //                                     info.type(),
              //                                     info.partnerType());
            }
          }
          info.processIndexToProductID().emplace_back(processIndex,
                                                      bd->productID());
          if (initializeTheCache && branch) {
            tree_->AddBranchToCache(branch, kTRUE);
          }
        } else if (branchMapReader_.branchInRegistryOfAnyOpenedFile(
                     productID)) {
          addBranchData(
            bd->branchName(), processIndex, info, initializeTheCache);
        }
      }
      updateBranchDataIndexOrderedByHistory(info);
    }
    if (initializeTheCache) {
      tree_->StopCacheLearningPhase();
    }
    // This must be cleared because the BranchIDLists may be different
    // in different input files.
    branchDataMissingSet_.clear();
    updateEvent();
  }

  void
  SubRunDataGetterHelper::initializeTTreeCache()
  {
    tree_->SetCacheSize();
    tree_->AddBranchToCache(subrunNavigator_->subrunAuxiliaryBranch(), kTRUE);
    for (auto const& info : infoVector_) {
      for (auto const& i : info.processIndexToProductID()) {
        auto const productID = i.second;
        auto& branchData = *branchDataMap_.at(productID);
        if (auto branch = branchData.branch()) {
          tree_->AddBranchToCache(branch, kTRUE);
        }
      }
    }
    tree_->StopCacheLearningPhase();
  }

  void
  SubRunDataGetterHelper::updateEvent()
  {
    initializedForProcessHistory_ = false;
  }

  void
  SubRunDataGetterHelper::initializeForProcessHistory() const
  {
    initializedForProcessHistory_ = true;
    // Do nothing if the process names in the process history are the same
    auto const processHistoryID = historyGetter_->processHistoryID();
    if (processHistoryID == previousProcessHistoryID_) {
      return;
    }
    previousProcessHistoryID_ = processHistoryID;
    art::ProcessHistory const& processHistory =
      historyGetter_->processHistory();
    if (previousProcessHistoryNames_.size() == processHistory.size()) {
      bool same = true;
      auto iPrevious = previousProcessHistoryNames_.begin();
      for (auto i = processHistory.begin(), iEnd = processHistory.end();
           i != iEnd;
           ++i, ++iPrevious) {
        if (i->processName() != *iPrevious) {
          same = false;
          break;
        }
      }
      if (same) {
        return;
      }
    }
    previousProcessHistoryNames_.clear();
    // update for the new process history
    orderedProcessIndexes_.clear();
    for (auto const& processConfig : processHistory) {
      std::string const& processName = processConfig.processName();
      previousProcessHistoryNames_.push_back(processName);
      auto itFind = processNameToProcessIndex_.find(processName);
      if (itFind == processNameToProcessIndex_.end()) {
        addProcess(processName);
        itFind = processNameToProcessIndex_.find(processName);
      }
      orderedProcessIndexes_.push_back(itFind->second);
    }
    for (auto const& info : infoVector_) {
      updateBranchDataIndexOrderedByHistory(info);
    }
  }

  void
  SubRunDataGetterHelper::addProcess(std::string const& processName) const
  {
    unsigned int const processIndex = processNames_.size();
    processNames_.push_back(processName);
    processNameToProcessIndex_[processName] = processIndex;
    for (auto& info : infoVector_) {
      auto pd = branchMapReader_.productDescription(info, processName);
      if (pd == nullptr) {
        continue;
      }
      assert(branchMapReader_.branchInRegistryOfAnyOpenedFile(pd->productID()));
      addBranchData(pd->branchName(), processIndex, info);
    }
  }

  void
  SubRunDataGetterHelper::addBranchData(std::string branchName,
                                  unsigned int const processIndex,
                                  InfoForTypeLabelInstance const& info,
                                  bool const initializeTheCache) const
  {
    auto branch = tree_->GetBranch(branchName.c_str());
    if (branch == nullptr) {
      return;
    }
    auto pd =
      branchMapReader_.productDescription(info, processNames_[processIndex]);
    assert(pd);
    auto const productID = pd->productID();
    if (info.isAssns()) {
      auto const& processName = processNames_[processIndex];
      auto tClass = getTClass(info, processName);

      art::TypeID const typeIDInDescription{tClass->GetTypeInfo()};
      // branchDataMap_[productID] =
      //   std::make_unique<AssnsBranchData>(typeIDInDescription,
      //                                     tClass,
      //                                     branch,
      //                                     subrunNavigator_,
      //                                     this,
      //                                     move(branchName),
      //                                     info.type(),
      //                                     info.partnerType());
    } else {
      std::cout << "Fuck" << std::endl;
      // branchDataMap_[productID] =
      //   std::make_unique<BranchData>(info.type(),
      //                                info.tClass(),
      //                                branch,
      //                                subrunNavigator_,
      //                                this,
      //                                move(branchName));
    }
    info.processIndexToProductID().emplace_back(processIndex, productID);
    if (initializeTheCache && branch) {
      tree_->AddBranchToCache(branch, kTRUE);
    }
  }

  TClass*
  SubRunDataGetterHelper::getTClass(InfoForTypeLabelInstance const& info,
                              std::string const& processName) const
  {
    auto bd = branchMapReader_.productDescription(info, processName);
    if (bd == nullptr) {
      throw art::Exception{art::errors::LogicError,
                           "SubRunDataGetterHelper::getTClass: "}
        << "TBranch exists but no BranchDescription in ProductRegistry.\n"
        << "This shouldn't be possible. For type " << info.type().className()
        << "\n";
    }

    art::detail::AssnsStreamer::init_streamer(bd->producedClassName());
    auto tClass = TClass::GetClass(bd->wrappedName().c_str());
    if (tClass == nullptr) {
      throw art::Exception{art::errors::DictionaryNotFound,
                           "SubRunDataGetterHelper::getTClass: "}
        << "Missing dictionary for wrapped Assns class.\n"
        << bd->wrappedName() << "\n";
    }
    return tClass;
  }

  InfoForTypeLabelInstance&
  SubRunDataGetterHelper::getInfoForTypeLabelInstance(
    art::TypeID const& type,
    std::string const& label,
    std::string const& instance) const
  {
    if (label.empty()) {
      throw art::Exception(art::errors::LogicError)
        << "getValidHandle was passed an empty module label. Not allowed.\n";
    }
    TypeLabelInstanceKey const key{type, label, instance};
    auto itFind = infoMap_.find(key);
    if (itFind == cend(infoMap_)) {
      addTypeLabelInstance(type, label, instance);
      itFind = infoMap_.find(key);
    }
    return infoVector_[itFind->second];
  }

  void
  SubRunDataGetterHelper::addTypeLabelInstance(art::TypeID const& type,
                                         std::string const& label,
                                         std::string const& instance) const
  {
    dictChecker_.checkDictionaries(art::uniform_type_name(type.typeInfo()),
                                   true);
    dictChecker_.reportMissingDictionaries();
    unsigned int infoIndex = infoVector_.size();
    infoVector_.emplace_back(type, label, instance);
    insertIntoInfoMap(type, label, instance, infoIndex);
    InfoForTypeLabelInstance const& info = infoVector_[infoIndex];
    if (info.isAssns()) {
      insertIntoInfoMap(info.partnerType(), label, instance, infoIndex);
    }
    unsigned int processIndex{};
    for (auto const& processName : processNames_) {
      std::string branchName{buildBranchName(info, processName)};
      art::ProductID const productID{branchName};
      if (branchMapReader_.branchInRegistryOfAnyOpenedFile(productID)) {
        addBranchData(move(branchName), processIndex, info);
      }
      ++processIndex;
    }
    updateBranchDataIndexOrderedByHistory(info);
  }

  void
  SubRunDataGetterHelper::insertIntoInfoMap(art::TypeID const& type,
                                      std::string const& label,
                                      std::string const& instance,
                                      unsigned int const infoIndex) const
  {
    TypeLabelInstanceKey const newKey{type, label, instance};
    infoMap_[newKey] = infoIndex;
  }

  art::EDProduct const*
  SubRunDataGetterHelper::readProduct(art::ProductID const productID,
                                art::TypeID const& type) const
  {
    auto const& branchData = branchDataMap_.at(productID);
    return branchData->uniqueProduct_(type);
  }

  art::ProductID
  SubRunDataGetterHelper::getMaybeValidProductID(
    std::vector<IndexProductIDPair> const& processIndexToProductID,
    unsigned int const processIndex) const
  {
    auto itBranchDataIndex = lower_bound(
      processIndexToProductID.cbegin(),
      processIndexToProductID.cend(),
      std::make_pair(processIndex, art::ProductID::invalid()),
      [](auto const& l, auto const& r) { return l.first < r.first; });
    if (itBranchDataIndex != processIndexToProductID.cend() &&
        itBranchDataIndex->first == processIndex) {
      return itBranchDataIndex->second;
    }
    return art::ProductID::invalid();
  }

  BranchData const*
  SubRunDataGetterHelper::getBranchData(InfoForTypeLabelInstance const& info,
                                  unsigned int const processIndex) const
  {
    auto const pid =
      getMaybeValidProductID(info.processIndexToProductID(), processIndex);
    return pid.isValid() ? branchDataMap_.at(pid).get() : nullptr;
  }

  void
  SubRunDataGetterHelper::updateBranchDataIndexOrderedByHistory(
    InfoForTypeLabelInstance const& info) const
  {
    info.productIDsOrderedByHistory().clear();
    if (info.productIDsOrderedByHistory().capacity() <
        orderedProcessIndexes_.size()) {
      info.productIDsOrderedByHistory().reserve(orderedProcessIndexes_.size());
    }
    for (auto processIndex : orderedProcessIndexes_) {
      auto const& orderedProcessName = processNames_[processIndex];
      auto pd = branchMapReader_.productDescription(info, orderedProcessName);
      if (pd == nullptr) {
        continue;
      }
      auto branchData = getBranchData(info, processIndex);
      if (branchData && branchData->branch() != nullptr) {
        info.productIDsOrderedByHistory().emplace_back(pd->productID());
      }
    }
  }

  BranchData const*
  SubRunDataGetterHelper::getBranchData(art::BranchDescription const& desc) const
  {
    if (!initializedForProcessHistory_) {
      initializeForProcessHistory();
    }
    TClass* tClass = TClass::GetClass(desc.wrappedName().c_str());
    if (tClass == nullptr) {
      throw art::Exception(art::errors::DictionaryNotFound,
                           "SubRunDataGetterHelper::getBranchData: ")
        << "Missing dictionary for wrapped class.\n"
        << desc.wrappedName() << "\n";
    }

    auto const& typeInfoOfWrapper = *tClass->GetTypeInfo();
    art::TypeID const type{typeInfoOfWrapper};
    auto const& info = getInfoForTypeLabelInstance(
      type, desc.moduleLabel(), desc.productInstanceName());
    auto itProcess = processNameToProcessIndex_.find(desc.processName());
    if (itProcess != cend(processNameToProcessIndex_)) {
      unsigned int const processIndex = itProcess->second;
      auto branchData = getBranchData(info, processIndex);
      if (branchData && branchData->branch() != nullptr) {
        return branchData;
      }
    }
    return nullptr;
  }

  art::BranchDescription const&
  SubRunDataGetterHelper::getProductDescription(art::ProductID const productID) const
  {
    auto pd = branchMapReader_.productDescription(productID);
    if (pd == nullptr) {
      throw art::Exception{art::errors::ProductNotFound,
                           "SubRunDataGetterHelper::getProductDescription: "}
        << "No product description could be found for ProductID " << productID
        << ".\n";
    }
    return *pd;
  }

  art::EDProductGetter const*
  SubRunDataGetterHelper::getEDProductGetter_(art::ProductID const& productID) const
  {
    auto itFind = branchDataMap_.find(productID);
    if (itFind != cend(branchDataMap_)) {
      return itFind->second.get();
    }

    if (branchDataMissingSet_.find(productID) != cend(branchDataMissingSet_)) {
      return &invalidBranchData_;
    }

    auto const pd = branchMapReader_.productDescription(productID);
    if (pd == nullptr) {
      return &invalidBranchData_;
    }

    auto branchData = getBranchData(*pd);
    if (branchData == nullptr) {
      branchDataMissingSet_.insert(productID);
      return &invalidBranchData_;
    }
    return branchData;
  }

} // namespace gallery
