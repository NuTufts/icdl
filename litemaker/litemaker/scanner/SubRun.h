#ifndef gallery_SubRun_h
#define gallery_SubRun_h

// ====================================================================
// Main interface to users. It uses the DataGetterHelper and
// SubRunNavigator to iterate over subruns in a set of input files and
// find products in them.
// ====================================================================

#include "canvas/Persistency/Common/EDProduct.h"
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Provenance/ProcessHistoryID.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Persistency/Provenance/ProductToken.h"
#include "canvas/Persistency/Provenance/ProvenanceFwd.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib/container_algorithms.h"
#include "gallery/Handle.h"
#include "gallery/ValidHandle.h"

#include "SubRunNavigator.h"
#include "SubRunDataGetterHelper.h"

#include <cassert>
#include <memory>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

class TFile;
class TTree;

namespace gallery {

  // The gallery::SubRun provides read-only access to subrun data
  // products. It also provides the means to iterate through one or
  // more art/ROOT files, reading the subruns in those files; in this
  // sense, the gallery::SubRun is a kind of iterator. Only a valid
  // SubRun (as defined by SubRun::isValid()) can be used to access
  // subrun data products.
  class SubRun {
  public:
    // gallery::SubRun objects can be neither copied nor assigned. They
    // can be move-copied and move-assigned.

    // Construct an SubRun that will read subruns from each of the
    // named art/ROOT files, in order.
    explicit SubRun(std::vector<std::string> const& fileNames,
                   bool useTTreeCache = true,
                   unsigned int subrunsToLearnUsedBranches = 7);

    // Get access to a data product of type PROD, using a ValidHandle.
    template <typename PROD>
    gallery::ValidHandle<PROD> getValidHandle(art::InputTag const&) const;

    // Get access to a data product of type PROD, using a Handle.
    // Only if the return value is 'true' is the Handle valid.
    template <typename PROD>
    bool getByLabel(art::InputTag const&, Handle<PROD>& result) const;

    // Get access to all data products of type PROD, using a
    // std::vector<Handle>.
    template <typename PROD>
    void getManyByType(std::vector<Handle<PROD>>& result) const;

    // Return all input tags corresponding to products of type PROD.
    // This function call *does not* read or provide access to any
    // products.
    template <typename PROD>
    std::vector<art::InputTag> getInputTags() const;

    template <typename PROD>
    std::vector<art::ProductToken<PROD>> getProductTokens() const;

    art::SubRunAuxiliary const& subrunAuxiliary() const;

    // Return the product description if it is present.
    art::BranchDescription const& getProductDescription(art::ProductID) const;

    // Return the number of subruns in the currently-open file.
    long long numberOfSubRunsInFile() const;

    // Return the current entry number (the entry number of the SubRuns
    // tree in the current art/ROOT file).
    long long subrunEntry() const;

    // Return the index of the current file.
    long long fileEntry() const;

    // Go to the entry with the given index. If the index is
    // out-of-bounds, or negative, or if the SubRun is not suitable for
    // random access, an exception will be thrown.
    void goToEntry(long long entry);

    // Return true if the SubRun can be used to access data products.
    bool isValid() const;

    // Return true if we are at the end of the sequence of subruns
    // through which we are iterating.
    bool atEnd() const;

    // Go to the first subrun of the sequence we are to traverse.
    void toBegin();
    void first();

    // Go to the next subrun in the sequence.
    SubRun& operator++();
    void next();

    // Throws an exception if the SubRun was constructed with more than
    // one filename, or if we are already at the beginning of the
    // sequence. Otherwise, go to the previous subrun in the sequence.
    SubRun& operator--();
    void previous();

    TFile* getTFile() const;
    TTree* getTTree() const;

    template <typename T>
    using HandleT = Handle<T>;

  private:
    ProductWithID getByLabel(std::type_info const& typeInfoOfWrapper,
                             art::InputTag const& inputTag) const;

    std::vector<ProductWithID> getManyByType(
      std::type_info const& typeInfoOfWrapper) const;

    [[noreturn]] void throwProductNotFoundException(
      std::type_info const& typeInfo,
      art::InputTag const& tag) const;

    std::shared_ptr<art::Exception const> makeProductNotFoundException(
      std::type_info const& typeInfo,
      art::InputTag const& tag) const;

    void checkForEnd() const;
    void updateAfterSubRunChange(long long oldFileEntry);

    bool randomAccessOK_;
    std::unique_ptr<SubRunNavigator> subrunNavigator_;
    std::unique_ptr<SubRunDataGetterHelper> dataGetterHelper_;

    bool useTTreeCache_;
    unsigned int subrunsToLearnUsedBranches_;
    unsigned int subrunsProcessed_{};
  };
} // namespace gallery

template <typename PROD>
inline gallery::ValidHandle<PROD>
gallery::SubRun::getValidHandle(art::InputTag const& inputTag) const
{
  checkForEnd();
  std::type_info const& typeInfoOfWrapper{typeid(art::Wrapper<PROD>)};
  auto res = getByLabel(typeInfoOfWrapper, inputTag);
  auto edProduct = res.first;

  auto ptrToWrapper = dynamic_cast<art::Wrapper<PROD> const*>(edProduct);
  if (ptrToWrapper == nullptr) {
    throwProductNotFoundException(typeid(PROD), inputTag);
  }

  auto product = ptrToWrapper->product();
  return ValidHandle<PROD>{product, res.second};
}

template <typename PROD>
inline bool
gallery::SubRun::getByLabel(art::InputTag const& inputTag,
                           Handle<PROD>& result) const
{
  checkForEnd();
  if (inputTag.empty()) {
    result = Handle<PROD>{makeProductNotFoundException(typeid(PROD), inputTag)};
    return false;
  }

  std::type_info const& typeInfoOfWrapper{typeid(art::Wrapper<PROD>)};
  auto res = getByLabel(typeInfoOfWrapper, inputTag);
  auto edProduct = res.first;

  auto ptrToWrapper = dynamic_cast<art::Wrapper<PROD> const*>(edProduct);

  if (ptrToWrapper == nullptr) {
    result = Handle<PROD>{makeProductNotFoundException(typeid(PROD), inputTag)};
    return false;
  }
  auto product = ptrToWrapper->product();
  result = Handle<PROD>{product, res.second};
  return true;
}

template <typename PROD>
inline void
gallery::SubRun::getManyByType(std::vector<Handle<PROD>>& result) const
{
  std::type_info const& typeInfoOfWrapper{typeid(art::Wrapper<PROD>)};
  auto products = getManyByType(typeInfoOfWrapper);
  std::vector<Handle<PROD>> tmp;
  cet::transform_all(products, back_inserter(tmp), [](auto const& pr) {
    auto product = pr.first;
    auto wrapped_product = dynamic_cast<art::Wrapper<PROD> const*>(product);
    assert(wrapped_product != nullptr);
    auto user_product = wrapped_product->product();
    assert(user_product != nullptr);
    return Handle<PROD>{user_product, pr.second};
  });
  swap(tmp, result);
}

template <typename PROD>
inline std::vector<art::InputTag>
gallery::SubRun::getInputTags() const
{
  std::type_info const& typeInfoOfWrapper{typeid(art::Wrapper<PROD>)};
  return dataGetterHelper_->getInputTags(typeInfoOfWrapper);
}

template <typename PROD>
inline std::vector<art::ProductToken<PROD>>
gallery::SubRun::getProductTokens() const
{
  std::vector<art::ProductToken<PROD>> result;
  auto const tags = getInputTags<PROD>();
  cet::transform_all(tags, back_inserter(result), [](auto const& tag) {
    return art::ProductToken<PROD>{tag};
  });
  return result;
}

inline void
gallery::SubRun::first()
{
  return toBegin();
}

#endif /* gallery_SubRun_h */

// Local Variables:
// mode: c++
// End:
