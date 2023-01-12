#include "SubRunHistoryGetter.h"
#include "SubRunNavigator.h"

#include <type_traits>

namespace gallery {

  SubRunHistoryGetter::SubRunHistoryGetter(SubRunNavigator const* subrunNavigator)
    : subrunNavigator_(subrunNavigator)
  {
    static_assert(std::is_nothrow_destructible<SubRunHistoryGetter>::value,
                  "SubRunHistoryGetter is not nothrow destructible");
  }

  art::ProcessHistoryID const&
  SubRunHistoryGetter::processHistoryID() const
  {
    return subrunNavigator_->processHistoryID();
  }

  art::ProcessHistory const&
  SubRunHistoryGetter::processHistory() const
  {
    return subrunNavigator_->processHistory();
  }

  art::History const&
  SubRunHistoryGetter::history() const
  {
    return subrunNavigator_->history();
  }
} // namespace gallery
