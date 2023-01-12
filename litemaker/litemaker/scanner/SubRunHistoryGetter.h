#ifndef gallery_SubRunHistoryGetter_h
#define gallery_SubRunHistoryGetter_h

#include "canvas/Persistency/Provenance/ProcessHistoryID.h"

namespace art {
  class History;
  class ProcessHistory;
} // namespace art

namespace gallery {

  class SubRunNavigator;

  class SubRunHistoryGetter {
  public:
    SubRunHistoryGetter(SubRunNavigator const*);
    virtual ~SubRunHistoryGetter() = default;

    SubRunHistoryGetter(SubRunHistoryGetter const&) = delete;
    SubRunHistoryGetter const& operator=(SubRunHistoryGetter const&) = delete;

    virtual art::ProcessHistoryID const& processHistoryID() const;
    virtual art::ProcessHistory const& processHistory() const;
    virtual art::History const& history() const;

  private:
    SubRunNavigator const* subrunNavigator_;
  };
} // namespace gallery
#endif /* gallery_SubRunHistoryGetter_h */

// Local Variables:
// mode: c++
// End:
