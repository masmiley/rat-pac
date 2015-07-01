/**
 * @class RAT::DS::MCSample
 * @brief A sampled charge/time for a hit on an electronics channel
 */

#ifndef __RAT_DS_MCSample__
#define __RAT_DS_MCSample__

#include <Rtypes.h>
#include <TObject.h>

namespace RAT {
  namespace DS {

class MCSample : public TObject {
public:
  MCSample() {}
  virtual ~MCSample() {}

  virtual Double_t GetHitTime() const { return hitTime; }
  virtual void SetHitTime(Double_t _hitTime) { hitTime = _hitTime; }

  virtual Double_t GetCharge() const { return charge; }
  virtual void SetCharge(Double_t _charge) { charge = _charge; }

private:
  Double_t hitTime;
  Double_t charge;

  ClassDef(MCSample, 1)
};

  }  // namespace DS
}  // namespace RAT

#endif  // __RAT_DS_MCSample__

