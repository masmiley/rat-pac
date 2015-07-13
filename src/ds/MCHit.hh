/**
 * @class RAT::DS::MCHit
 * @brief A hit electronics channel
 */

#ifndef __RAT_DS_MCHit__
#define __RAT_DS_MCHit__

#include <Rtypes.h>
#include <TObject.h>
#include <vector>
#include <RAT/DS/MCSample.hh>

namespace RAT {
  namespace DS {

class MCHit : public TObject {
public:
  MCHit() {}
  virtual ~MCHit() {}

  virtual Int_t GetPMTID() const { return pmtID; }
  virtual void SetPMTID(Int_t _pmtID) { pmtID = _pmtID; }

  virtual Int_t GetCrate() const { return crate; }
  virtual void SetCrate(Int_t _crate) { crate = _crate; }

  virtual Int_t GetCard() const { return card; }
  virtual void SetCard(Int_t _card) { card = _card; }

  virtual Int_t GetChannel() const { return channel; }
  virtual void SetChannel(Int_t _channel) { channel = _channel; }

  virtual MCSample* GetMCSample(Int_t i) { return &sample.at(i); }
  virtual Int_t GetMCSampleCount() const { return sample.size(); }
  virtual MCSample* AddNewMCSample() {
    sample.resize(sample.size() + 1);
    return &sample.back();
  }
  virtual void PruneMCSample() { sample.resize(0); }

private:
  // ID number of hit PMT
  Int_t pmtID;

  // Electronics number of hit channel
  Int_t crate;
  Int_t card;
  Int_t channel;

  std::vector<MCSample> sample;
  
  ClassDef(MCHit, 1)
};

  }  // namespace DS
}  // namespace RAT

#endif  // __RAT_DS_MCHit__

