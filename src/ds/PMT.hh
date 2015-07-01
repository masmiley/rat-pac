/**
 * @class PMT
 * Data Structure: PMT in triggered event
 * 
 * This represents a PMT in a detector event.
 */

#ifndef __RAT_DS_PMT__
#define __RAT_DS_PMT__

#include <Rtypes.h>
#include <RAT/DS/Sample.hh>

namespace RAT {
  namespace DS {

class PMT : public TObject {
public:
  PMT() : TObject() {}
  virtual ~PMT() {}

  /** ID number of PMT */
  virtual void SetID(Int_t _id) { this->id = _id; }
  virtual Int_t GetID() { return id; }

  /** Conventience method to get total charge in the waveform */
  virtual inline Float_t GetCharge();

  /** Digitized charge and time samples */
  virtual Sample* GetSample(Int_t i) { return &sample[i]; }
  virtual Int_t GetSampleCount() const { return sample.size(); }
  virtual Sample* AddNewSample() {
    sample.resize(sample.size() + 1);
    return &sample.back();
  }
  virtual void PruneSample() { sample.resize(0); }

  ClassDef(PMT, 1);

protected:
  Int_t id;
  Float_t charge;
  Float_t time;
  std::vector<Sample> sample;
};


inline Float_t PMT::GetCharge() {
  Float_t charge = 0;
  std::vector<Sample>::iterator it;
  for (it=sample.begin(); it!=sample.end(); it++) {
    charge += it->GetCharge();
  }
  return charge;
}

  } // namespace DS
} // namespace RAT

#endif

