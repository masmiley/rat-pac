/**
 * @class DS::MC
 * Data Structure: Monte Carlo information
 *
 * @author Stan Seibert <sseibert@hep.upenn.edu>
 *
 * This structure contains the "Monte Carlo Truth" for the event.
 * Intermediate simulated quantities, like electronics response, can
 * also go here.  Only once the data finally looks like "real life"
 * should it go into the EV class.
 */

#ifndef __RAT_DS_MC__
#define __RAT_DS_MC__

#include <vector>
#include <algorithm>
#include <TObject.h>
#include <TTimeStamp.h>
#include <RAT/DS/MCParticle.hh>
#include <RAT/DS/MCTrack.hh>
#include <RAT/DS/MCSummary.hh>
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/MCHit.hh>

namespace RAT {
  namespace DS {

class MC : public TObject {
public:
  MC() : TObject() {}
  virtual ~MC() {}

  /**
   * Functor to compare an MCParticle name to a given string.
   *
   * Used for filtering the MCTrack list by particle type.
   */
  class ParticleNameEqualityFunctor {
  public:
    ParticleNameEqualityFunctor(std::string _name) : name(_name) {}
    bool operator()(MCTrack& track) {
      return track.GetParticleName() == this->name;
    }
  private:
    std::string name;
  };
  
  /** Event number. */
  virtual int GetID() const { return id; }
  virtual void SetID(int _id) { id = _id; }
  
  /** Absolute time of event. */
  virtual TTimeStamp GetUTC() const { return utc; }
  virtual void SetUTC(const TTimeStamp& _utc) { utc = _utc; }
      
  /** Initial particles in event */
  virtual MCParticle* GetMCParticle(Int_t i) { return &particle[i]; }
  virtual int GetMCParticleCount() const { return particle.size(); }
  virtual MCParticle* AddNewMCParticle() {
    particle.resize(particle.size() + 1);
    return &particle.back();
  }
  virtual void PruneMCParticle() { particle.resize(0); }

  /** Parent particles of interaction */
  virtual MCParticle* GetMCParent(Int_t i) { return &parent[i]; }
  virtual int GetMCParentCount() const { return parent.size(); }
  virtual MCParticle* AddNewMCParent() {
    parent.resize(parent.size() + 1);
    return &parent.back();
  }
  virtual void PruneMCParent() { parent.resize(0); }

  /** Particle tracks in detector.
   *
   *  Not filled unless the command

@verbatim 
/tracking/storeTrajectory 1 
@endverbatim

   *  is run in the user macro.
   */
  virtual MCTrack* GetMCTrack(int i) { return &track[i]; }
  virtual int GetMCTrackCount() const { return track.size(); }
  virtual MCTrack* AddNewMCTrack() {
    track.resize(track.size() + 1);
    return &track.back();
  }
  virtual void PruneMCTrack() { track.resize(0); };
  virtual void PruneMCTrack(const std::string& particleName) {
    ParticleNameEqualityFunctor pnef(particleName);
    std::remove_if(track.begin(), track.end(), pnef);
  }

  /** List of PMTs which had at least one photoelectron generated */
  virtual MCPMT* GetMCPMT(int i) { return &pmt[i]; }
  virtual int GetMCPMTCount() const { return pmt.size(); }
  virtual MCPMT* AddNewMCPMT() {
    pmt.resize(pmt.size() + 1);
    return &pmt.back();
  };
  virtual void PrunePMT() { pmt.resize(0); }

  /** List of hit channels */
  virtual MCHit* GetMCHit(int i) { return &hit[i]; }
  virtual int GetMCHitCount() const { return hit.size(); }
  virtual MCHit* AddNewMCHit() {
    hit.resize(hit.size() + 1);
    return &hit.back();
  };
  virtual void PruneMCHit() { hit.resize(0); }

  // Sort by hit time
  typedef bool (*MCHitSortFcn)(DS::MCHit a, DS::MCHit b);
  virtual inline void SortMCHit(MCHitSortFcn fcn) {
    std::sort(hit.begin(), hit.end(), fcn);
  }

  /** Total number of photoelectrons generated in this event */
  virtual int GetNumPE() const { return numPE; }
  virtual void SetNumPE(int _numPE) { numPE = _numPE; }

  /** Total number of "dark" photoelectrons generated in this event */
  virtual int GetNumDark() const { return numDarkHits; }
  virtual void SetNumDark(int _numDarkHits) { numDarkHits = _numDarkHits; }

  /** Summary of detector level properties of this event */
  virtual MCSummary* GetMCSummary() {
    if (summary.empty()) {
      summary.resize(1);
    }
    return &summary[0];
  }
  virtual bool ExistMCSummary() const { return !summary.empty(); }
  virtual void PruneMCSummary() { summary.resize(0); }

  ClassDef(MC, 1)
    
protected:
  int id;  //< ID number for this physics event
  int numPE;  //< Number of detected photoelectrons
  int numDarkHits;  //< Number of "dark" noise hits
  TTimeStamp utc;  //< Time of the event (UTC)
  std::vector<MCSummary> summary;  //< Collected summary data
  std::vector<MCParticle> particle;  //< List of primary particles
  std::vector<MCTrack> track;  //< List of tracks, if enabled
  std::vector<MCPMT> pmt;  //< List of MC truth PMTs
  std::vector<MCHit> hit;  //< List of MC truth electronics hits
  std::vector<MCParticle> parent;
};

  } // namespace DS
} // namespace RAT

#endif

