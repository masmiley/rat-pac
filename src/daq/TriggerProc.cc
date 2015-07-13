#include <algorithm>
#include <string>
#include <iostream>
#include <cfloat>
#include <vector>
#include <RAT/DB.hh>
#include <RAT/Log.hh>
#include <RAT/TimeUtil.hh>
#include <RAT/TriggerPulse.hh>
#include <RAT/TriggerSum.hh>
#include <RAT/DS/MCHit.hh>
#include <RAT/DS/MCSample.hh>
#include <RAT/DS/EV.hh>
#include <RAT/TriggerProc.hh>

namespace RAT {

inline bool Cmp_TPulse_TimeAscending(const TriggerPulse* a, 
                                     const TriggerPulse* b) {
  return a->GetStartTime() < b->GetStartTime();
}


TriggerProc::TriggerProc() : Processor("trigger") {
  fEventCount = 0;
  DBLinkPtr ldaq = DB::Get()->GetLink("DAQ");

  fTrigLockout = ldaq->GetD("lockout");
  fNHITPulseWidth = ldaq->GetD("nhit_width");
  fNThreshold = ldaq->GetD("nhit_thresh");

  info << dformat(" NHIT Threshold       : %3.1f hits\n", fNThreshold)
       << dformat(" Trigger Lockout Time : %3.1f ns\n", fTrigLockout)
       << dformat(" NHIT Pulse Width     : %5.2f ns\n", fNHITPulseWidth);
}


TriggerProc::Result TriggerProc::DSEvent(DS::Root* ds) {
  Log::Assert(ds->ExistMC(), "TriggerProc: No MC information found");

  // We need to build up arrays of TriggerPulse objects, sort them, and
  // then sum them to decide if we have a trigger.  
  TriggerSum nsum;

  int hitCount = ds->GetMC()->GetMCHitCount();
  for (int ihit=0; ihit<hitCount; ihit++) {
    DS::MCHit* hit = ds->GetMC()->GetMCHit(ihit);
    int sampleCount = hit->GetMCSampleCount();
    double lastTime = -DBL_MAX;
   
    // Now loop over samples, and create a trigger pulse each time there is
    // a new sample more than PulseWidthDB from the last (ie, allow retriggers)
    for (int isample=0; isample<sampleCount; isample++) {
      DS::MCSample* sample = hit->GetMCSample(isample);
      double time = sample->GetHitTime();
      if (time > lastTime + fNHITPulseWidth) {
        // Nhit trigger pulses
        SquareTriggerPulse* npulse = new SquareTriggerPulse();
        npulse->SetStartTime(sample->GetHitTime());
        npulse->SetHeight(1.0);
        npulse->SetWidth(fNHITPulseWidth);

        nsum.pulse.push_back(npulse);
        lastTime = sample->GetHitTime();
      }
    }
  }

  // Now we need to loop through the created trigger pulses, checking whether
  // the sum exceeds threshold. If so, create an event and fill the relevant
  // trigger information.
  //
  // We start by sorting the trigger pulses in time, which will make things 
  // faster.  We only look for a trigger if there have been trigger pulses
  // created.
  if (!nsum.pulse.empty()) {
    std::sort(nsum.pulse.begin(), nsum.pulse.end(), Cmp_TPulse_TimeAscending);
    double lastTime = -DBL_MAX;

    // Keep track of last hit included in an event
    lastTime = nsum.pulse[0]->GetStartTime();
    for (size_t i=0; i<nsum.pulse.size(); i++) {
      double sorttime = nsum.pulse[i]->GetStartTime();
      if (sorttime > lastTime) {
        float sum = nsum.GetHeight(i); // Gives height for sum up to pulse i
        if (sum > fNThreshold) {
          DS::EV* ev = ds->AddNewEV();
          ev->SetID(fEventCount);
          ev->SetUTC(AddNanoseconds(ds->GetMC()->GetUTC(),
                                    nsum.pulse[i]->GetStartTime()));
          ev->SetEventTime(nsum.pulse[i]->GetStartTime());
         
          // Step to first hit after trigger lockout time
          lastTime = nsum.pulse[i]->GetStartTime() + fTrigLockout;

          fEventCount++;
        }
      }
    }
  }

  return OK;
}

}  // namespace RAT

