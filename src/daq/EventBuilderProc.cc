#include <algorithm>
#include <cmath>
#include <vector>
#include <RAT/DB.hh>
#include <RAT/Log.hh>
#include <RAT/TimeUtil.hh>
#include <RAT/TriggerPulse.hh>
#include <RAT/TriggerSum.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/MCHit.hh>
#include <RAT/DS/MCSample.hh>
#include <RAT/DS/MCHit.hh>
#include <G4ThreeVector.hh>
#include <RAT/EventBuilderProc.hh>

namespace RAT {

inline bool Cmp_MCHit_TimeAscending(DS::MCHit a, DS::MCHit b) {
  return a.GetMCSample(0)->GetHitTime() < b.GetMCSample(0)->GetHitTime();
}


EventBuilderProc::EventBuilderProc() : Processor("eventbuilder") {
  DBLinkPtr ldaq = DB::Get()->GetLink("DAQ");
  
  fTriggerDelay = ldaq->GetD("trigger_delay");
  fTriggerGate = ldaq->GetD("trigger_gate");

  info << dformat(" Trigger Delay : %5.2f ns\n", fTriggerDelay)
       << dformat(" Trigger Gate  : %5.2f ns\n", fTriggerGate);
}


Processor::Result EventBuilderProc::DSEvent(DS::Root* ds) {
  Log::Assert(ds->ExistMC(), "EventBuilderProc: No MC information found.");
 
 // First sort the MCHits
  ds->GetMC()->SortMCHit(Cmp_MCHit_TimeAscending);
  // for (int ihit = 0; ihit < ds->GetMC()->GetMCHitCount(); ihit++) {
  //  for (int isample = 0; isample < ds->GetMC()->GetMCHit(ihit)->GetMCSampleCount(); isample++) {
  //    info << "Hit: " << ihit << " Sample: " << isample << " Sample time: " 
  //         << ds->GetMC()->GetMCHit(ihit)->GetMCSample(isample)->GetHitTime() << "\n";
  //  }
  //}
  // Now loop over all triggers, collecting hits and samples
  for (int iev=0; iev<ds->GetEVCount(); iev++) {
    DS::EV* ev = ds->GetEV(iev);
    info <<"iev: " << iev <<  "Event time: " << ev->GetUTC().GetNanoSec() << "\n";

    //Set DeltaT for ev
    if (iev == 0) {
      ev->SetDeltaT(0);
    }
    else {
      ev->SetDeltaT(TimeDifference(ev->GetUTC(), ds->GetEV(iev-1)->GetUTC()));
      //info << "EV time difference: " << TimeDifference(ev->GetUTC(), ds->GetEV(iev-1)->GetUTC()) << "\n";
    }

    ev->SetTotalCharge(0);
    double triggerStopTime = ev->GetEventTime() +  fTriggerDelay;
   
    if (ds->GetMC()->GetMCHitCount() == 0) {
      continue;
    }

    // Loop over hits. firstSampleTime tells us the time of the first sample
    // in a given MCHit.
    double firstSampleTime = ds->GetMC()->GetMCHit(0)->GetMCSample(0)->GetHitTime();
    info << "First sample time retrieved: " << firstSampleTime << " \n";
    int ihit = 0;
    info <<"Hit count: " << ds->GetMC()->GetMCHitCount() << " triggerStopTime: " <<  triggerStopTime << " Difference: " 
         << triggerStopTime -  firstSampleTime<< " \n Hits: " << ds->GetMC()->GetMCHitCount() << " \n \n";
    while (ihit < ds->GetMC()->GetMCHitCount() &&
           firstSampleTime < triggerStopTime) {
      DS::MCHit* hit = ds->GetMC()->GetMCHit(ihit);

      if (ihit % 20 == 0) {
        info << "Looping over hits. At: " << ihit << "\n";
      }
      // Advance through samples on this hit until we find the first one 
      // (if any) that is after the last trigger lockout time.
      double nextSampleTime;
      int imcsample = 0;
      while (imcsample < hit->GetMCSampleCount()) {
        DS::MCSample* mcsample = hit->GetMCSample(imcsample);

        // If sample is inside gate, create a PMT data object and fill the
        // time and charge samples, as well as PMT information..
        nextSampleTime = mcsample->GetHitTime();

        if ((triggerStopTime - nextSampleTime) < fTriggerGate &&
             nextSampleTime < triggerStopTime) {
          DS::PMT* pmt = ev->AddNewPMT();
          pmt->SetID(hit->GetPMTID());

          // Loop over the rest of the MC samples, filling PMT data samples
          // until we run out of the trigger gate.
          while (imcsample < hit->GetMCSampleCount() &&
                   triggerStopTime - hit->GetMCSample(imcsample)->GetHitTime() > 0) {
            DS::Sample* sample = pmt->AddNewSample();
            double mctime = mcsample->GetHitTime();
            double mccharge = mcsample->GetCharge();
            sample->SetTime(DigitizeT(mctime, triggerStopTime));
            sample->SetCharge(DigitizeQ(mccharge));
            ev->SetTotalCharge(ev->GetTotalCharge() + sample->GetCharge());
            imcsample++;
          }
        }
        else {
          imcsample++;
        }
      }

     // info << "Setting first sample time again";
      if (ihit < ds->GetMC()->GetMCHitCount()) {
        firstSampleTime = hit->GetMCSample(0)->GetHitTime();
      }
     // info << "Successful set time";
      ihit++;
    }
  }

  return OK;
}

// Eventually this will have to return a 12-bit integer,
// in ADC or TDC counts
double EventBuilderProc::DigitizeT(double hitTime, double stopTime) {
  return stopTime - hitTime;
}


// Eventually this will have to return a 12-bit integer,
// in ADC counts
float EventBuilderProc::DigitizeQ(float hitCharge) {
  return hitCharge;
}

}  // namespace RAT

