#include <cmath>
#include <algorithm>
#include <vector>
#include <CLHEP/Random/RandGauss.h>
#include <RAT/Log.hh>
#include <RAT/DB.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/MCSample.hh>
#include <RAT/DS/MCHit.hh>
#include <RAT/FrontEndProc.hh>

namespace RAT {

FrontEndProc::FrontEndProc() : Processor("frontend") {
  DBLinkPtr ldaq = DB::Get()->GetLink("DAQ");

  fSamplingTime = ldaq->GetD("sampling_time"); 
  fHitThreshold = ldaq->GetD("chan_thresh");
  fPulseWidth = ldaq->GetD("pulse_width");
  fNoiseAmpl = ldaq->GetD("noise_amplitude"); 
  fGDelay = ldaq->GetD("gate_delay"); 

  info << dformat(" \n \n")
       << dformat(" DAQ Constants Loaded:\n")
       << dformat(" PMT Pulse Width              : %3.1f ns\n", fPulseWidth)
       << dformat(" PMT Channel Integration Time : %5.2f ns\n", fSamplingTime)
       << dformat(" Channel Gate Delay           : %3.1f ns\n", fGDelay)
       << dformat(" PMT Channel Hit Threshold    : %3.2f pe\n", fHitThreshold)
       << dformat(" Hi Freq. Channel Noise       : %3.2f pe\n", fNoiseAmpl);
}


Processor::Result FrontEndProc::DSEvent(DS::Root* ds) {
  Log::Assert(ds->ExistMC(), "FrontEndProc: No MC information found.");
  bool isHit = false;
  // Outer loop is over all hit MCPMTs
  int pmtCount = ds->GetMC()->GetMCPMTCount();
  for (int ipmt=0; ipmt<pmtCount; ipmt++) {
    DS::MCPMT* pmt = ds->GetMC()->GetMCPMT(ipmt);
    DS::MCHit* hit = NULL;

    // For each hit PMT, loop over the photoelectrons. At the time
    // of each photoelectron, we need to add up all the pulseheights that 
    // contribute to the charge sum seen by the discriminator
    int photonCount = pmt->GetMCPhotonCount();
    for (int iphoton=0; iphoton<photonCount; iphoton++) {
      DS::MCPhoton* photon = pmt->GetMCPhoton(iphoton);
      double timeNow = photon->GetHitTime();
     if (!isHit) {
      info << dformat("PMT count: %d Photon count: %d PMT: %d Photon: %d \n", pmtCount, photonCount, ipmt, iphoton) 
           << dformat("timeNow: %f \n", timeNow);
      }
      // Now we check how many pulses before our current pulse contribute
      // to the total seen by the front end discriminator
      int sumIndex = iphoton;
      double deltaT = 0;
      double charge = photon->GetCharge();
     // if (!isHit) {
     //  info  << dformat("Before loop sumIndex : %d  \n", sumIndex) << "Charge: " << charge << " Photon get charge: " << photon->GetCharge() << "\n";
     // }
      while (deltaT <= fPulseWidth && sumIndex > 0) {
         sumIndex--;
         charge += pmt->GetMCPhoton(sumIndex)->GetCharge();
         deltaT = timeNow - pmt->GetMCPhoton(sumIndex)->GetHitTime();
      }
     //if (!isHit) {
     //   info << dformat("After loop  sumIndex : %d \n", sumIndex)
     //        << dformat("Charge: %f PE", charge) << dformat("     deltaT: %f ns\n", deltaT);
     // }

      // Now add some noise to the sum, and check it against threshold
      float qfuzzed = charge + fNoiseAmpl * CLHEP::RandGauss::shoot();
      if (qfuzzed >= fHitThreshold) {
        // We have a hit above threshold
        if (!hit) {
          // Fill hit info if PMT not already hit
          hit = ds->GetMC()->AddNewMCHit();
          hit->SetPMTID(pmt->GetID());
          hit->SetCrate(0);  // FIXME Made up for now
          hit->SetCard(0);
          hit->SetChannel(0);
        }

        // Now need to collect photoelectrons which belong to this sample
        // so step back in time a little, and then go forward a full 
        // sampling time.
        DS::MCSample* sample = hit->AddNewMCSample();

        // Time of photon which crossed threshold
        sample->SetHitTime(timeNow);
        if (ds->GetMC()->GetMCHitCount() == 1) {
          info << "Sample time: " <<  sample->GetHitTime()<< "\n";
          isHit = true;
        }
        double timeStartGate = timeNow - fGDelay;
        //if (photonCount > 1) {
        //   info << "Time Now: " << timeNow << " Gate start: " << timeStartGate << "Difference: " << timeNow - timeStartGate << " \n";
        //}
        // Loop backwards until we're at the beginning of the gate
        int hitIndex = iphoton;
        double pmttime = timeNow;
        if (hitIndex > 0) {
          while (pmttime >= timeStartGate && hitIndex > 0) {
            hitIndex--;
            pmttime = pmt->GetMCPhoton(hitIndex)->GetHitTime();
          }

          // One step forward so we're inside gate
          hitIndex++;
        }

        pmttime = pmt->GetMCPhoton(hitIndex)->GetHitTime();
        sample->SetCharge(0);
       // if (photonCount > 1) {
        // info << "pmttime: " << pmttime << " timeStartGate + fSamplingTime: " << timeStartGate + fSamplingTime << "\n hitIndex: " << hitIndex
        //      << " Photon count " << pmt->GetMCPhotonCount() - 1 << " \n";
       // }
        // Now loop forward until the integration gate closes
        while (pmttime <= timeStartGate + fSamplingTime &&
               hitIndex < pmt->GetMCPhotonCount()) {// - 1) {
          sample->SetCharge(sample->GetCharge() +
                            pmt->GetMCPhoton(hitIndex)->GetCharge()); 
          hitIndex++;
          if (hitIndex < pmt->GetMCPhotonCount()) {
           pmttime = pmt->GetMCPhoton(hitIndex)->GetHitTime();
          }
        }

        // Skip to the last hit sampled
        iphoton = hitIndex;
      }
    }
  }

  return OK;
}

}  // namespace RAT

