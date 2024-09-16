/**
 *  @copyright Copyright 2020 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file EventCategorizer.cpp
 */

#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>
#include <JPetMCHit/JPetMCHit.h>
#include "EventCategorizerTools.h"
#include "EventCategorizer.h"
#include "HitFinderTools.h"
#include <iostream>

using namespace jpet_options_tools;
using namespace std;

EventCategorizer::EventCategorizer(const char* name): JPetUserTask(name) {}

const double EventCategorizer::kUnknownEventType = 66666666;

EventCategorizer::~EventCategorizer() {}

bool EventCategorizer::init()
{
  INFO("Event categorization started.");
  // Parameter for back to back categorization
  if (isOptionSet(fParams.getOptions(), kBack2BackSlotThetaDiffParamKey)){
    fB2BSlotThetaDiff = getOptionAsFloat(fParams.getOptions(), kBack2BackSlotThetaDiffParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kBack2BackSlotThetaDiffParamKey.c_str(), fB2BSlotThetaDiff
    ));
  }
  // Parameter for scattering determination
  if (isOptionSet(fParams.getOptions(), kScatterTOFTimeDiffParamKey)) {
    fScatterTOFTimeDiff = getOptionAsFloat(fParams.getOptions(), kScatterTOFTimeDiffParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kScatterTOFTimeDiffParamKey.c_str(), fScatterTOFTimeDiff
    ));
  }
  // Parameters for deexcitation TOT cut
  if (isOptionSet(fParams.getOptions(), kDeexTOTCutMinParamKey)) {
    fDeexTOTCutMin = getOptionAsFloat(fParams.getOptions(), kDeexTOTCutMinParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kDeexTOTCutMinParamKey.c_str(), fDeexTOTCutMin
    ));
  }
  if (isOptionSet(fParams.getOptions(), kDeexTOTCutMaxParamKey)) {
    fDeexTOTCutMax = getOptionAsFloat(fParams.getOptions(), kDeexTOTCutMaxParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kDeexTOTCutMaxParamKey.c_str(), fDeexTOTCutMax
    ));
  }
  if (isOptionSet(fParams.getOptions(), kMaxTimeDiffParamKey)) {
    fMaxTimeDiff = getOptionAsFloat(fParams.getOptions(), kMaxTimeDiffParamKey);
  } else {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kMaxTimeDiffParamKey.c_str(), fMaxTimeDiff));
  }
  // Getting bool for saving histograms
  if (isOptionSet(fParams.getOptions(), kSaveControlHistosParamKey)) {
    fSaveControlHistos = getOptionAsBool(fParams.getOptions(), kSaveControlHistosParamKey);
  }
  if (isOptionSet(fParams.getOptions(), kTOTCalculationType)) {
    fTOTCalculationType = getOptionAsString(fParams.getOptions(), kTOTCalculationType);
  } else {
    WARNING("No TOT calculation option given by the user. Using standard sum.");
  }


  // Input events type
  fOutputEvents = new JPetTimeWindow("JPetEvent");
  // Initialise hisotgrams
  if(fSaveControlHistos) initialiseHistograms();
  return true;
}

bool EventCategorizer::exec()
{
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {
    auto timeWindowMC = dynamic_cast<const JPetTimeWindowMC* const>(fEvent);
    vector<JPetEvent> events;
    for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++) {
      const auto& event_temp = dynamic_cast<const JPetEvent&>(timeWindow->operator[](i));
      auto hits = event_temp.getHits();

      hits.erase(std::remove_if(hits.begin(), hits.end(), [](auto hit) { return fabs(hit.getPosZ()) > 23; }), hits.end()); //remove hits with |z| > 23 cm

      JPetEvent event;
      event.setHits(hits);
      if(event.getHits().size()==0) continue;
      // Check types of current event

      fEnergy.clear();
      fTime.clear();
      fPosX.clear();
      fPosY.clear();
      fPosZ.clear();
      fHitType.clear();
      fVtxIndex.clear();

      for(const auto& hit : event.getHits()){
        fPosX.push_back(hit.getPosX());
        fPosY.push_back(hit.getPosY());
        fPosZ.push_back(hit.getPosZ());
        //fEnergy.push_back(hit.getEnergy());
        fTime.push_back(hit.getTime());
        if (timeWindowMC) {
          auto mcHit = timeWindowMC->getMCHit<JPetMCHit>(hit.getMCindex());
          fEnergy.push_back(mcHit.getEnergy());
          fHitType.push_back(mcHit.getGenGammaMultiplicity());
          fVtxIndex.push_back(mcHit.getMCVtxIndex());
        } else {
          fHitType.push_back(kUnknownEventType);
          fVtxIndex.push_back(0);
        }
      }
      if(fPosX.size() != event.getHits().size()) continue;

      bool isAccidental = EventCategorizerTools::checkForAccidental(
        event, getStatistics(), fSaveControlHistos, fVtxIndex
      );

      int ind = 0;
      vector<int> fScatterIndex={};
      vector<int> fPromptIndex={};
      vector<int> fSignalIndex={};
      
      for(auto i : fHitType){
        auto number = i%10;
        auto n = 1;
        if(i%100!=0 && number!=0 && i>10){
          fIsScattered=kTRUE;
          fScatterIndex.push_back(ind);}
        if(number==1){
          fContainsPrompt=kTRUE;
          fPromptIndex.push_back(ind);}
        else if(number==2)
          fIsPickOff=kTRUE;
        else if(number==3){
          fIsOPs=kTRUE;
          fSignalIndex.push_back(ind);}
        else if(number==0){
          fIsSecondary=kTRUE;
          while(number==0&&i>0){
            n = n*10;
            number = (i/n)%10;
            if(i%n*100 == 0 && number!=0 && i/n>10){
              fIsScattered=kTRUE;
              if(std::count(fScatterIndex.begin(), fScatterIndex.end(), ind)==0) fScatterIndex.push_back(ind);}
            if(number==1)
              fContainsPrompt=kTRUE;
            else if(number==2)
              fIsPickOff=kTRUE;
            else if(number==3){
              fIsOPs=kTRUE;
              fSignalIndex.push_back(ind);}
          }
	}
      ind++;
      }
      vector<bool> fHitInfo = {fContainsPrompt, fIsOPs, fIsPickOff, fIsScattered, fIsSecondary};

      bool is2Gamma = EventCategorizerTools::checkFor2Gamma(
        event, getStatistics(), fSaveControlHistos, fB2BSlotThetaDiff, fMaxTimeDiff
      );
      bool isPrompt = EventCategorizerTools::checkForPrompt(
        event, getStatistics(), fSaveControlHistos, fDeexTOTCutMin, fDeexTOTCutMax, fTOTCalculationType
      );
      bool isScattered = EventCategorizerTools::checkForScatter(
        event, getStatistics(), fSaveControlHistos, fScatterTOFTimeDiff, fTOTCalculationType, fHitInfo
      );
      bool isScattered_minimalDifference = EventCategorizerTools::checkForScatter_minimalDifference(
        event, getStatistics(), fSaveControlHistos, fScatterTOFTimeDiff, fTOTCalculationType, fHitInfo
      );
      bool is3Gamma = 0;
      if(isPrompt&&!isScattered){
        is3Gamma = EventCategorizerTools::checkFor3Gamma(
          event, getStatistics(), fSaveControlHistos
        );
       }
      JPetEvent newEvent = event;
      if(is2Gamma) newEvent.addEventType(JPetEventType::k2Gamma);
      if(is3Gamma) newEvent.addEventType(JPetEventType::k3Gamma);
      if(isPrompt) newEvent.addEventType(JPetEventType::kPrompt);
      if(isScattered) newEvent.addEventType(JPetEventType::kScattered);

      if(fSaveControlHistos){
      if(!isAccidental&&fIsOPs&&fContainsPrompt&&!fIsScattered&&!fIsSecondary){
      double average_time = 0;
	    for(auto i:fSignalIndex){
            average_time += fTime.at(i);
            getStatistics().fillHistogram("3Gamma_true_Energy", fEnergy.at(i));
            getStatistics().fillHistogram("3Gamma_true_Z_Energy", fEnergy.at(i), fPosZ.at(i));
            getStatistics().fillHistogram("3Gamma_true_XYpos", fPosX.at(i), fPosY.at(i));
          }
          if(fSignalIndex.size()==3 && fPromptIndex.size()==1){
            double lifetime = average_time/3. - fTime.at(fPromptIndex.at(0)); 
            getStatistics().fillHistogram("Lifetime_3Gamma_true", lifetime);
          }
	    }
       for(auto i:fScatterIndex){
         getStatistics().fillHistogram("Scattered_true_Energy", fEnergy.at(i));
         getStatistics().fillHistogram("Scattered_true_Z_Energy", fEnergy.at(i), fPosZ.at(i));
         getStatistics().fillHistogram("Scattered_true_XYpos", fPosX.at(i), fPosY.at(i));
       }
       for(auto i:fPromptIndex){
         if(!fIsScattered&&!fIsSecondary){
            getStatistics().fillHistogram("Prompt_true_Energy", fEnergy.at(i));
            getStatistics().fillHistogram("Prompt_true_Z_Energy", fEnergy.at(i), fPosZ.at(i));
            getStatistics().fillHistogram("Prompt_true_XYpos", fPosX.at(i), fPosY.at(i));
       }
       }
	
        for(auto hit : event.getHits()){
          getStatistics().fillHistogram("All_XYpos", hit.getPosX(), hit.getPosY());
	  getStatistics().fillHistogram("All_Energy", hit.getEnergy());
          getStatistics().fillHistogram("All_Z_Energy", hit.getEnergy(),  hit.getPosZ());
          getStatistics().fillHistogram("All_mult",event.getHits().size());	  
	}
      }

    fIsSecondary = kFALSE;
    fIsScattered = kFALSE;
    fContainsPrompt = kFALSE;
    fIsPickOff = kFALSE;
    fIsOPs = kFALSE;

    events.push_back(newEvent);
    }
    saveEvents(events);
  } else { return false; }
  return true;
}


bool EventCategorizer::terminate()
{
  INFO("Event categorization completed.");
  return true;
}

void EventCategorizer::saveEvents(const vector<JPetEvent>& events)
{
  for (const auto& event : events) { fOutputEvents->add<JPetEvent>(event); }
}

void EventCategorizer::initialiseHistograms(){

  // General histograms
  getStatistics().createHistogramWithAxes(
    new TH2D("All_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("All_Z_Energy", "Hit position Z vs Energy", 200, 0, 2e3, 201, -23.25, 23.25),
    "Energy [keV]", "Hit Z position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("All_Energy", "Energy", 200, 0, 2e3),
    "Energy [keV]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("All_mult", "Multiplicity", 11, -0.5, 10.5),
    "Multiplicity", "Number of Hits"
  );


  // Histograms for 2Gamma category
  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_Zpos", "Z-axis position of 2 gamma hits", 201, -23.25, 23.25),
    "Z axis position [cm]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_DLOR", "Delta LOR distance", 100, -0.25, 49.25),
    "Delta LOR [cm]", "Counts"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_ThetaDiff", "Angle difference of 2 gamma hits ", 181, -0.5, 180.5),
    "Hits theta diff [deg]", "Counts"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_TimeDiff", "Time difference of 2 gamma hits", 200, -10100.0, 99900.0),
    "Time Difference [ps]", "Number of Hit Pairs"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_Dist", "B2B hits distance", 150, -0.5, 149.5),
    "Distance [cm]", "Number of Hit Pairs"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Annih_TOF", "Annihilation pairs Time of Flight", 201, -3015.0, 3015.0),
    "Time of Flight [ps]", "Number of Annihilation Pairs"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("AnnihPoint_XY", "XY position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "X position [cm]", "Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("AnnihPoint_ZX", "ZX position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Z position [cm]", "X position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("AnnihPoint_ZY", "ZY position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Z position [cm]", "Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Annih_DLOR", "Delta LOR distance of annihilation photons", 100, -0.25, 49.25),
    "Delta LOR [cm]", "Counts"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_Energy", "Energy", 200, 0, 2e3),
    "Energy [keV]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("2Gamma_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("2Gamma_Z_Energy", "Hit position Z vs Energy", 200, 0, 2e3, 201, -23.25, 23.25),
    "Energy [keV]", "Hit Z position [cm]"
  );


  // Histograms for 3Gamama category
  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_Angles", "Relative angles - transformed", 250, -0.5, 249.5, 20, -0.5, 199.5),
    "Relative angle 1-2", "Relative angle 2-3"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_Angles_after_cut", "Relative angles - transformed after cut", 250, -0.5, 249.5, 200, -0.5, 199.5),
    "Relative angle 1-2", "Relative angle 2-3"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_Angles3D", "Relative angles - transformed", 250, -0.5, 249.5, 200, -0.5, 199.5),
    "Relative angle 1-2", "Relative angle 2-3"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_Angles3D_after_cut", "Relative angles - transformed after cut", 250, -0.5, 249.5, 200, -0.5, 199.5),
    "Relative angle 1-2", "Relative angle 2-3"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("3Gamma_Energy", "Energy", 200, 0, 2e3),
    "Energy [keV]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("3Gamma_Energy_vertex_cut", "Energy", 200, 0, 2e3),
    "Energy [keV]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("3Gamma_true_Energy", "Energy", 200, 0, 2e3),
    "Energy [keV]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_true_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_Z_Energy", "Hit position Z vs Energy", 200, 0, 2e3, 201, -23.25, 23.25),
    "Energy [keV]", "Hit Z position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_true_Z_Energy", "Hit position Z vs Energy", 200, 0, 2e3, 201, -23.25, 23.25),
    "Energy [keV]", "Hit Z position [cm]"
  );
  
  getStatistics().createHistogramWithAxes(
    new TH1D("Lifetime_3Gamma_true", "Lifetime of 3 gamma events",
    200, -1700000, 1700000),
    "Lifetime [ps]", "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Lifetime_3Gamma", "Lifetime of 3 gamma events",
    200, -1700000, 1700000),
    "Lifetime [ps]", "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Lifetime_3Gamma_vertex_cut", "Lifetime of 3 gamma events",
    200, -1700000, 1700000),
    "Lifetime [ps]", "Number of Events"
  );


  getStatistics().createHistogramWithAxes(
    new TH2D("decay point XY", "transverse position of the o-Ps->3g decay point", 100, -50., 50., 100, -50., 50.),
    "X [cm]", "Y [cm]"
	);

  getStatistics().createHistogramWithAxes(
    new TH2D("decay point XZ", "position of the o-Ps->3g decay point in XZ", 100, -50., 50.,100, -50., 50.),
    "Z [cm]", "X [cm]"
	);

getStatistics().createHistogramWithAxes(
    new TH1D("decay point time", "time of the o-Ps->3g decay point", 250, -2500., 2500.),
    "t [ps]", "Number of events"
	);

getStatistics().createHistogramWithAxes(
    new TH1D("decay point time within detector", "time of the o-Ps->3g decay point", 500, -50000., 50000.),
    "t [ps]", "Number of events"
        );

  getStatistics().createHistogramWithAxes(
    new TH1D("decay point time within detector flag", "time of the o-Ps->3g decay point", 500, -50000., 50000.),
    "t [ps]", "Number of events"
        );

  getStatistics().createHistogramWithAxes(
    new TH2D("decay point XY all", "transverse position of the o-Ps->3g decay point", 100, -50., 50., 100, -50., 50.),
    "X [cm]", "Y [cm]"
        );

  getStatistics().createHistogramWithAxes(
    new TH2D("decay point XZ all", "position of the o-Ps->3g decay point in XZ", 100, -50., 50.,100, -50., 50.),
    "Z [cm]", "X [cm]"
        );

  getStatistics().createHistogramWithAxes(
    new TH2D("decay point XY within detector", "transverse position of the o-Ps->3g decay point", 110, -55., 55., 110, -55., 55.),
    "X [cm]", "Y [cm]"
        );

  getStatistics().createHistogramWithAxes(
    new TH2D("decay point XZ within detector", "position of the o-Ps->3g decay point in XZ", 100, -50., 50.,110, -55., 55.),
    "Z [cm]", "X [cm]"
        );


  getStatistics().createHistogramWithAxes(
    new TH2D("decay point XY within detector flag", "transverse position of the o-Ps->3g decay point", 110, -55., 55., 110, -55., 55.),
    "X [cm]", "Y [cm]"
        );

  getStatistics().createHistogramWithAxes(
    new TH2D("decay point XZ within detector flag", "position of the o-Ps->3g decay point in XZ", 100, -50., 50.,110, -55., 55.),
    "Z [cm]", "X [cm]"
        );

getStatistics().createHistogramWithAxes(
    new TH1D("decay point time all", "time of the o-Ps->3g decay point", 5000, -50000., 50000.),
    "t [ps]", "Number of events"
	);

getStatistics().createHistogramWithAxes(
    new TH1D("decay point error", "error of the o-Ps->3g decay point", 5, -0.5, 4.5),
    "error", "Number of events"
	);

  
getStatistics().createHistogramWithAxes(
    new TH1D("decay point error within detector", "error of the o-Ps->3g decay point", 5, -0.5, 4.5),
    "error", "Number of events"
	);

getStatistics().createHistogramWithAxes(
    new TH1D("decay point error within detector flag", "error of the o-Ps->3g decay point", 5, -0.5, 4.5),
    "error", "Number of events"
	);


  // Histograms for scattering category
  getStatistics().createHistogramWithAxes(
    new TH1D("ScatterTOF_TimeDiff", "Difference of Scatter TOF and hits time difference",
    3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Scat_TOF - time diff [ps]", "Number of Hit Pairs"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("smallest_difference_ScatterTOF_TimeDiff", "Difference of Scatter TOF and hits time difference",
    3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Scat_TOF - time diff [ps]", "Number of Hit Pairs"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("ScatterTOF_TimeDiff_before", "Difference of Scatter TOF and hits time difference",
    3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Scat_TOF - time diff [ps]", "Number of Hit Pairs"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_TOF_vs_TimeDiff", "Scatter TOF vs hits time difference",
    fScatterTOFTimeDiff/10., -500, 3.0*fScatterTOFTimeDiff-0.5,
    fScatterTOFTimeDiff/10., -500, 3.0*fScatterTOFTimeDiff-0.5),
    "Time diff [ps]", "Scat_TOF [ps]"
  );

    getStatistics().createHistogramWithAxes(
    new TH2D("smallest_difference_Scatt_TOF_vs_TimeDiff", "Scatter TOF vs hits time difference",
    fScatterTOFTimeDiff/10., -500, 3.0*fScatterTOFTimeDiff-0.5,
    fScatterTOFTimeDiff/10., -500, 3.0*fScatterTOFTimeDiff-0.5),
    "Time diff [ps]", "Scat_TOF [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_TOF_vs_TimeDiff_before", "Scatter TOF vs hits time difference",
    fScatterTOFTimeDiff/10., -500, 3.0*fScatterTOFTimeDiff-0.5,
    fScatterTOFTimeDiff/10., -500, 3.0*fScatterTOFTimeDiff-0.5),
    "Time diff [ps]", "Scat_TOF [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_TOF_vs_TimeDiff_before_scatter", "Scatter TOF vs hits time difference",
    fScatterTOFTimeDiff/10., -500, 3.0*fScatterTOFTimeDiff-0.5,
    fScatterTOFTimeDiff/10., -500, 3.0*fScatterTOFTimeDiff-0.5),
    "Time diff [ps]", "Scat_TOF [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_TOF_vs_TimeDiff_before_prompt", "Scatter TOF vs hits time difference",
    fScatterTOFTimeDiff/10., -500, 3.0*fScatterTOFTimeDiff-0.5,
    fScatterTOFTimeDiff/10., -500, 3.0*fScatterTOFTimeDiff-0.5),
    "Time diff [ps]", "Scat_TOF [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_TOF_vs_TimeDiff_before_signal", "Scatter TOF vs hits time difference",
    fScatterTOFTimeDiff/10., -500, 3.0*fScatterTOFTimeDiff-0.5,
    fScatterTOFTimeDiff/10., -500, 3.0*fScatterTOFTimeDiff-0.5),
    "Time diff [ps]", "Scat_TOF [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_Dist_vs_M", "Difference of Scatter TOF and hits time difference vs distance",
    120, 0, 120, 3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Distance [cm]", "Scat_TOF - time diff [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("smallest_difference_Scatt_Dist_vs_M", "Difference of Scatter TOF and hits time difference vs distance",
    120, 0, 120, 3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Distance [cm]", "Scat_TOF - time diff [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_Dist_vs_M_before", "Difference of Scatter TOF and hits time difference vs distance",
    120, 0, 120, 3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Distance [cm]", "Scat_TOF - time diff [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_Dist_vs_M_before_signal", "Difference of Scatter TOF and hits time difference vs distance",
    120, 0, 120, 3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Distance [cm]", "Scat_TOF - time diff [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_Dist_vs_M_before_prompt", "Difference of Scatter TOF and hits time difference vs distance",
    120, 0, 120, 3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Distance [cm]", "Scat_TOF - time diff [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_Dist_vs_M_before_scatter", "Difference of Scatter TOF and hits time difference vs distance",
    120, 0, 120, 3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Distance [cm]", "Scat_TOF - time diff [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_Dist_2D_vs_M", "Difference of Scatter TOF and hits time difference vs distance",
    50, 0, 50, 3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Distance [cm]", "Scat_TOF - time diff [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_Dist_2D_vs_M_before", "Difference of Scatter TOF and hits time difference vs distance",
    50, 0, 50, 3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Distance [cm]", "Scat_TOF - time diff [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_Dist_2D_vs_M_before_scatter", "Difference of Scatter TOF and hits time difference vs distance",
    50, 0, 50, 3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Distance [cm]", "Scat_TOF - time diff [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("ScatterAngle_PrimaryTOT", "Angle of scattering vs. TOT of primary hits",
    200, -0.5, 199.5, 200, -100.0, 39900.0),
    "Scattering Angle", "TOT of primary hit [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("smallest_difference_ScatterAngle_PrimaryTOT", "Angle of scattering vs. TOT of primary hits",
    200, -0.5, 199.5, 200, -100.0, 39900.0),
    "Scattering Angle", "TOT of primary hit [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("ScatterAngle_PrimaryTOT_before", "Angle of scattering vs. TOT of primary hits",
    200, -0.5, 199.5, 200, -100.0, 39900.0),
    "Scattering Angle", "TOT of primary hit [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("ScatterAngle_ScatterTOT", "Angle of scattering vs. TOT of scattered hits",
    200, -0.5, 199.5, 200, -100.0, 39900.0),
    "Scattering Angle", "TOT of scattered hit [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("smallest_difference_ScatterAngle_ScatterTOT", "Angle of scattering vs. TOT of scattered hits",
    200, -0.5, 199.5, 200, -100.0, 39900.0),
    "Scattering Angle", "TOT of scattered hit [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("ScatterAngle_ScatterTOT_before", "Angle of scattering vs. TOT of scattered hits",
    200, -0.5, 199.5, 200, -100.0, 39900.0),
    "Scattering Angle", "TOT of scattered hit [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Scattered_Energy", "Energy", 200, 0, 2e3),
    "Energy [keV]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("smallest_difference_Scattered_Energy", "Energy", 200, 0, 2e3),
    "Energy [keV]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Scattered_true_Energy", "Energy", 200, 0, 2e3),
    "Energy [keV]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scattered_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scattered_true_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("smallest_difference_Scattered_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scattered_Z_Energy", "Hit position Z vs Energy", 200, 0, 2e3, 201, -23.25, 23.25),
    "Energy [keV]", "Hit Z position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("smallest_difference_Scattered_Z_Energy", "Hit position Z vs Energy", 200, 0, 2e3, 201, -23.25, 23.25),
    "Energy [keV]", "Hit Z position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scattered_true_Z_Energy", "Hit position Z vs Energy", 200, 0, 2e3, 201, -23.25, 23.25),
    "Energy [keV]", "Hit Z position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Hit_multiplicity_scatt0", "Multiplicity", 11, -0.5, 10.5),
    "Multiplicity", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Hit_multiplicity_scatt1", "Multiplicity", 11, -0.5, 10.5),
    "Multiplicity", "Number of Hits"
  );

  // Histograms for deexcitation
  getStatistics().createHistogramWithAxes(
    new TH1D("Deex_TOT_cut", "TOT of all hits with deex cut (30,50) ns", 200, 24950.0, 54950.0),
    "Energy [keV]", "Number of Hits"
  );

  // Histograms for accidentals

  getStatistics().createHistogramWithAxes(
    new TH1D("Acc_Energy", "Energy", 200, 0, 2e3),
    "Energy [keV]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Acc_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Acc_Z_Energy", "Hit position Z vs Energy", 200, 0, 2e3, 201, -23.25, 23.25),
    "Energy [keV]", "Hit Z position [cm]"
  );

  // Histograms for prompts

  getStatistics().createHistogramWithAxes(
    new TH1D("Prompt_true_Energy", "Energy", 200, 0, 2e3),
    "Energy [keV]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Prompt_true_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Prompt_true_Z_Energy", "Hit position Z vs Energy", 200, 0, 2e3, 201, -23.25, 23.25),
    "Energy [keV]", "Hit Z position [cm]"
  );


}
