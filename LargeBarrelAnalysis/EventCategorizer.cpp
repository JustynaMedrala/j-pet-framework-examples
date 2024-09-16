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


  // Parameters for 2 annihilation TOT cut
  if (isOptionSet(fParams.getOptions(), kToTCut2AnniMinParamKey)) {
    fToTCut2AnniMin = getOptionAsFloat(fParams.getOptions(), kToTCut2AnniMinParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kToTCut2AnniMinParamKey.c_str(), fToTCut2AnniMin
    ));
  }
  if (isOptionSet(fParams.getOptions(), kToTCut2AnniMaxParamKey)) {
    fToTCut2AnniMax = getOptionAsFloat(fParams.getOptions(), kToTCut2AnniMaxParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kToTCut2AnniMaxParamKey.c_str(), fToTCut2AnniMax
    ));
  }

  // Parameters for 3 annihilation TOT cut
  if (isOptionSet(fParams.getOptions(), kToTCut3AnniMinParamKey)) {
    fToTCut3AnniMin = getOptionAsFloat(fParams.getOptions(), kToTCut3AnniMinParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kToTCut3AnniMinParamKey.c_str(), fToTCut3AnniMin
    ));
  }
  if (isOptionSet(fParams.getOptions(), kToTCut3AnniMaxParamKey)) {
    fToTCut3AnniMax = getOptionAsFloat(fParams.getOptions(), kToTCut3AnniMaxParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kToTCut2AnniMaxParamKey.c_str(), fToTCut3AnniMax
    ));
  }

  // Parameters for scatters TOT cut
  if (isOptionSet(fParams.getOptions(), kToTCutScattMaxParamKey)) {
    fToTCutScattMax = getOptionAsFloat(fParams.getOptions(), kToTCutScattMaxParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kToTCutScattMaxParamKey.c_str(), fToTCutScattMax
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

  if (isOptionSet(fParams.getOptions(), kHistoLimit)) {
    fHistoLimit = getOptionAsFloat(fParams.getOptions(), kHistoLimit);
  } else {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kHistoLimit.c_str(), fHistoLimit));
  }

  if (isOptionSet(fParams.getOptions(), kEnergy_unit)) {
    fEnergy_unit = getOptionAsString(fParams.getOptions(), kEnergy_unit);
  } else {
    WARNING("No historgam title.");
  }

  


  // Input events type
  fOutputEvents = new JPetTimeWindow("JPetEvent");
  // Initialise hisotgrams
  if(fSaveControlHistos) initialiseHistograms();
  return true;
}

bool EventCategorizer::exec()
{
  
  JPetTimeWindowMC *timeWindowMC = dynamic_cast<JPetTimeWindowMC *const>(fEvent);
  if (timeWindowMC) {
    fIsMC = true;
    INFO("The input file is MC.");
  } else {
    INFO("The input file is DATA.");
  }

  n_oPs_true = 0;
  n_oPs_prompt_true = 0;
  n_oPs_scatter = 0;
  n_oPs_prompt_scatter = 0;
  n_prompt_true = 0;
  n_random = 0;

  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {
    vector<JPetEvent> events;
    for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++) {
      auto& event_temp = dynamic_cast<const JPetEvent&>(timeWindow->operator[](i));
      auto hits = event_temp.getHits();

      hits.erase(std::remove_if(hits.begin(), hits.end(), [](auto hit) { return fabs(hit.getPosZ()) > 23; }), hits.end()); //remove hits with |z| > 23 cm
      
      JPetEvent event;
      event.setHits(hits);

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
        //cout<<"Hit type: ";
        if (fIsMC) {
          auto mcHit = timeWindowMC->getMCHit<JPetMCHit>(hit.getMCindex());
          fEnergy.push_back(mcHit.getEnergy());
          fHitType.push_back(mcHit.getGenGammaMultiplicity());
          fVtxIndex.push_back(mcHit.getMCVtxIndex());
          //cout<<mcHit.getGenGammaMultiplicity()<<", ";
        } else {
          fHitType.push_back(kUnknownEventType);
          fVtxIndex.push_back(0);
        }
        //cout<<endl;
      }
      if(fPosX.size() != event.getHits().size()) continue;
      if(fIsMC) {
          bool isAccidental = false;
          if(fVtxIndex.size() != 0) isAccidental = EventCategorizerTools::checkForAccidental(event, getStatistics(), 
                      fSaveControlHistos, fVtxIndex);
          vector<int> fScatterIndex={};
          vector<int> fPromptIndex={};
          vector<int> fSignalIndex={};
          EventCategorizerTools::HitTypeMC(fHitType, fIsScattered, fContainsPrompt, fIsPickOff, fIsOPs, 
                                           fIsSecondary,  fScatterIndex, fPromptIndex, fSignalIndex);
          vector<bool> fHitInfo = {fContainsPrompt, fIsOPs, fIsPickOff, fIsScattered, fIsSecondary};
        if(fSaveControlHistos){
          if(fIsOPs&&(fIsScattered||fIsSecondary)){
            n_oPs_scatter = 1;
          }
          if(isAccidental){
            n_random = 1;
          }
          if(fIsOPs&&fContainsPrompt&&(fIsScattered||fIsSecondary)){
            n_oPs_prompt_scatter = 1;
          }
          if(fContainsPrompt&&(fIsScattered||fIsSecondary)){
            n_prompt_scatter = 1;
          }
          if(fIsOPs&&fContainsPrompt&&!fIsScattered&&!fIsSecondary&!isAccidental){
          n_oPs_true = 1;
          double average_time = 0;
          for(auto i:fSignalIndex){
                average_time += fTime.at(i);
                getStatistics().fillHistogram("3Gamma_true_Energy", fEnergy.at(i));
                getStatistics().fillHistogram("3Gamma_true_Z_Energy", fEnergy.at(i), fPosZ.at(i));
                getStatistics().fillHistogram("3Gamma_true_XYpos", fPosX.at(i), fPosY.at(i));
              }
              if(fSignalIndex.size()==3 && fPromptIndex.size()==1){
                n_oPs_prompt_true = 1;
                double lifetime = average_time/3. - fTime.at(fPromptIndex.at(0)); 
                getStatistics().fillHistogram("Lifetime_3Gamma_true", lifetime);
              }
          }
          for(auto i:fScatterIndex){
            n_scattered = 1;
            getStatistics().fillHistogram("Scattered_true_Energy", fEnergy.at(i));
            getStatistics().fillHistogram("Scattered_true_Z_Energy", fEnergy.at(i), fPosZ.at(i));
            getStatistics().fillHistogram("Scattered_true_XYpos", fPosX.at(i), fPosY.at(i));
          }
          for(auto i:fPromptIndex){
            if(!fIsScattered&&!fIsSecondary){
                n_prompt_true = 1;
                getStatistics().fillHistogram("Prompt_true_Energy", fEnergy.at(i));
                getStatistics().fillHistogram("Prompt_true_Z_Energy", fEnergy.at(i), fPosZ.at(i));
                getStatistics().fillHistogram("Prompt_true_XYpos", fPosX.at(i), fPosY.at(i));
            }
          }
        } 

      }

      // Check types of current event
      bool is2Gamma = EventCategorizerTools::checkFor2Gamma(
        event, getStatistics(), fSaveControlHistos, fB2BSlotThetaDiff, fMaxTimeDiff, fDeexTOTCutMin, fDeexTOTCutMax, 
        fToTCut2AnniMin, fToTCut2AnniMax
      );
      bool is3Gamma = false;
      bool isPrompt = EventCategorizerTools::checkForPrompt(
        event, getStatistics(), fSaveControlHistos, fDeexTOTCutMin, fDeexTOTCutMax, fTOTCalculationType
      );
      bool isScattered = EventCategorizerTools::checkForScatter(
        event, getStatistics(), fSaveControlHistos, fScatterTOFTimeDiff, fTOTCalculationType, fToTCutScattMax
      );
      bool isScattered_minimalDifference = EventCategorizerTools::checkForScatter_minimalDifference(
        event, getStatistics(), fSaveControlHistos, fScatterTOFTimeDiff, fTOTCalculationType, fToTCutScattMax
      );
      if(isPrompt&&!isScattered){
        is3Gamma = EventCategorizerTools::checkFor3Gamma(
          event, getStatistics(), fSaveControlHistos,  fDeexTOTCutMin, fDeexTOTCutMax, fToTCut3AnniMin, fToTCut3AnniMax
        );
       }
      JPetEvent newEvent = event;
      if(is2Gamma) newEvent.addEventType(JPetEventType::k2Gamma);
      if(is3Gamma) newEvent.addEventType(JPetEventType::k3Gamma);
      if(isPrompt) newEvent.addEventType(JPetEventType::kPrompt);
      if(isScattered) newEvent.addEventType(JPetEventType::kScattered);

      if(fSaveControlHistos){
	      double energy = 0;
        double energy_syn = 0;
        for(auto hit : event.getHits()){
          getStatistics().fillHistogram("All_XYpos", hit.getPosX(), hit.getPosY());
	        getStatistics().fillHistogram("All_ToT_syn", hit.getEnergy());
          getStatistics().fillHistogram("All_ToT", HitFinderTools::calculateTOT(hit,HitFinderTools::getTOTCalculationType(fTOTCalculationType)));
          getStatistics().fillHistogram("All_Z_ToT", hit.getEnergy(),  hit.getPosZ());
          getStatistics().fillHistogram("All_mult",event.getHits().size());	  
          if(event.getHits().size()==4){ 
            energy_syn += hit.getEnergy();
            energy += HitFinderTools::calculateTOT(hit,HitFinderTools::getTOTCalculationType(fTOTCalculationType));}
        }
      if(event.getHits().size()==4){ 
        getStatistics().fillHistogram("All_ToT_sum", energy);
        getStatistics().fillHistogram("All_ToT_sum_syn", energy_syn);}
    }

    fIsSecondary = kFALSE;
    fIsScattered = kFALSE;
    fContainsPrompt = kFALSE;
    fIsPickOff = kFALSE;
    fIsOPs = kFALSE;

    events.push_back(newEvent);
    }
    //saveEvents(events);
  } else { return false; }
  for(int i = 0; i<n_oPs_true; i++) getStatistics().fillHistogram("evt_types", 1, 1);
  for(int i = 0; i<n_oPs_scatter; i++) getStatistics().fillHistogram("evt_types", 1, 2);
  for(int i = 0; i<n_oPs_prompt_true; i++) getStatistics().fillHistogram("evt_types", 2, 1);
  for(int i = 0; i<n_oPs_prompt_scatter; i++) getStatistics().fillHistogram("evt_types", 2, 2);
  for(int i = 0; i<n_prompt_true; i++) getStatistics().fillHistogram("evt_types", 3, 1);
  for(int i = 0; i<n_prompt_scatter; i++) getStatistics().fillHistogram("evt_types", 3, 2);
  for(int i = 0; i<n_random; i++) getStatistics().fillHistogram("evt_types", 4, 1);
  //cout<<n_oPs_true<<", "<<n_oPs_prompt_true<<", "<<n_oPs_scatter<<", "<<n_prompt_true<<", "<<n_random<<endl;
  return true;
}


bool EventCategorizer::terminate()
{
  INFO("Event categorization completed.");
  return true;
}

void EventCategorizer::saveEvents(const std::vector<JPetEvent>& events)
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
    new TH2D("All_Z_ToT", "Hit position Z vs ToT", 200, 0, fHistoLimit, 201, -25, 25),
    fEnergy_unit.c_str(), "Hit Z position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("All_ToT", fEnergy_unit.c_str(), 200, 0, fHistoLimit),
    fEnergy_unit.c_str(), "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("All_ToT_syn", fEnergy_unit.c_str(), 200, fHistoLimit*1e-3, fHistoLimit),
    fEnergy_unit.c_str(), "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("All_ToT_sum", "#sum ToT", 400, -1*fHistoLimit, 3*fHistoLimit),
    fEnergy_unit.c_str(), "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("All_ToT_sum_syn", "#sum ToT", 400, -1*fHistoLimit, 3*fHistoLimit),
    fEnergy_unit.c_str(), "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("All_mult", "Multiplicity", 11, -0.5, 10.5),
    "Multiplicity", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(new TH2D("evt_types", "Categories of events", 4, 0.5, 4.5, 2, 0.5, 2.5), " ", " ");
  vector<pair<unsigned, string>> binLabels_x = {make_pair(1, "oPs"), make_pair(2, "oPs+Prompt"), make_pair(3, "Prompt"),
                                              make_pair(4, "Randoms")};
  vector<pair<unsigned, string>> binLabels_y = {make_pair(1, "true"), make_pair(2, "true + scattered")};
  getStatistics().setHistogramBinLabel("evt_types", getStatistics().AxisLabel::kXaxis, binLabels_x);
  getStatistics().setHistogramBinLabel("evt_types", getStatistics().AxisLabel::kYaxis, binLabels_y);

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
    new TH1D("2Gamma_ToT", fEnergy_unit.c_str(), 200, 0, fHistoLimit),
    fEnergy_unit.c_str(), "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("2Gamma_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("2Gamma_Z_ToT", "Hit position Z vs ToT", 200, 0, fHistoLimit, 201, -23.25, 23.25),
    fEnergy_unit.c_str(), "Hit Z position [cm]"
  );


  // Histograms for 3Gamama category
  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_Angles", "Relative angles - transformed", 250, -0.5, 249.5, 200, -0.5, 199.5),
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
    new TH1D("3Gamma_ToT", fEnergy_unit.c_str(), 200, 0, fHistoLimit),
    fEnergy_unit.c_str(), "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("3Gamma_true_Energy", "Energy", 200, 0, fHistoLimit),
    fEnergy_unit.c_str(), "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("3Gamma_ToT_vertex_cut", fEnergy_unit.c_str(), 200, 0, fHistoLimit),
    fEnergy_unit.c_str(), "Number of Hits"
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
    new TH2D("3Gamma_Z_ToT", "Hit position Z vs ToT", 200, 0, fHistoLimit, 201, -23.25, 23.25),
    fEnergy_unit.c_str(), "Hit Z position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_true_Z_Energy", "Hit position Z vs Energy", 200, 0, fHistoLimit, 201, -23.25, 23.25),
    fEnergy_unit.c_str(), "Hit Z position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Lifetime_3Gamma", "Lifetime of 3 gamma events",
    200, -1700000, 1700000),
    "Lifetime [ps]", "Number of Events"
  );

    getStatistics().createHistogramWithAxes(
    new TH1D("Lifetime_3Gamma_true", "Lifetime of 3 gamma events",
    200, -1700000, 1700000),
    "Lifetime [ps]", "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Lifetime_3Gamma_zoom", "Lifetime of 3 gamma events",
    200, -100000, 100000),
    "Lifetime [ps]", "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Lifetime_3Gamma_vertex_reconstruction", "Lifetime of 3 gamma events",
    200, -1700000, 1700000),
    "Lifetime [ps]", "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Lifetime_3Gamma_vertex_reconstruction_zoom", "Lifetime of 3 gamma events",
    200, -100000, 100000),
    "Lifetime [ps]", "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Lifetime_3Gamma_vertex_cut", "Lifetime of 3 gamma events",
    200, -1700000, 1700000),
    "Lifetime [ps]", "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Lifetime_3Gamma_vertex_cut_zoom", "Lifetime of 3 gamma events",
    200, -100000, 100000),
    "Lifetime [ps]", "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Lifetime_3Gamma_vertex_cut_larger", "Lifetime of 3 gamma events",
    200, -1700000, 1700000),
    "Lifetime [ps]", "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Lifetime_3Gamma_vertex_cut_larger_zoom", "Lifetime of 3 gamma events",
    200, -100000, 100000),
    "Lifetime [ps]", "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_AnnihPoint_XY", "XY position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "X position [cm]", "Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_AnnihPoint_ZX", "ZX position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Z position [cm]", "X position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_AnnihPoint_ZY", "ZY position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Z position [cm]", "Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("decay point XY", "transverse position of the o-Ps->3g decay point", 200, -50., 50., 200, -50., 50.),
    "X [cm]", "Y [cm]"
	);

  getStatistics().createHistogramWithAxes(
    new TH2D("decay point XZ", "position of the o-Ps->3g decay point in XZ", 200, -50., 50.,200, -50., 50.),
    "Z [cm]", "X [cm]"
	);

getStatistics().createHistogramWithAxes(
    new TH1D("decay point time", "time of the o-Ps->3g decay point", 250, -20000., 20000.),
    "t [ps]", "Number of events"
	);

getStatistics().createHistogramWithAxes(
    new TH1D("decay point error", "error of the o-Ps->3g decay point", 5, -0.5, 4.5),
    "error", "Number of events"
	);

getStatistics().createHistogramWithAxes(
    new TH2D("decay point XY all", "transverse position of the o-Ps->3g decay point", 200, -50., 50., 200, -50., 50.),
    "X [cm]", "Y [cm]"
        );

getStatistics().createHistogramWithAxes(
    new TH2D("decay point XZ all", "position of the o-Ps->3g decay point in XZ", 200, -50., 50.,200, -50., 50.),
    "Z [cm]", "X [cm]"
        );

getStatistics().createHistogramWithAxes(
    new TH1D("decay point time all", "time of the o-Ps->3g decay point", 500, -50000., 50000.),
    "t [ps]", "Number of events"
        );

getStatistics().createHistogramWithAxes(
    new TH2D("decay point XY vertex cut", "transverse position of the o-Ps->3g decay point", 200, -50., 50., 200, -50., 50.),
    "X [cm]", "Y [cm]"
        );

getStatistics().createHistogramWithAxes(
    new TH2D("decay point XZ vertex cut", "position of the o-Ps->3g decay point in XZ", 200, -50., 50.,200, -50., 50.),
    "Z [cm]", "X [cm]"
        );

getStatistics().createHistogramWithAxes(
    new TH1D("decay point time vertex cut", "time of the o-Ps->3g decay point", 500, -50000., 50000.),
    "t [ps]", "Number of events"
        );

getStatistics().createHistogramWithAxes(
    new TH2D("decay point XY vertex cut larger", "transverse position of the o-Ps->3g decay point", 200, -50., 50., 200, -50., 50.),
    "X [cm]", "Y [cm]"
        );

getStatistics().createHistogramWithAxes(
    new TH2D("decay point XZ vertex cut larger", "position of the o-Ps->3g decay point in XZ", 200, -50., 50.,200, -50., 50.),
    "Z [cm]", "X [cm]"
        );

getStatistics().createHistogramWithAxes(
    new TH1D("decay point time vertex cut larger", "time of the o-Ps->3g decay point", 500, -50000., 50000.),
    "t [ps]", "Number of events"
        );

getStatistics().createHistogramWithAxes(
    new TH1D("DOP all", "Distance between anihhilation plane and decay point", 101, -0.5, 5.5),
    "DOP [cm]", "Number of Hits"
        );

getStatistics().createHistogramWithAxes(
    new TH1D("DOP", "Distance between anihhilation plane and decay point", 101, -0.5, 5.5),
    "DOP [cm]", "Number of Hits"
        );

getStatistics().createHistogramWithAxes(
    new TH1D("DOP vertex cut", "Distance between anihhilation plane and decay point", 101, -0.5, 5.5),
    "DOP [cm]", "Number of Hits"
        );

getStatistics().createHistogramWithAxes(
    new TH1D("DOP vertex cut larger", "Distance between anihhilation plane and decay point", 101, -0.5, 5.5),
    "DOP [cm]", "Number of Hits"
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
    new TH1D("ScatterTOF_TimeDiff_all", "Difference of Scatter TOF and hits time difference",
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
    new TH2D("Scatt_TOF_vs_TimeDiff_all", "Scatter TOF vs hits time difference for all events with 4 hits",
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
    new TH2D("Scatt_Dist_vs_M", "Difference of Scatter TOF and hits time difference vs distance",
    120, 0, 120, 3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Distance [cm]", "Scat_TOF - time diff [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_Dist_vs_M_all", "Difference of Scatter TOF and hits time difference vs distance for all events with 4 hits",
    120, 0, 120, 3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Distance [cm]", "Scat_TOF - time diff [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_Dist_vs_M_before_without_ID_cut", "Difference of Scatter TOF and hits time difference vs distance",
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
    new TH2D("Scatt_Dist_2D_vs_M", "Difference of Scatter TOF and hits time difference vs distance",
    50, 0, 50, 3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff-0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Distance [cm]", "Scat_TOF - time diff [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatt_Dist_2D_vs_M_before", "Difference of Scatter TOF and hits time difference vs distance",
    50, 0, 50, 3.0*fScatterTOFTimeDiff, -3.0*fScatterTOFTimeDiff+0.5, 3.0*fScatterTOFTimeDiff-0.5),
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
    new TH1D("Scatter_ToT", fEnergy_unit.c_str(), 200, -1*fHistoLimit*1e-3, fHistoLimit),
    fEnergy_unit.c_str(), "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Scatter_ToT_all", "ToT for all events with 4 hits", 200, -1*fHistoLimit*1e-3, fHistoLimit),
    fEnergy_unit.c_str(), "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("smallest_difference_Scatter_ToT", fEnergy_unit.c_str(), 200, 0, fHistoLimit),
    fEnergy_unit.c_str(), "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Scatter_ToT_sum", "#sum ToT", 400, -1*fHistoLimit, 3*fHistoLimit),
    fEnergy_unit.c_str(), "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Scatter_ToT_sum_syn", "#sum ToT", 400, -1*fHistoLimit, 3*fHistoLimit),
    fEnergy_unit.c_str(), "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("smallest_difference_Scatter_ToT_sum", "#sum ToT", 400, -1*fHistoLimit, 3*fHistoLimit),
    fEnergy_unit.c_str(), "Number of Events"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("smallest_difference_Scatter_ToT_sum_syn", "#sum ToT", 400, -1*fHistoLimit, 3*fHistoLimit),
    fEnergy_unit.c_str(), "Number of Events"
  );


  getStatistics().createHistogramWithAxes(
    new TH2D("Scatter_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("smallest_difference_Scatter_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scatter_Z_ToT", "Hit position Z vs ToT", 200, 0, fHistoLimit, 201, -23.25, 23.25),
    fEnergy_unit.c_str(), "Hit Z position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("smallest_difference_Scatter_Z_ToT", "Hit position Z vs ToT", 200, 0, fHistoLimit, 201, -23.25, 23.25),
    fEnergy_unit.c_str(), "Hit Z position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Scattered_true_Energy", "Energy", 200, 0, fHistoLimit),
    fEnergy_unit.c_str(), "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scattered_true_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Scattered_true_Z_Energy", "Hit position Z vs Energy", 200, 0, 2e3, 201, -23.25, 23.25),
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
    new TH1D("Deex_TOT_cut", "TOT of all hits with deex cut", 200, 0.5*fHistoLimit, fHistoLimit),
    fEnergy_unit.c_str(), "Number of Hits"
  );

  // Histograms for accidentals

  getStatistics().createHistogramWithAxes(
    new TH1D("Acc_Energy", fEnergy_unit.c_str(), 200, 0, fHistoLimit),
    fEnergy_unit.c_str(), "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Acc_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("Acc_Z_Energy", "Hit position Z vs ToT", 200, 0, fHistoLimit, 201, -23.25, 23.25),
    fEnergy_unit.c_str(), "Hit Z position [cm]"
  );

}
