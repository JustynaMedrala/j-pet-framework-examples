/**
 *  @copyright Copyright 2024 The J-PET Framework Authors. All rights reserved.
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

#include "EventCategorizer.h"
#include "CalibrationTools.h"
#include "EventCategorizerTools.h"
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>

using namespace jpet_options_tools;
using namespace std;

EventCategorizer::EventCategorizer(const char* name) : JPetUserTask(name) {}

EventCategorizer::~EventCategorizer() {}

bool EventCategorizer::init()
{
  INFO("Event categorization started.");

  std::cout<<"init"<<std::endl;

  // Reading user parameters
  if (isOptionSet(fParams.getOptions(), kEventTimeParamKey))
  {
    fEventTimeWindow = getOptionAsDouble(fParams.getOptions(), kEventTimeParamKey);
  }

  if (isOptionSet(fParams.getOptions(), kEventTimeZoomParamKey))
  {
    fEventTimeWindow_zoom = getOptionAsDouble(fParams.getOptions(), kEventTimeZoomParamKey);
  }

  if (isOptionSet(fParams.getOptions(), k2gThetaDiffParamKey))
  {
    f2gThetaDiff = getOptionAsDouble(fParams.getOptions(), k2gThetaDiffParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", k2gThetaDiffParamKey.c_str(), f2gThetaDiff));
  }

  if (isOptionSet(fParams.getOptions(), k2gTimeDiffParamKey))
  {
    f2gTimeDiff = getOptionAsDouble(fParams.getOptions(), k2gTimeDiffParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", k2gTimeDiffParamKey.c_str(), f2gTimeDiff));
  }

  // 3 gamma selection
  if (isOptionSet(fParams.getOptions(), k3gMinRelAngleParamKey))
  {
    f3gMinRelAngle = getOptionAsDouble(fParams.getOptions(), k3gMinRelAngleParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", k3gMinRelAngleParamKey.c_str(), f3gMinRelAngle));
  }

  // Reading ToT cut values
  if (isOptionSet(fParams.getOptions(), kToTCutAnniMinParamKey))
  {
    fToTCutAnniMin = getOptionAsDouble(fParams.getOptions(), kToTCutAnniMinParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kToTCutAnniMinParamKey.c_str(), fToTCutAnniMin));
  }

  if (isOptionSet(fParams.getOptions(), kToTCutAnniMaxParamKey))
  {
    fToTCutAnniMax = getOptionAsDouble(fParams.getOptions(), kToTCutAnniMaxParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kToTCutAnniMaxParamKey.c_str(), fToTCutAnniMax));
  }

  if (isOptionSet(fParams.getOptions(), kToTCutDeexMinParamKey))
  {
    fToTCutDeexMin = getOptionAsDouble(fParams.getOptions(), kToTCutDeexMinParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kToTCutDeexMinParamKey.c_str(), fToTCutDeexMin));
  }

  if (isOptionSet(fParams.getOptions(), kToTCutDeexMaxParamKey))
  {
    fToTCutDeexMax = getOptionAsDouble(fParams.getOptions(), kToTCutDeexMaxParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kToTCutDeexMaxParamKey.c_str(), fToTCutDeexMax));
  }

  // For plotting ToT histograms
  if (isOptionSet(fParams.getOptions(), kToTHistoUpperLimitParamKey))
  {
    fToTHistoUpperLimit = getOptionAsDouble(fParams.getOptions(), kToTHistoUpperLimitParamKey);
  }

  // Cuts around source position
  if (isOptionSet(fParams.getOptions(), kSourceDistCutXYParamKey))
  {
    fSourceDistXYCut = getOptionAsDouble(fParams.getOptions(), kSourceDistCutXYParamKey);
  }
  else
  {
    WARNING(
        Form("No value of the %s parameter provided by the user. Using default value of %lf.", kSourceDistCutXYParamKey.c_str(), fSourceDistXYCut));
  }

  if (isOptionSet(fParams.getOptions(), kSourceDistCutZParamKey))
  {
    fSourceDistZCut = getOptionAsDouble(fParams.getOptions(), kSourceDistCutZParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kSourceDistCutZParamKey.c_str(), fSourceDistZCut));
  }

  // LOR cuts
  if (isOptionSet(fParams.getOptions(), kLORAngleCutParamKey))
  {
    fLORAngleCut = getOptionAsDouble(fParams.getOptions(), kLORAngleCutParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kLORAngleCutParamKey.c_str(), fLORAngleCut));
  }

  if (isOptionSet(fParams.getOptions(), kLORPosZCutParamKey))
  {
    fLORPosZCut = getOptionAsDouble(fParams.getOptions(), kLORPosZCutParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kLORPosZCutParamKey.c_str(), fLORPosZCut));
  }

  // Source position
  if (isOptionSet(fParams.getOptions(), kSourcePosXParamKey) && isOptionSet(fParams.getOptions(), kSourcePosYParamKey) &&
      isOptionSet(fParams.getOptions(), kSourcePosZParamKey))
  {
    auto x = getOptionAsDouble(fParams.getOptions(), kSourcePosXParamKey);
    auto y = getOptionAsDouble(fParams.getOptions(), kSourcePosYParamKey);
    auto z = getOptionAsDouble(fParams.getOptions(), kSourcePosZParamKey);
    fSourcePos.SetXYZ(x, y, z);
    INFO(Form("Source position is: %lf, %lf, %lf", x, y, z));
  }
  else
  {
    fSourcePos.SetXYZ(0.0, 0.0, 0.0);
    INFO("Source is positioned in (0, 0, 0).");
  }

  // Reading file with constants to property tree
  if (isOptionSet(fParams.getOptions(), kConstantsFileParamKey))
  {
    boost::property_tree::read_json(getOptionAsString(fParams.getOptions(), kConstantsFileParamKey), fConstansTree);
  }

  // Set the type of scatter test - default is simple parameter cut
  if (isOptionSet(fParams.getOptions(), kScatterTestTypeParamKey))
  {
    auto type = getOptionAsString(fParams.getOptions(), kScatterTestTypeParamKey);
    if (type == "simple_param")
    {
      fTestType = EventCategorizerTools::kSimpleParam;
    }
    else if (type == "min_max")
    {
      fTestType = EventCategorizerTools::kMinMaxParams;
    }
  }

  if (isOptionSet(fParams.getOptions(), kScatterTOFTimeDiffParamKey))
  {
    fScatterTOFTimeDiff = getOptionAsDouble(fParams.getOptions(), kScatterTOFTimeDiffParamKey);
  }

  if (isOptionSet(fParams.getOptions(), kScatterTimeMinParamKey))
  {
    fScatterTimeMin = getOptionAsDouble(fParams.getOptions(), kScatterTimeMinParamKey);
  }
  if (isOptionSet(fParams.getOptions(), kScatterTimeMaxParamKey))
  {
    fScatterTimeMax = getOptionAsDouble(fParams.getOptions(), kScatterTimeMaxParamKey);
  }

  if (isOptionSet(fParams.getOptions(), kScatterAngleMinParamKey))
  {
    fScatterAngleMin = getOptionAsDouble(fParams.getOptions(), kScatterAngleMinParamKey);
  }
  if (isOptionSet(fParams.getOptions(), kScatterAngleMaxParamKey))
  {
    fScatterAngleMax = getOptionAsDouble(fParams.getOptions(), kScatterAngleMaxParamKey);
  }

  // Time variable used as +- axis limits for histograms with time spectra
  if (isOptionSet(fParams.getOptions(), kMaxTimeDiffParamKey))
  {
    fMaxTimeDiff = getOptionAsDouble(fParams.getOptions(), kMaxTimeDiffParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kMaxTimeDiffParamKey.c_str(), fMaxTimeDiff));
  }

  if (isOptionSet(fParams.getOptions(), k3DOPParamKey))
  {
    f3gDOP = getOptionAsDouble(fParams.getOptions(), k3DOPParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kMaxTimeDiffParamKey.c_str(), fMaxTimeDiff));
  }

  // Getting bools for saving histograms
  if (isOptionSet(fParams.getOptions(), kSaveControlHistosParamKey))
  {
    fSaveControlHistos = getOptionAsBool(fParams.getOptions(), kSaveControlHistosParamKey);
  }
  if (isOptionSet(fParams.getOptions(), kSaveCalibHistosParamKey))
  {
    fSaveCalibHistos = getOptionAsBool(fParams.getOptions(), kSaveCalibHistosParamKey);
  }
  if (isOptionSet(fParams.getOptions(), kTrentoCalibrationParamKey))
  {
    fTrentoCalibHistos = getOptionAsBool(fParams.getOptions(), kTrentoCalibrationParamKey);
  }
  if (isOptionSet(fParams.getOptions(), kDataTypeParamKey))
  {
    fDataType = getOptionAsBool(fParams.getOptions(), kDataTypeParamKey);
  }

  if (fTrentoCalibHistos)
  {
    // Variable used in measurements with Trento setup - rotation of Z axis with respect of
    // vertical direction in Earth frame
    if (isOptionSet(fParams.getOptions(), kDetectorYRotation))
    {
      fDetectorYRotationDeg = getOptionAsDouble(fParams.getOptions(), kDetectorYRotation);
    }
    else
    {
      WARNING(
          Form("No value of the %s parameter provided by the user. Using default value of %lf.", kDetectorYRotation.c_str(), fDetectorYRotationDeg));
    }

    if (isOptionSet(fParams.getOptions(), kCosmicMaxThetaDeg))
    {
      fCosmicMaxThetaDiffDeg = getOptionAsDouble(fParams.getOptions(), kCosmicMaxThetaDeg);
    }
    else
    {
      WARNING(
          Form("No value of the %s parameter provided by the user. Using default value of %lf.", kCosmicMaxThetaDeg.c_str(), fCosmicMaxThetaDiffDeg));
    }
  }

  std::cout<<"init histo"<<std::endl;
  // Initialise hisotgrams
  if (fSaveControlHistos)
  {
    initialiseHistograms(fDataType);
  }

  if (fSaveCalibHistos)
  {
    initialiseCalibrationHistograms(fTrentoCalibHistos);
  }

  // Input events type
  fOutputEvents = new JPetTimeWindow("JPetEvent");

  return true;
}

bool EventCategorizer::exec()
{
  std::cout<<"exec"<<std::endl;
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent))
  {
    std::cout<<"loop"<<std::endl;
    vector<JPetEvent> events;
    for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++)
    {
      const auto& event = dynamic_cast<const JPetEvent&>(timeWindow->operator[](i));


      vector<int> bad_ID = {206, 219, 222, 232, 258, 271, 284, 293, 297, 300, 302, 310, 313, 336, 349, 361, 375, 378, 388, 414, 427, 440, 
                            449, 452, 456, 458, 466, 469, 492, 505};

      // Categorization of the events
      //bool is1Gamma = EventCategorizerTools::checkFor1Gamma(event, fToTCut1AnniMin, fToTCut1AnniMax, fToTCut1AnniMin_larger, fToTCut1AnniMax_larger, 
      //                                                      fToTCutDeexMin, fToTCutDeexMax, getStatistics(), fSaveControlHistos);
      
      /*bool is2Gamma = EventCategorizerTools::checkFor2Gamma(event, getStatistics(), fSaveControlHistos, f2gThetaDiff, f2gTimeDiff, fToTCutAnniMin,
                                                            fToTCutAnniMax, fSourcePos, fTestType, fScatterTOFTimeDiff, fScatterTimeMin,
                                                            fScatterTimeMax, fScatterAngleMin, fScatterAngleMax);

      bool is3Gamma = EventCategorizerTools::checkFor3Gamma(event, f3gMinRelAngle, f2gTimeDiff, fToTCut3AnniMin, fToTCut3AnniMax, getStatistics(), 
                                                            fSaveControlHistos);
      */
      /*bool isLifetime2Gamma = EventCategorizerTools::checkFor2GammaLifetime(
          event, getStatistics(), fSaveControlHistos, f2gThetaDiff, f2gTimeDiff, fToTCutAnniMin, fToTCutAnniMax, fToTCutDeexMin, fToTCutDeexMax,
          fSourcePos, fTestType, fScatterTOFTimeDiff, fScatterTimeMin, fScatterTimeMax, fScatterAngleMin, fScatterAngleMax);
      */
      std::cout<<"isLifetime2Gamma"<<std::endl;
      bool isLifetime2Gamma = EventCategorizerTools::checkFor2GammaLifetime(
          event, bad_ID, getStatistics(), fSaveControlHistos, f2gThetaDiff, f2gTimeDiff, f2gDOP, fToTCutAnniMin, fToTCutAnniMax, fToTCutDeexMin, fToTCutDeexMax,
          fSourcePos, fTestType, fScatterTOFTimeDiff, fScatterTimeMin, fScatterTimeMax, fScatterAngleMin, fScatterAngleMax);
      
      /*bool isLifetime2Gamma_good = EventCategorizerTools::checkFor2GammaLifetime_exactly_3hits(
          event, bad_ID, getStatistics(), fSaveControlHistos, f2gThetaDiff, f2gTimeDiff, fToTCutAnniMin, fToTCutAnniMax, fToTCutDeexMin, fToTCutDeexMax,
          fSourcePos, fTestType, fScatterTOFTimeDiff, fScatterTimeMin, fScatterTimeMax, fScatterAngleMin, fScatterAngleMax);


      bool isLifetime3Gamma = EventCategorizerTools::checkFor3GammaLifetime(
          event, f3gMinRelAngle, f3gTimeDiff, getStatistics(), fSaveControlHistos, fToTCut3AnniMin, fToTCut3AnniMax, fToTCutDeexMin, fToTCutDeexMax, fTestType, 
          fScatterTOFTimeDiff, fScatterTimeMin, fScatterTimeMax, fScatterAngleMin, fScatterAngleMax);
      */
      std::cout<<"isLifetime3Gamma"<<std::endl;
      bool isLifetime3Gamma = EventCategorizerTools::checkFor3GammaLifetime(
          event, bad_ID, f3gMinRelAngle, f3gMinRelPhi, f3gMinDist, f3gTimeDiff, f3gDOP, getStatistics(), fSaveControlHistos, fToTCut3AnniMin, fToTCut3AnniMax, fToTCutDeexMin, fToTCutDeexMax, fSourcePos, fTestType, 
          fScatterTOFTimeDiff, fScatterTimeMin, fScatterTimeMax, f3gScatterAngleMin, f3gScatterAngleMax);

      std::cout<<"Unknow"<<std::endl;
      JPetEvent newEvent = event;

      /*if (isLifetime2Gamma_good)
      {
        getStatistics().fillHistogram("stats_2gamma_masking", 2);
      }

      if (isLifetime3Gamma_good)
      {
        getStatistics().fillHistogram("stats_3gamma_masking", 2);
      }

      if (is1Gamma)
      {
        newEvent.addEventType(JPetEventType::kPrompt);
        getStatistics().fillHistogram("evt_types", 6);
      }
      if (is2Gamma)
      {
        newEvent.addEventType(JPetEventType::k2Gamma);
        getStatistics().fillHistogram("evt_types", 2);
      }
      if (is3Gamma)
      {
        newEvent.addEventType(JPetEventType::k3Gamma);
        getStatistics().fillHistogram("evt_types", 3);
      }*/
      if (isLifetime2Gamma)
      {
        newEvent.addEventType(JPetEventType::kPrompt);
        getStatistics().fillHistogram("evt_types", 4);
      }

      if (isLifetime3Gamma)
      {
        newEvent.addEventType(JPetEventType::kPrompt);
        getStatistics().fillHistogram("evt_types", 5);
      }

      if (newEvent.isOnlyTypeOf(JPetEventType::kUnknown))
      {
        getStatistics().fillHistogram("evt_types", 1);
      }
      else
      {
        // Saving the event only if it was categorized
        events.push_back(newEvent);
      }
    }
    //saveEvents(events); do not save the tree
  }
  else
  {
    return false;
  }
  return true;
}

bool EventCategorizer::terminate()
{
  INFO("Event categorization completed.");
  return true;
}

void EventCategorizer::saveEvents(const vector<JPetEvent>& events)
{
  for (const auto& event : events)
  {
    fOutputEvents->add<JPetEvent>(event);
  }
}

void EventCategorizer::initialiseHistograms(bool dataType)
{
  int minScinID = 200;
  int maxScinID = 512;

  std::string energy_units = "", energy_hist_title = "";

  if(!dataType){ energy_units = "Time over Threshold [ps]";
  energy_hist_title = "average ToT scaled";}
  if(dataType){ energy_units = "Energy [keV]";
  energy_hist_title = "Energy";}

  // Event categories
  getStatistics().createHistogramWithAxes(new TH1D("evt_types", "Categories of events", 6, 0.5, 6.5), " ", "Number of events");
  vector<pair<unsigned, string>> binLabels = {make_pair(1, "Unknown"), make_pair(2, "2 gamma"), make_pair(3, "3 gamma"),
                                              make_pair(4, "2 gamma + prompt"), make_pair(5, "3 gamma + prompt"), make_pair(6, "1 gamma + prompt")};
  getStatistics().setHistogramBinLabel("evt_types", getStatistics().AxisLabel::kXaxis, binLabels);

  /*getStatistics().createHistogramWithAxes(new TH1D("1g_tot", "average ToT scaled", 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("1g_prompt_tot", "average ToT scaled", 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("1g_tot_tot", "average ToT scaled", 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("1g_tot_prompt_tot", "average ToT scaled", 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("1g_tot_larger_tot", "average ToT scaled", 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("1g_tot_larger_prompt_tot", "average ToT scaled", 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");
*/
  // Histograms for 2 gamama events
    getStatistics().createHistogramWithAxes(new TH1D("none_2g_tot", ("2 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");
                                    
    getStatistics().createHistogramWithAxes(new TH1D("hits_2g_tot", ("2 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");
                                        
    getStatistics().createHistogramWithAxes(new TH1D("2g_tot", ("2 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

    getStatistics().createHistogramWithAxes(new TH1D("2g_ID_all", "ID of all hits", 314, 199.5, 513.5),
                                           "ID", "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("2g_theta", "2 gamma event - flight vectors theta", 181, -0.5, 180.5), "Angle [degree]",
                                          "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH2D("2g_timeDiff_theta", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff, 181, -0.5, 180.5),
      "Time  Difference between annihilation hits [ps]", "Angle [degree]");

  getStatistics().createHistogramWithAxes(new TH1D("2g_DOP", "2 gamma event - distance between source and the reconstructed annihilation point", 100, 0.0,100.0),
                                          "Distance [cm]", "Number of Hits");

  getStatistics().createHistogramWithAxes(
      new TH1D("2g_timeDiff", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH1D("2g_masking_timeDiff", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

    getStatistics().createHistogramWithAxes(
      new TH2D("2g_timeDiff_ID", "Time difference between hits", maxScinID - minScinID + 1, minScinID - 0.5, maxScinID + 0.5, 201, -fMaxTimeDiff, fMaxTimeDiff),
      "ScinID", "Time  Difference between annihilation hits [ps]");

  getStatistics().createHistogramWithAxes(new TH2D("2g_time_annih_12", "Time of registration", 500, 0, 5e7, 500, 0, 5e7),
                                          "annih_1 registration time [ps]", "annih_2 registration time [ps]");

  getStatistics().createHistogramWithAxes(new TH2D("2g_time_hit_12", "Time of registration", 500, 0, 5e7, 500, 0, 5e7),
                                          "time of the first hit in the event [ps]", "time of the second hit in the event [ps]");

    getStatistics().createHistogramWithAxes(new TH1D("2g_ID_peak", "ID of hits in time difference range [3.5, 5] ns", 314, 199.5, 513.5),
                                          "ID", "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH2D("2g_xy_peak", "XY position of hits (bin 0.5 cm)", 242, -60.5, 60.5, 242, -60.5, 60.5),
                                          "X position [cm]", "Y position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("2g_time_annih_12_peak", "Time of registration", 500, 0, 5e7, 500, 0, 5e7),
                                          "annih_1 registration time [ps]", "annih_2 registration time [ps]");

  getStatistics().createHistogramWithAxes(new TH2D("2g_zx_peak", "ZX position of hits (bin 0.5 cm)", 242, -60.5, 60.5, 242, -60.5, 60.5),
                                          "Z position [cm]", "X position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("2g_xy_source_peak", "XY position of annihilation point (bin 0.5 cm)", 242, -60.5, 60.5, 242, -60.5, 60.5),
                                          "X position [cm]", "Y position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("2g_zx_source_peak", "ZX position of annihilation point (bin 0.5 cm)",242, -60.5, 60.5, 242, -60.5, 60.5),
                                          "Z position [cm]", "X position [cm]");

  getStatistics().createHistogramWithAxes(new TH1D("2g_ID_remain", "ID of all hits", 314, 199.5, 513.5),
                                           "ID", "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("2g_scatter_test_time", "Scatter Test - Time Difference", 401, -10000.0, 10000.0),
                                          "annih1 - annih2, tDiff - d/c [ps]", "Number of pairs");
  
  getStatistics().createHistogramWithAxes(
      new TH1D("scatter_2g_tot", ("2 gamma event after theta cut - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit), energy_units.c_str(),
      "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("scatter_2g_theta", "2 gamma event after ToT cut - theta between flight vectors", 181, -0.5, 180.5),
                                          "Angle [degree]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH1D("scatter_2g_timeDiff", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH2D("scatter_2g_xy_peak", "XY position of hits (bin 0.5 cm)", 242, -60.5, 60.5, 242, -60.5, 60.5),
                                          "X position [cm]", "Y position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("scatter_2g_zx_peak", "ZX position of hits (bin 0.5 cm)", 242, -60.5, 60.5, 242, -60.5, 60.5),
                                          "Z position [cm]", "X position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("scatter_2g_xy_source_peak", "XY position of annihilation point (bin 0.5 cm)", 242, -60.5, 60.5, 242, -60.5, 60.5),
                                          "X position [cm]", "Y position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("scatter_2g_zx_source_peak", "ZX position of annihilation point (bin 0.5 cm)", 242, -60.5, 60.5, 242, -60.5, 60.5),
                                          "Z position [cm]", "X position [cm]");

  getStatistics().createHistogramWithAxes(new TH1D("scatter_2g_scatter_test_time", "Scatter Test - Time Difference", 401, -10000.0, 10000.0),
                                          "annih1 - annih2, tDiff - d/c [ps]", "Number of pairs");

getStatistics().createHistogramWithAxes(
      new TH1D("tdiff_2g_tot", ("2 gamma event after theta cut - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit), energy_units.c_str(),
      "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("tdiff_2g_theta", "2 gamma event after ToT cut - theta between flight vectors", 181, -0.5, 180.5),
                                          "Angle [degree]", "Number of Hit Pairs");
                                  
  getStatistics().createHistogramWithAxes(
      new TH1D("tdiff_2g_timeDiff", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");
    
  getStatistics().createHistogramWithAxes(new TH1D("tdiff_2g_scatter_test_time", "Scatter Test - Time Difference", 401, -10000.0, 10000.0),
                                          "annih1 - annih2, tDiff - d/c [ps]", "Number of pairs");

    getStatistics().createHistogramWithAxes(
      new TH1D("theta_2g_tot", ("2 gamma event after theta cut - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit), energy_units.c_str(),
      "Number of Hit Pairs");

    getStatistics().createHistogramWithAxes(new TH1D("theta_2g_ID_all", "ID of all hits", 314, 199.5, 513.5),
                                           "ID", "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("theta_2g_theta", "2 gamma event after ToT cut - theta between flight vectors", 181, -0.5, 180.5),
                                          "Angle [degree]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH2D("theta_2g_timeDiff_theta", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff, 181, -0.5, 180.5),
      "Time  Difference between annihilation hits [ps]", "Angle [degree]");

  getStatistics().createHistogramWithAxes(
      new TH1D("theta_2g_timeDiff", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH2D("theta_2g_time_annih_12", "Time of registration", 500, 0, 5e7, 500, 0, 5e7),
                                          "t1 [ps]", "t2 [ps]");

    getStatistics().createHistogramWithAxes(new TH1D("theta_2g_ID_peak", "ID of hits in time difference range [3.5, 5] ns", 314, 199.5, 513.5),
                                           "ID", "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH2D("theta_2g_xy_first_peak", "XY position of hits (bin 0.5 cm)", 242, -60.5, 60.5, 242, -60.5, 60.5),
                                          "X position [cm]", "Y position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("theta_2g_xy_second_peak", "XY position of hits (bin 0.5 cm)", 242, -60.5, 60.5, 242, -60.5, 60.5),
                                          "X position [cm]", "Y position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("theta_2g_zx_first_peak", "ZX position of hits (bin 0.5 cm)", 242, -60.5, 60.5, 242, -60.5, 60.5),
                                          "Z position [cm]", "X position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("theta_2g_zx_second_peak", "ZX position of hits (bin 0.5 cm)", 242, -60.5, 60.5, 242, -60.5, 60.5),
                                          "Z position [cm]", "X position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("theta_2g_xy_source_peak", "XY position of annihilation point (bin 0.5 cm)", 242, -60.5, 60.5, 242, -60.5, 60.5),
                                          "X position [cm]", "Y position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("theta_2g_zx_source_peak", "ZX position of annihilation point (bin 0.5 cm)", 242, -60.5, 60.5, 242, -60.5, 60.5),
                                          "Z position [cm]", "X position [cm]");

  getStatistics().createHistogramWithAxes(new TH1D("theta_2g_scatter_test_time", "Scatter Test - Time Difference", 401, -10000.0, 10000.0),
                                          "annih1 - annih2, tDiff - d/c [ps]", "Number of pairs");

getStatistics().createHistogramWithAxes(new TH1D("ap_2g_tot", ("2 gamma event after theta cut - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit), energy_units.c_str(),
      "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("ap_2g_theta", "Annihilation pairs flight vectors theta", 181, -0.5, 180.5), "Angle [degree]",
                                          "Number of Hit Pairs");
                                  
  getStatistics().createHistogramWithAxes(
      new TH1D("ap_2g_timeDiff", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH2D("ap_xy", "XY position of annihilation point (bin 0.5 cm)", 202, -50.5, 50.5, 202, -50.5, 50.5),
                                          "X position [cm]", "Y position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("ap_zx", "ZX position of annihilation point (bin 0.5 cm)", 202, -50.5, 50.5, 202, -50.5, 50.5),
                                          "Z position [cm]", "X position [cm]");


  getStatistics().createHistogramWithAxes(new TH2D("ap_zy", "ZY position of annihilation point (bin 0.5 cm)", 202, -50.5, 50.5, 202, -50.5, 50.5),
                                          "Z position [cm]", "Y position [cm]");

  getStatistics().createHistogramWithAxes(
      new TH3D("ap_pos", "Position of the annihilation point (bin 0.5 cm)", 202, -50.5, 50.5, 202, -50.5, 50.5, 202, -50.5, 50.5), "Z [cm]", "X [cm]",
      "Y [cm]");

  getStatistics().createHistogramWithAxes(
      new TH2D("ap_xy_zoom", "XY position of annihilation point (bin 0.2 cm)", 45, -4.5, 4.5, 45, -4.5, 4.5), "X position [cm]",
      "Y position [cm]");

  getStatistics().createHistogramWithAxes(
      new TH2D("ap_zx_zoom", "ZX position of annihilation point (bin 0.2 cm)", 45, -4.5, 4.5, 45, -4.5, 4.5), "Z position [cm]",
      "X position [cm]");

  
  getStatistics().createHistogramWithAxes(
      new TH2D("ap_zy_zoom", "ZY position of annihilation point (bin 0.2 cm)", 45, -4.5, 4.5, 45, -4.5, 4.5), "Z position [cm]",
      "Y position [cm]");

  getStatistics().createHistogramWithAxes(
      new TH3D("ap_pos_zoom", "Position of the annihilation point (bin 0.2 cm)", 45, -4.5, 4.5, 45, -4.5, 4.5, 45, -4.5, 4.5), "Z [cm]",
      "X [cm]", "Y [cm]");

  getStatistics().createHistogramWithAxes(new TH1D("ap_2g_scatter_test_time", "Scatter Test - Time Difference", 401, -10000.0, 10000.0),
                                          "annih1 - annih2, tDiff - d/c [ps]", "Number of pairs");

  //Histograms for 2gamma + 1 prompt category

    getStatistics().createHistogramWithAxes(
      new TH1D("ap_2g_tot_lifetime", ("2 gamma event after theta cut - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit), energy_units.c_str(),
      "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("ap_2g_theta_lifetime", "Annihilation pairs flight vectors theta", 181, -0.5, 180.5), "Angle [degree]",
                                          "Number of Hit Pairs");
                                  
  getStatistics().createHistogramWithAxes(
      new TH1D("ap_2g_timeDiff_lifetime", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH2D("ap_xy_lifetime", "XY position of annihilation point (bin 0.5 cm)", 202, -50.5, 50.5, 202, -50.5, 50.5),
                                          "X position [cm]", "Y position [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("ap_zx_lifetime", "ZX position of annihilation point (bin 0.5 cm)", 202, -50.5, 50.5, 202, -50.5, 50.5),
                                          "Z position [cm]", "X position [cm]");


  getStatistics().createHistogramWithAxes(new TH2D("ap_zy_lifetime", "ZY position of annihilation point (bin 0.5 cm)", 202, -50.5, 50.5, 202, -50.5, 50.5),
                                          "Z position [cm]", "Y position [cm]");

  getStatistics().createHistogramWithAxes(
      new TH3D("ap_pos_lifetime", "Position of the annihilation point (bin 0.5 cm)", 202, -50.5, 50.5, 202, -50.5, 50.5, 202, -50.5, 50.5), "Z [cm]", "X [cm]",
      "Y [cm]");

  getStatistics().createHistogramWithAxes(
      new TH2D("ap_xy_zoom_lifetime", "XY position of annihilation point (bin 0.2 cm)", 45, -4.5, 4.5, 45, -4.5, 4.5), "X position [cm]",
      "Y position [cm]");

  getStatistics().createHistogramWithAxes(
      new TH2D("ap_zx_zoom_lifetime", "ZX position of annihilation point (bin 0.2 cm)", 45, -4.5, 4.5, 45, -4.5, 4.5), "Z position [cm]",
      "X position [cm]");

  
  getStatistics().createHistogramWithAxes(
      new TH2D("ap_zy_zoom_lifetime", "ZY position of annihilation point (bin 0.2 cm)", 45, -4.5, 4.5, 45, -4.5, 4.5), "Z position [cm]",
      "Y position [cm]");

  getStatistics().createHistogramWithAxes(
      new TH3D("ap_pos_zoom_lifetime", "Position of the annihilation point (bin 0.2 cm)", 45, -4.5, 4.5, 45, -4.5, 4.5, 45, -4.5, 4.5), "Z [cm]",
      "X [cm]", "Y [cm]");

  getStatistics().createHistogramWithAxes(new TH1D("ap_2g_scatter_test_time_lifetime", "Scatter Test - Time Difference", 401, -10000.0, 10000.0),
                                          "annih1 - annih2, tDiff - d/c [ps]", "Number of pairs");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_2g_prompt", "Time difference of 2 gamma pair decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");
                            
  getStatistics().createHistogramWithAxes(new TH1D("lifetime_2g_prompt_zoom", "Time difference of 2 gamma pair decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_tdiff_2g_prompt", "Time difference of 2 gamma pair decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_tdiff_2g_prompt_zoom", "Time difference of 2 gamma pair decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_scatter_2g_prompt", "Time difference of 2 gamma pair decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_scatter_2g_prompt_zoom", "Time difference of 2 gamma pair decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_theta_2g_prompt", "Time difference of 2 gamma pair decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_theta_2g_prompt_zoom", "Time difference of 2 gamma pair decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_ap_2g_prompt", "Time difference of 2 gamma pair decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_ap_2g_prompt_zoom", "Time difference of 2 gamma pair decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 


  // Histograms for scattering category
  getStatistics().createHistogramWithAxes(new TH1D("scatter_test_dist_abs", "Scatter Test - Distance Difference", 201, 0.0, 120.0), "Dist Diff [cm]",
                                          "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("scatter_test_dist_rel", "Scatter Test - Distance Difference", 201, -120.0, 120.0),
                                          "Dist Diff [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("scatter_test_time_abs", "Scatter Test - Time Difference", 201, 0.0, 10000.0), 
                                          "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("scatter_test_time_rel", "Scatter Test - Time Difference", 201, -5000.0, 5000.0), 
                                          "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("scatter_test_abs_pass", "Passed Scatter Test - Time Difference", 201, 0.0, 10000.0),
                                          "Time Diff [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("scatter_test_abs_fail", "Failed Scatter Test - Time Difference", 201, 0.0, 10000.0),
                                          "Time Diff [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("scatter_test_rel_pass", "Passed Scatter Test - Time Difference", 201, -5000.0, 5000.0),
                                          "Time Diff [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("scatter_test_rel_fail", "Failed Scatter Test - Time Difference", 201, -5000.0, 5000.0),
                                          "Time Diff [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH2D("scatter_angle_time", "Scatter angle vs. scatter test measure", 201, -4000.0, 6000.0, 181, -0.5, 180.5), "Time Diff [ps]",
      "Scatter angle");

  getStatistics().createHistogramWithAxes(
      new TH2D("scatter_angle_time_small", "Scatter angle vs. scatter test measure", 201, -4000.0, 6000.0, 41, 139.5, 180.5), "Time Diff [ps]",
      "Scatter angle");

  getStatistics().createHistogramWithAxes(
      new TH2D("scatter_angle_time_pass", "Passed Scatter angle vs. scatter test measure", 201, -4000.0, 6000.0, 181, -0.5, 180.5), "Time Diff [ps]",
      "Scatter angle");

  getStatistics().createHistogramWithAxes(
      new TH2D("scatter_angle_time_fail", "Failed Scatter angle vs. scatter test measure", 201, -4000.0, 6000.0, 181, -0.5, 180.5), "Time Diff [ps]",
      "Scatter angle");

  // Histograms for 3 gamma events

/*  getStatistics().createHistogramWithAxes(new TH1D("3g_tot_all", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("3g_tot", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");


  getStatistics().createHistogramWithAxes(new TH1D("tot_3g_tot", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("tdiff_3g_tot", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("theta_3g_tot", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("ap_3g_tot", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");
  
  getStatistics().createHistogramWithAxes(
      new TH2D("3g_rel_angles", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH2D("tot_3g_rel_angles", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH2D("tdiff_3g_rel_angles", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH2D("theta_3g_rel_angles", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH2D("ap_3g_rel_angles", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(new TH2D("3g_rel_angles_sel",
                                          "Sum vs. difference of two smallest relative angles in 3 gamma event - after cut", 250, 0.0, 250,
                                          200, 0.0, 200.0), "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH1D("3g_timeDiff", "Maximal time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH1D("tot_3g_timeDiff", "Maximal time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH1D("tdiff_3g_timeDiff", "Maximal time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH1D("theta_3g_timeDiff", "Maximal time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH1D("ap_3g_timeDiff", "Maximal time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");
*/
  //3 gamma + prompt events
  getStatistics().createHistogramWithAxes(new TH1D("3g_stats_multi_prompt", "Number of prompts in event", 20, -0.5, 19.5),
                                          "Number of prompts", "Multiplicity");

  getStatistics().createHistogramWithAxes(new TH1D("3g_stats_multi_annihilations", "Number of annihilations in event", 20, -0.5, 19.5),
                                          "Number of 3 annihilations", "Multiplicity");

  getStatistics().createHistogramWithAxes(new TH1D("3g_stats_multi", "Number of hits in event", 20, -0.5, 19.5),
                                          "Number of hits", "Multiplicity");

  getStatistics().createHistogramWithAxes(new TH1D("3g_cut_stats", "Categories of events", 6, 0.5, 6.5), " ", "Number of events");
  binLabels = {make_pair(1, "All"), make_pair(2, "DOP"), make_pair(3, "time difference"),
                                              make_pair(4, "theta"), make_pair(5, "phi"), make_pair(6, "all")};
  getStatistics().setHistogramBinLabel("3g_cut_stats", getStatistics().AxisLabel::kXaxis, binLabels);

                                          
  getStatistics().createHistogramWithAxes(new TH1D("3g_scatter_test_time_lifetime", "Scatter Test - Time Difference", 801, -20000.0, 20000.0), 
                                          "Scatter test [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("DOP_3g_scatter_test_time_lifetime", "Scatter Test - Time Difference", 801, -20000.0, 20000.0), 
                                          "Scatter test [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("tdiff_3g_scatter_test_time_lifetime", "Scatter Test - Time Difference", 801, -20000.0, 20000.0), 
                                          "Scatter test [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("theta_3g_scatter_test_time_lifetime", "Scatter Test - Time Difference", 801, -20000.0, 20000.0), 
                                          "Scatter test [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("phi_3g_scatter_test_time_lifetime", "Scatter Test - Time Difference", 801, -20000.0, 20000.0), 
                                          "Scatter test [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("phi0_3g_scatter_test_time_lifetime", "Scatter Test - Time Difference", 801, -20000.0, 20000.0), 
                                          "Scatter test [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("dist_3g_scatter_test_time_lifetime", "Scatter Test - Time Difference", 801, -20000.0, 20000.0), 
                                          "Scatter test [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("vtx_3g_scatter_test_time_lifetime", "Scatter Test - Time Difference", 801, -20000.0, 20000.0), 
                                          "Scatter test [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("ap_all_3g_scatter_test_time_lifetime", "Scatter Test - Time Difference", 801, -20000.0, 20000.0), 
                                          "Scatter test [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("ap_3g_scatter_test_time_lifetime", "Scatter Test - Time Difference", 801, -20000.0, 20000.0), 
                                          "Scatter test [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("3g_scatter_test_time_prompt_lifetime", "Scatter Test - Time Difference", 801, -20000.0, 20000.0), 
                                          "Scatter test [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("3g_dist_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("DOP_3g_dist_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("tdiff_3g_dist_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("theta_3g_dist_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("phi_3g_dist_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("phi0_3g_dist_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("dist_3g_dist_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("vtx_3g_dist_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("ap_all_3g_dist_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("ap_3g_dist_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");


  getStatistics().createHistogramWithAxes(new TH1D("3g_dist2D_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("DOP_3g_dist2D_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("tdiff_3g_dist2D_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("theta_3g_dist2D_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("phi_3g_dist2D_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("phi0_3g_dist2D_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("dist_3g_dist2D_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");
  
  getStatistics().createHistogramWithAxes(new TH1D("vtx_3g_dist2D_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("ap_all_3g_dist2D_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("ap_3g_dist2D_lifetime", "Distance Difference", 111, -10.0, 100.0), 
                                          "Distance between annihilation hits [cm]", "Number of Hit Pairs");


  getStatistics().createHistogramWithAxes(new TH1D("3g_phi_lifetime", "#phi Difference", 361, -180.0, 180.0), 
                                          "#phi between annihilation hits [deg]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("DOP_3g_phi_lifetime", "#phi Difference", 361, -180.0, 180.0), 
                                          "#phi between annihilation hits [deg]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("tdiff_3g_phi_lifetime", "#phi Difference", 361, -180.0, 180.0), 
                                          "#phi between annihilation hits [deg]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("theta_3g_phi_lifetime", "#phi Difference", 361, -180.0, 180.0), 
                                          "#phi between annihilation hits [deg]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("phi_3g_phi_lifetime", "#phi Difference", 361, -180.0, 180.0), 
                                          "#phi between annihilation hits [deg]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("phi0_3g_phi_lifetime", "#phi Difference", 361, -180.0, 180.0), 
                                          "#phi between annihilation hits [deg]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("dist_3g_phi_lifetime", "#phi Difference", 361, -180.0, 180.0), 
                                          "#phi between annihilation hits [deg]", "Number of Hit Pairs");
  
  getStatistics().createHistogramWithAxes(new TH1D("vtx_3g_phi_lifetime", "#phi Difference", 361, -180.0, 180.0), 
                                          "#phi between annihilation hits [deg]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("ap_all_3g_phi_lifetime", "#phi Difference", 361, -180.0, 180.0), 
                                          "#phi between annihilation hits [deg]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(new TH1D("ap_3g_phi_lifetime", "#phi Difference", 361, -180.0, 180.0), 
                                          "#phi between annihilation hits [deg]", "Number of Hit Pairs");


  getStatistics().createHistogramWithAxes(new TH2D("3g_annihilation_point_xy_lifetime", "Annihilation Point - XY", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Y [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("DOP_3g_annihilation_point_xy_lifetime", "Annihilation Point - XY", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Y [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("tdiff_3g_annihilation_point_xy_lifetime", "Annihilation Point - XY", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Y [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("theta_3g_annihilation_point_xy_lifetime", "Annihilation Point - XY", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Y [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("phi_3g_annihilation_point_xy_lifetime", "Annihilation Point - XY", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Y [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("phi0_3g_annihilation_point_xy_lifetime", "Annihilation Point - XY", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Y [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("dist_3g_annihilation_point_xy_lifetime", "Annihilation Point - XY", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Y [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("vtx_3g_annihilation_point_xy_lifetime", "Annihilation Point - XY", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Y [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("ap_all_3g_annihilation_point_xy_lifetime", "Annihilation Point - XY", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Y [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("ap_3g_annihilation_point_xy_lifetime", "Annihilation Point - XY", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Y [cm]");


  getStatistics().createHistogramWithAxes(new TH2D("3g_annihilation_point_xz_lifetime", "Annihilation Point - XZ", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Z [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("DOP_3g_annihilation_point_xz_lifetime", "Annihilation Point - XZ", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Z [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("tdiff_3g_annihilation_point_xz_lifetime", "Annihilation Point - XZ", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Z [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("theta_3g_annihilation_point_xz_lifetime", "Annihilation Point - XZ", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Z [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("phi_3g_annihilation_point_xz_lifetime", "Annihilation Point - XZ", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Z [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("phi0_3g_annihilation_point_xz_lifetime", "Annihilation Point - XZ", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Z [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("dist_3g_annihilation_point_xz_lifetime", "Annihilation Point - XZ", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Z [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("vtx_3g_annihilation_point_xz_lifetime", "Annihilation Point - XZ", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Z [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("ap_all_3g_annihilation_point_xz_lifetime", "Annihilation Point - XZ", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Z [cm]");

  getStatistics().createHistogramWithAxes(new TH2D("ap_3g_annihilation_point_xz_lifetime", "Annihilation Point - XZ", 201, -100.0, 100.0, 201, -100.0, 100.0), 
                                          "X [cm]", "Z [cm]");


  getStatistics().createHistogramWithAxes(new TH1D("3g_tot_all_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("3g_tot_mask_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("hits_3g_tot_all_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("hits_3g_tot_mask_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("3g_tot_exactly_4_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("3g_tot_prompts_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("3g_tot_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("DOP_3g_tot_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("tdiff_3g_tot_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("theta_3g_tot_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("phi_3g_tot_lifetime",( "3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("phi0_3g_tot_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("dist_3g_tot_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("vtx_3g_tot_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("ap_3g_stats_events", "Number of events classified as o-Ps", 10, -0.5, 9.5),
                                          "Number of events", "Multiplicity");

  getStatistics().createHistogramWithAxes(new TH1D("ap_3g_stats_prompts", "Number of prompts in events classified as o-Ps", 10, -0.5, 9.5),
                                          "Number of prompts", "Multiplicity");

  getStatistics().createHistogramWithAxes(new TH1D("ap_3g_stats_annihilations", "Number of annihilation groups in events classified as o-Ps", 10, -0.5, 9.5),
                                          "Number of annihilation groups", "Multiplicity");

  getStatistics().createHistogramWithAxes(new TH1D("ap_all_3g_tot_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("ap_3g_tot_lifetime", ("3 gamma event - "+energy_hist_title).c_str(), 201, 0.0, fToTHistoUpperLimit),
                                          energy_units.c_str(), "Number of Hits");


  getStatistics().createHistogramWithAxes(new TH1D("3g_DOP_lifetime", "3 gamma event - distance between source and the plane", 100, 0.0,100.0),
                                          "Distance [cm]", "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("DOP_3g_DOP_lifetime", "3 gamma event - distance between source and the plane", 100, 0.0,100.0),
                                          "Distance [cm]", "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("theta_3g_DOP_lifetime", "3 gamma event - distance between source and the plane", 100, 0.0,100.0),
                                          "Distance [cm]", "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("phi_3g_DOP_lifetime", "3 gamma event - distance between source and the plane", 100, 0.0,100.0),
                                          "Distance [cm]", "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("phi0_3g_DOP_lifetime", "3 gamma event - distance between source and the plane", 100, 0.0,100.0),
                                          "Distance [cm]", "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("tdiff_3g_DOP_lifetime", "3 gamma event - distance between source and the plane", 100, 0.0,100.0),
                                          "Distance [cm]", "Number of Hits");
  
  getStatistics().createHistogramWithAxes(new TH1D("dist_3g_DOP_lifetime", "3 gamma event - distance between source and the plane", 100, 0.0,100.0),
                                          "Distance [cm]", "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("vtx_3g_DOP_lifetime", "3 gamma event - distance between source and the plane", 100, 0.0,100.0),
                                          "Distance [cm]", "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("ap_all_3g_DOP_lifetime", "3 gamma event - distance between source and the plane", 100, 0.0,100.0),
                                          "Distance [cm]", "Number of Hits");

  getStatistics().createHistogramWithAxes(new TH1D("ap_3g_DOP_lifetime", "3 gamma event - distance between source and the plane", 100, 0.0,100.0),
                                          "Distance [cm]", "Number of Hits");

  getStatistics().createHistogramWithAxes(
      new TH2D("3g_rel_angles_lifetime", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH2D("DOP_3g_rel_angles_lifetime", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH2D("tdiff_3g_rel_angles_lifetime", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH2D("theta_3g_rel_angles_lifetime", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH2D("phi_3g_rel_angles_lifetime", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH2D("phi0_3g_rel_angles_lifetime", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH2D("dist_3g_rel_angles_lifetime", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH2D("vtx_3g_rel_angles_lifetime", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH2D("ap_all_3g_rel_angles_lifetime", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH2D("ap_3g_rel_angles_lifetime", "Sum vs. difference of two smallest relative angles in 3 gamma event", 250, 0.0, 250, 200, 0.0, 200.0),
      "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(new TH2D("3g_rel_angles_sel_lifetime",
                                          "Sum vs. difference of two smallest relative angles in 3 gamma event - after cut", 250, 0.0, 250,
                                          200, 0.0, 200.0), "ang1+ang2 [deg]", "ang2-ang1 [deg]");

  getStatistics().createHistogramWithAxes(
      new TH1D("3g_timeDiff_lifetime", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH1D("DOP_3g_timeDiff_lifetime", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH1D("tdiff_3g_timeDiff_lifetime", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH1D("dist_3g_timeDiff_lifetime", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH1D("vtx_3g_timeDiff_lifetime", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");


  getStatistics().createHistogramWithAxes(
      new TH1D("theta_3g_timeDiff_lifetime", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH1D("phi_3g_timeDiff_lifetime", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH1D("phi0_3g_timeDiff_lifetime", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH1D("ap_all_3g_timeDiff_lifetime", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");

  getStatistics().createHistogramWithAxes(
      new TH1D("ap_3g_timeDiff_lifetime", "Time difference between hits", 201, -fMaxTimeDiff, fMaxTimeDiff),
      "Time  Difference between annihilation hits [ps]", "Number of Hit Pairs");



  getStatistics().createHistogramWithAxes(new TH1D("lifetime_3g_prompt", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow),"Time Diff [ps]", "Number of events");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_3g_prompt_zoom", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_DOP_3g_prompt", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_tdiff_3g_prompt", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_tdiff_3g_prompt_zoom", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_theta_3g_prompt", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_theta_3g_prompt_zoom", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_phi_3g_prompt", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_phi_3g_prompt_zoom", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_phi0_3g_prompt", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_phi0_3g_prompt_zoom", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 


  getStatistics().createHistogramWithAxes(new TH1D("lifetime_vtx_3g_prompt", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_vtx_3g_prompt_zoom", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_dist_3g_prompt", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_dist_3g_prompt_zoom", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_ap_all_3g_prompt", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_ap_all_3g_prompt_zoom", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_ap_3g_prompt", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_ap_3g_prompt_zoom", "Time difference of 3 gamma hits decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 

/*  getStatistics().createHistogramWithAxes(new TH1D("lifetime_1g_prompt", "Time difference of 1 gamma hits decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_1g_prompt_zoom", "Time difference of 1 gamma hits decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_tot_1g_prompt", "Time difference of 1 gamma hits decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");
                                        
  getStatistics().createHistogramWithAxes(new TH1D("lifetime_tot_1g_prompt_zoom", "Time difference of 1 gamma hits decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_tot_larger_1g_prompt", "Time difference of 1 gamma hits decay time and prompt emmission time", 201,
                                          -1.5 * fEventTimeWindow, 1.5*fEventTimeWindow), "Time Diff [ps]", "Number of events");

  getStatistics().createHistogramWithAxes(new TH1D("lifetime_tot_larger_1g_prompt_zoom", "Time difference of 1 gamma hits decay time and prompt emmission time", 201,
                                          -1*fEventTimeWindow_zoom, 10*fEventTimeWindow_zoom), "Time Diff [ps]", "Number of events"); 
*/
}

void EventCategorizer::initialiseCalibrationHistograms(bool includeTrento)
{
  auto minScinID = getParamBank().getScins().begin()->first;
  auto maxScinID = getParamBank().getScins().rbegin()->first;

  // Synchronization of TOF with annihilaion-deexcitation pairs
  getStatistics().createHistogramWithAxes(new TH2D("tdiff_anni_scin", "A-D time difference for annihilation hit per scin", maxScinID - minScinID + 1,
                                                   minScinID - 0.5, maxScinID + 0.5, 200, -fMaxTimeDiff, fMaxTimeDiff),
                                          "Scin ID", "Time diffrence [ps]");

  getStatistics().createHistogramWithAxes(new TH2D("tdiff_deex_scin", "A-D time difference for deex hit per scin", maxScinID - minScinID + 1,
                                                   minScinID - 0.5, maxScinID + 0.5, 200, -fMaxTimeDiff, fMaxTimeDiff),
                                          "Scin ID", "Time diffrence [ps]");

  // Time walk histograms
  double revToTLimit = 0.000000030;

  getStatistics().createHistogramWithAxes(
      new TH2D("time_walk_ab_tdiff", "AB TDiff vs. reversed ToT", 200, -fMaxTimeDiff / 2.0, fMaxTimeDiff / 2.0, 200, -revToTLimit, revToTLimit),
      "AB Time Difference [ps]", "Reversed ToT [1/ps]");

  getStatistics().createHistogramWithAxes(
      new TH2D("time_walk_tof", "TOF vs. reversed ToT", 200, -fMaxTimeDiff / 2.0, fMaxTimeDiff / 2.0, 200, -revToTLimit, revToTLimit),
      "Time of Flight [ps]", "Reversed ToT [1/ps]");

  getStatistics().createHistogramWithAxes(new TH3D("time_walk_ab_tdiff_scin", "AB TDiff vs. reversed ToT per Scintillator", 200, -fMaxTimeDiff / 2.0,
                                                   fMaxTimeDiff / 2.0, 200, -revToTLimit, revToTLimit, maxScinID - minScinID + 1, minScinID - 0.5,
                                                   maxScinID + 0.5),
                                          "AB Time Difference [ps]", "Reversed ToT [1/ps]", "Scintillator ID");

  getStatistics().createHistogramWithAxes(new TH3D("time_walk_tof_scin", "TOF vs. reversed ToT per Scintillator", 200, -fMaxTimeDiff / 2.0,
                                                   fMaxTimeDiff / 2.0, 200, -revToTLimit, revToTLimit, maxScinID - minScinID + 1, minScinID - 0.5,
                                                   maxScinID + 0.5),
                                          "Time of Flight [ps]", "Reversed ToT [1/ps]", "Scintillator ID");

  // Fine channel offset calibration - using only 2g events
  auto minChannelID = getParamBank().getChannels().begin()->first;
  auto maxChannelID = getParamBank().getChannels().rbegin()->first;

  getStatistics().createHistogramWithAxes(new TH2D("evtcat_channel_offsets", "Offset of Channel in Matrix vs. Channel ID in annihilation hits",
                                                   maxChannelID - minChannelID + 1, minChannelID - 0.5, maxChannelID + 0.5, 200, -3000.0, 3000.0),
                                          "Channel ID", "Offset");

  if (includeTrento)
  {
    // Cosmic ToF - histograms for Trento setup
    for (int scinID = 201; scinID <= 226; ++scinID)
    {
      getStatistics().createHistogramWithAxes(new TH2D(Form("cosmic_tof_tdiff_scin_%d_all", scinID),
                                                       Form("Time of Flight between hits from scin ID %d and from layer below", scinID), 26, 200.5,
                                                       226.5, 200, -fMaxTimeDiff, fMaxTimeDiff),
                                              "Scintillator ID", "Time of Flight [ps]");

      getStatistics().createHistogramWithAxes(new TH2D(Form("cosmic_tof_tdiff_scin_%d_cut", scinID),
                                                       Form("Time of Flight between hits from scin ID %d and from layer below", scinID), 26, 200.5,
                                                       226.5, 200, -fMaxTimeDiff, fMaxTimeDiff),
                                              "Scintillator ID", "Time of Flight [ps]");

      getStatistics().createHistogramWithAxes(new TH2D(Form("cosmic_tof_offset_scin_%d_all", scinID),
                                                       Form("Time of Flight between hits from scin ID %d and from layer below", scinID), 26, 200.5,
                                                       226.5, 200, -fMaxTimeDiff, fMaxTimeDiff),
                                              "Scintillator ID", "Time of Flight [ps]");

      getStatistics().createHistogramWithAxes(new TH2D(Form("cosmic_tof_offset_scin_%d_cut", scinID),
                                                       Form("Time of Flight between hits from scin ID %d and from layer below", scinID), 26, 200.5,
                                                       226.5, 200, -fMaxTimeDiff, fMaxTimeDiff),
                                              "Scintillator ID", "Time of Flight [ps]");
    }

    getStatistics().createHistogramWithAxes(new TH1D("cosmic_hits_x_diff_all", "X-Position difference of two comsic hits", 120, -30.0, 30.0),
                                            "positon diff [cm]", "Number of pairs");
    getStatistics().createHistogramWithAxes(
        new TH1D("cosmic_hits_x_diff_cut", "X-Position difference of two comsic hits after angle cut", 120, -30.0, 30.0), "positon diff [cm]",
        "Number of pairs");

    getStatistics().createHistogramWithAxes(new TH1D("cosmic_hits_z_diff_all", "Z-Position difference of two comsic hits", 120, -30.0, 30.0),
                                            "positon diff [cm]", "Number of pairs");
    getStatistics().createHistogramWithAxes(
        new TH1D("cosmic_hits_z_diff_cut", "Z-Position difference of two comsic hits after angle cut", 120, -30.0, 30.0), "positon diff [cm]",
        "Number of pairs");

    getStatistics().createHistogramWithAxes(new TH1D("cosmic_hits_y_diff_all", "Y-Position difference of two comsic hits", 120, -30.0, 30.0),
                                            "positon diff [cm]", "Number of pairs");
    getStatistics().createHistogramWithAxes(
        new TH1D("cosmic_hits_y_diff_cut", "Y-Position difference of two comsic hits after angle cut", 120, -30.0, 30.0), "positon diff [cm]",
        "Number of pairs");

    getStatistics().createHistogramWithAxes(new TH1D("cosmic_hits_theta_xz_all", "Theta of two comsic hits", 360, -180.0, 180.0), "theta [deg]",
                                            "Number of pairs");
    getStatistics().createHistogramWithAxes(new TH1D("cosmic_hits_theta_xz_cut", "Theta of two comsic hits", 360, -180.0, 180.0), "theta [deg]",
                                            "Number of pairs");

    getStatistics().createHistogramWithAxes(new TH1D("cosmic_hits_theta_xy_all", "Theta of two comsic hits", 360, -180.0, 180.0), "theta [deg]",
                                            "Number of pairs");
    getStatistics().createHistogramWithAxes(new TH1D("cosmic_hits_theta_xy_cut", "Theta of two comsic hits", 360, -180.0, 180.0), "theta [deg]",
                                            "Number of pairs");
  }
}
