/**
 *  @copyright Copyright 2021 The J-PET Framework Authors. All rights reserved.
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
 *  @file EventAnalyzer.h
 */

#ifndef EVENTANALYZER_H
#define EVENTANALYZER_H

#include <JPetEvent/JPetEvent.h>
#include <JPetUserTask/JPetUserTask.h>
#include "../ModularDetectorAnalysis/EventCategorizerTools.h"

class EventAnalyzer : public JPetUserTask
{
public:
  EventAnalyzer(const char* name);
  virtual ~EventAnalyzer();
  virtual bool init() override;
  virtual bool exec() override;
  virtual bool terminate() override;

  double fEventTimeWindow = 5000.0;
  double fScatterTOFTimeDiff = 2000.0;
  double fScatterTimeMin = -5000.0;
  double fScatterTimeMax = 0.0;
  double fScatterAngleMin = 160.0;
  double fScatterAngleMax = 180.0;
  double fMaxTimeDiff = 15000.0;
  double f2gThetaDiff = 3.0;
  double f2gTimeDiff = 2000.0;
  double f3gMinRelAngle = 185.0;
  double fToTCutAnniMin = 150000.0;
  double fToTCutAnniMax = 250000.0;
  double fToTCutDeexMin = 270000.0;
  double fToTCutDeexMax = 370000.0;
  double fToTHistoUpperLimit = 200000.0;
  double fLORAngleCut = 5.0;
  double fLORPosZCut = 5.0;
  double fSourceDistXYCut = 5.0;
  double fSourceDistZCut = 10.0;
  double fDetectorYRotationDeg = 60.0;
  double fCosmicMaxThetaDiffDeg = 3.0;
  TVector3 fSourcePos;
  EventCategorizerTools::ScatterTestType fTestType = EventCategorizerTools::kSimpleParam;

  bool fSaveControlHistos = true;
  bool fSaveCalibHistos = false;
  bool fTrentoCalibHistos = false;

  static int n_acc;
  static int n_all;

protected:
  void fillResolutionHistograms(const JPetEvent& event, const JPetTimeWindowMC* tw);
  bool fIsMC = false;
};
#endif /* !EVENTANALYZER_H */
