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
 *  @file EventCategorizer.h
 */

#ifndef EVENTCATEGORIZER_H
#define EVENTCATEGORIZER_H

#include <JPetUserTask/JPetUserTask.h>
#include "EventCategorizerTools.h"
#include <JPetEvent/JPetEvent.h>
#include <JPetHit/JPetHit.h>
#include <vector>
#include <map>


class JPetWriter;

/**
 * @brief User Task categorizing Events
 *
 * Task attempts to add types of events to each event. Each category/type
 * has separate method for checking, if current event fulfills set of conditions.
 * These methods are defined in tools class. More than one type can be added to an event.
 * Set of controll histograms are created, unless the user decides not to produce them.
 */
class EventCategorizer : public JPetUserTask{
public:
	EventCategorizer(const char * name);
	virtual ~EventCategorizer();
	virtual bool init() override;
	virtual bool exec() override;
	virtual bool terminate() override;
	static const double kUnknownEventType;
protected:
	const std::string kBack2BackSlotThetaDiffParamKey = "Back2Back_Categorizer_SlotThetaDiff_float";
	const std::string kScatterTOFTimeDiffParamKey = "Scatter_Categorizer_TOF_TimeDiff_float";
	const std::string kDeexTOTCutMinParamKey = "Deex_Categorizer_TOT_Cut_Min_float";
	const std::string kDeexTOTCutMaxParamKey = "Deex_Categorizer_TOT_Cut_Max_float";
	const std::string kToTCut2AnniMinParamKey = "2Anni_TOT_Cut_Min_float";
	const std::string kToTCut2AnniMaxParamKey = "2Anni_TOT_Cut_Max_float";
	const std::string kToTCut3AnniMinParamKey = "3Anni_TOT_Cut_Min_float";
	const std::string kToTCut3AnniMaxParamKey = "3Anni_TOT_Cut_Max_float";
	const std::string kToTCutScattMaxParamKey = "Scatt_TOT_Cut_Max_float";
	const std::string kMaxTimeDiffParamKey = "EventCategorizer_MaxTimeDiff_float";
	const std::string kSaveControlHistosParamKey = "Save_Control_Histograms_bool";
    const std::string kTOTCalculationType = "HitFinder_TOTCalculationType_std::string";
    const std::string kHistoLimit = "EventCategorizer_HistoLimit_float";
    const std::string kEnergy_unit = "EventCategorizer_kEnergy_unit_std::string";
	void saveEvents(const std::vector<JPetEvent>& event);
	double fScatterTOFTimeDiff = 2000.0;
	double fB2BSlotThetaDiff = 3.0;
	double fDeexTOTCutMin = 67000.0;
	double fDeexTOTCutMax = 122000.0;
    double fToTCut2AnniMin = 5000000.0;
    double fToTCut2AnniMax = 7500000.0;
    double fToTCut3AnniMin = 1000.0;
    double fToTCut3AnniMax = 6700.0;
	double fToTCutScattMax = 80000.0;
	double fMaxTimeDiff = 1000.;
	double fHistoLimit = 2e5;
	bool fSaveControlHistos = true;
    std::string fTOTCalculationType = "";
	std::string fEnergy_unit = "";
	void initialiseHistograms();


    std::vector<float> fPosX;
    std::vector<float> fPosY;
    std::vector<float> fPosZ;
    std::vector<float> fEnergy;
    std::vector<float> fTime;
	std::vector<unsigned int> fHitType;
    std::vector<unsigned int> fVtxIndex;
    bool fIsAcc = kFALSE;
	bool fIsMC = false;
	bool fIsOPs = kFALSE;
	bool fIsPickOff = kFALSE;
	bool fContainsPrompt = kFALSE;
	bool fIsScattered = kFALSE;
	bool fIsSecondary = kFALSE;
	int n_oPs_true = 0;
	int n_oPs_scatter = 0;
	int n_random = 0;
	int n_oPs_prompt_true = 0;
	int n_scattered = 0;
	int n_prompt_true = 0;
	int n_prompt_scatter = 0;
	int n_oPs_prompt_scatter = 0;
	};
#endif /* !EVENTCATEGORIZER_H */
