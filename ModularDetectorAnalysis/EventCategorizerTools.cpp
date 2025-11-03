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
 *  @file EventCategorizerTools.cpp
 */

#include "EventCategorizerTools.h"
#include <Hits/JPetPhysRecoHit/JPetPhysRecoHit.h>
#include <Hits/JPetMCRecoHit/JPetMCRecoHit.h>
#include <Math/DistFunc.h>
#include <TMath.h>
#include <TRandom.h>
#include <vector>
#include <algorithm> 
#include <random> 

using namespace std;

void EventCategorizerTools::selectAndCategorizeFourHits(const JPetEvent& event, const JPetTimeWindowMC* timeWindowMC,
  JPetStatistics& stats, double maxZ, double totCutAnniMin, double totCutAnniMax, double totCutDeexMin,
  double totCutDeexMax)
{
  const auto& allHits = event.getHits();

  // Select hits inside maxZ
  std::vector<const JPetMCRecoHit*> selectedHits;
  for (const auto* hit : allHits) {
    if (std::fabs(hit->getPosZ()) < maxZ) {
      selectedHits.push_back(dynamic_cast<const JPetMCRecoHit*>(hit));
    }
  }

  if (selectedHits.size() != 4) return; // Need exactly 4 hits

  // Classify hits by ToT into ops and prompt, also get raw hits for tag and vertex info
  std::vector<const JPetMCRecoHit*> selectedOpsHits;
  std::vector<const JPetMCRecoHit*> selectedPromptHits;
  std::vector<const JPetRawMCHit*> selectedOpsRawHits;
  std::vector<const JPetRawMCHit*> selectedPromptRawHits;

  for (const auto& hit : selectedHits) {
    double tot = calculateToT(hit);
    auto& rawHit = timeWindowMC->getMCHit<JPetRawMCHit>(hit->getMCindex());
    if (tot >= totCutAnniMin && tot <= totCutAnniMax) {
      selectedOpsHits.push_back(hit);
      selectedOpsRawHits.push_back(&rawHit);
    } else if (tot >= totCutDeexMin && tot <= totCutDeexMax) {
      selectedPromptHits.push_back(hit);
      selectedPromptRawHits.push_back(&rawHit);
    }
  }

  if (selectedOpsHits.size() != 3 || selectedPromptHits.size() != 1) return;

  // Check gamma tags:
  bool validOpsTags = true;
  for (const auto* rawHit : selectedOpsRawHits) {
    //std::cout<<rawHit->getGammaTag()<<", ";
    if (rawHit->getGammaTag() != 3) {
      validOpsTags = false;
      //break;
    }
  }
  
  bool validPromptTag = (selectedPromptRawHits[0]->getGammaTag() == 1);
  //std::cout<<selectedPromptRawHits[0]->getGammaTag()<<std::endl;

  if (!validOpsTags || !validPromptTag) {
    // Background case
    fillBackground(stats, selectedOpsHits, selectedPromptHits, timeWindowMC);
    //return;
  }

  // Shuffle ops hits for randomness
  std::random_device rd;
  std::default_random_engine rng(rd());
  std::shuffle(selectedOpsHits.begin(), selectedOpsHits.end(), rng);

  const auto* firstHit  = selectedOpsHits[0];
  const auto* secondHit = selectedOpsHits[1];
  const auto* thirdHit  = selectedOpsHits[2];

  // Check all ops and prompt hits have the same vertex (using raw hits)
  int commonVertex = selectedOpsRawHits[0]->getMCVtxIndex();
  bool sameVertex = true;
  for (const auto* rawHit : selectedOpsRawHits) {
    if (rawHit->getMCVtxIndex() != commonVertex) {
      sameVertex = false;
      break;
    }
  }
  int promptVertex = selectedPromptRawHits[0]->getMCVtxIndex();

  //std::cout<<"Same vertex: "<<selectedOpsRawHits[0]->getMCVtxIndex()<<", "<<selectedOpsRawHits[1]->getMCVtxIndex()<<", "<<selectedOpsRawHits[2]->getMCVtxIndex()<<", "<<promptVertex<<std::endl;

  if (!sameVertex || promptVertex != commonVertex) {
    // Background if vertices don't match
    fillBackground(stats, selectedOpsHits, selectedPromptHits, timeWindowMC);
    return;
  }

  // Calculate time differences between ops hits
  double time_hit1 = firstHit->getTime();
  double time_hit2 = secondHit->getTime();
  double time_hit3 = thirdHit->getTime();
  std::vector<double> time_Diffs = {time_hit1 - time_hit2, time_hit2 - time_hit3, time_hit3 - time_hit1};

  // Calculate relative angles
  std::vector<double> relativeAngles;
  relativeAngles.push_back(TMath::RadToDeg() * firstHit->getPos().Angle(secondHit->getPos()));
  relativeAngles.push_back(TMath::RadToDeg() * secondHit->getPos().Angle(thirdHit->getPos()));
  relativeAngles.push_back(TMath::RadToDeg() * thirdHit->getPos().Angle(firstHit->getPos()));
  std::sort(relativeAngles.begin(), relativeAngles.end());

  double transformedX = relativeAngles[0] + relativeAngles[1];
  double transformedY = relativeAngles[1] - relativeAngles[0];

  // Fill signal histograms
  for (const auto& hit : selectedOpsHits) {
    stats.fillHistogram("ap_3g_tot_oPs", EnergyToToT(hit->getEnergy()));
  }
  for (double timeDiff : time_Diffs) {
    stats.fillHistogram("ap_3g_timeDiff_oPs", timeDiff);
  }
  stats.fillHistogram("ap_3g_rel_angles_oPs", transformedX, transformedY);
  stats.fillHistogram("ap_3g_tot_prompt", EnergyToToT(selectedPromptHits[0]->getEnergy()));
  stats.fillHistogram("event_count", 1);
}


void EventCategorizerTools::fillBackground(JPetStatistics& stats, const std::vector<const JPetMCRecoHit*>& selectedOpsHits, const std::vector<const JPetMCRecoHit*>& selectedPromptHits,
  const JPetTimeWindowMC* timeWindowMC)
{
  // Fill background histograms and count
  stats.fillHistogram("event_count", 2);
  for (const auto& hit : selectedOpsHits) {
    stats.fillHistogram("ap_3g_tot_bkg_oPs", EnergyToToT(hit->getEnergy()));
  }
  for (const auto& hit : selectedPromptHits) {
    stats.fillHistogram("ap_3g_tot_bkg_prompt", EnergyToToT(hit->getEnergy()));
  }
  // You can add timeDiffs and relative angles if you want, but they require hit times and positions
}

void EventCategorizerTools::processMCEvent(const JPetEvent& event, const JPetTimeWindowMC* timeWindowMC, const TVector3& sourcePos, JPetStatistics& stats) {
  std::vector<const JPetMCRecoHit*> oPsHits;
  std::vector<const JPetMCRecoHit*> promptHits;
  std::vector<JPetRawMCHit> oPsRawHits;
  std::vector<JPetRawMCHit> promptRawHits;
  std::vector<int> fVtxValues;      // vertices from oPs hits (for triple check)
  std::vector<int> fTrueVtxValues;  // vertices confirmed by triple coincidence

  // Process all hits, fill oPsHits, promptHits, and vertex IDs
  processHits(event, timeWindowMC, stats, oPsHits, promptHits, oPsRawHits, promptRawHits, fVtxValues);

  // Identify vertices with at least 3 hits
  checkTripleVertex(fVtxValues, fTrueVtxValues);

  // Collect all hits and raw hits for random classification
  std::vector<const JPetMCRecoHit*> allHits;
  std::vector<JPetRawMCHit> allRawHits;
  for (uint i = 0; i < event.getHits().size(); i++) {
      auto hit = dynamic_cast<const JPetMCRecoHit*>(event.getHits().at(i));
      if (!hit) continue;
      allHits.push_back(hit);
      allRawHits.push_back(timeWindowMC->getMCHit<JPetRawMCHit>(hit->getMCindex()));
  }

  // Separate random hits (vertex not in true vertices)
  std::vector<const JPetMCRecoHit*> randomHits;
  std::vector<JPetRawMCHit> randomRawHits;

  for (size_t i = 0; i < allHits.size(); i++) {
      const auto& rawHit = allRawHits[i];
      if (std::find(fTrueVtxValues.begin(), fTrueVtxValues.end(), rawHit.getMCVtxIndex()) == fTrueVtxValues.end()) {
          // Hit with vertex ID NOT in true vertices → random
          randomHits.push_back(allHits[i]);
          randomRawHits.push_back(rawHit);
          stats.fillHistogram("random_tot_true", calculateToT(allHits[i]));
      }
  }

  // Fill histograms for true categorized hits (oPs, prompt)
  fillHistograms(stats, oPsHits, oPsRawHits, promptHits, promptRawHits, fTrueVtxValues, sourcePos);

}

// Extracts hit properties and categorizes them
void EventCategorizerTools::processHits(const JPetEvent& event, const JPetTimeWindowMC* timeWindowMC, JPetStatistics& stats, 
                                      std::vector<const JPetMCRecoHit*>& oPsHits, std::vector<const JPetMCRecoHit*>& promptHits, 
                                      std::vector<JPetRawMCHit>& oPsRawHits, std::vector<JPetRawMCHit>& promptRawHits, std::vector<int>& fVtxValues) 
{
  for (uint i = 0; i < event.getHits().size(); i++) {
      auto hit = dynamic_cast<const JPetMCRecoHit*>(event.getHits().at(i));
      if (!hit) continue;

      stats.fillHistogram("all_tot_true", calculateToT(hit));

      const JPetRawMCHit mcHit = timeWindowMC->getMCHit<JPetRawMCHit>(hit->getMCindex());
      bool fIsSecondary = false, fIsScattered = false, fContainsPrompt = false, fIsPickOff = false, fIsOPs = false;
      
      categorizeHit(mcHit, fIsSecondary, fIsScattered, fContainsPrompt, fIsPickOff, fIsOPs);

      if (fIsOPs && !fIsScattered && !fIsSecondary) {
          oPsHits.push_back(hit);
          oPsRawHits.push_back(mcHit);
          fVtxValues.push_back(mcHit.getMCVtxIndex());
      }
      if (fContainsPrompt && !fIsScattered && !fIsSecondary) {
          promptHits.push_back(hit);
          promptRawHits.push_back(mcHit);
      }


      if(fIsScattered || fIsSecondary) stats.fillHistogram("scatter_tot_true", calculateToT(hit));
      if(fContainsPrompt && !fIsScattered && !fIsSecondary) stats.fillHistogram("prompt_tot_true", calculateToT(hit));
      if(fIsPickOff) stats.fillHistogram("2g_tot_true", calculateToT(hit));
  }
}

// Determines hit category based on hit type
void EventCategorizerTools::categorizeHit(const JPetRawMCHit& mcHit, bool& fIsSecondary, bool& fIsScattered, bool& fContainsPrompt, 
                                          bool& fIsPickOff, bool& fIsOPs) 
{
  int hitType = mcHit.getGammaTag();
  int hitNumber = hitType % 10;
  int n = 1;

  if (hitType % 100 != 0 && hitNumber != 0 && hitType > 10) fIsScattered = true;
  if (hitNumber == 1) fContainsPrompt = true;
  else if (hitNumber == 2) fIsPickOff = true;
  else if (hitNumber == 3) fIsOPs = true;
  else if (hitNumber == 0) {
      fIsSecondary = true;
      while (hitNumber == 0 && hitType > 0) {
          n *= 10;
          hitNumber = (hitType / n) % 10;
          if (hitType % (n * 100) == 0 && hitNumber != 0 && hitType / n > 10)
              fIsScattered = true;
      }
  }
}

void EventCategorizerTools::checkTripleVertex(const std::vector<int>& fVtxValues, std::vector<int>& fTrueVtxValues) {
  if (fVtxValues.size() < 3) return;

  std::map<int, int> vtxCount;
  for (const auto& vtx : fVtxValues) vtxCount[vtx]++;
  
  for (const auto& entry : vtxCount) {
      if (entry.second >= 3) fTrueVtxValues.push_back((int)entry.first);
  }
}

// Fills histograms based on categorized hits
void EventCategorizerTools::fillHistograms(JPetStatistics& stats, const std::vector<const JPetMCRecoHit*>& oPsHits, const std::vector<JPetRawMCHit>& oPsRawHits, 
                                          const std::vector<const JPetMCRecoHit*>& promptHits, const std::vector<JPetRawMCHit>& promptRawHits, std::vector<int>& fTrueVtxValues, const TVector3& sourcePos) 
{
  bool vtxCheck, promptVtxCheck = false;
  int n_ops = 0;
  std::vector<int> true_oPsHits;      // Store indices of true oPs
  std::vector<int> true_promptHits;   // Store indices of true prompts
  for (int i = 0; i<oPsHits.size(); i++) {
      vtxCheck = false;
      promptVtxCheck = false;
      const auto& oPs = oPsHits[i];
      const auto& oPsRaw = oPsRawHits[i];
      stats.fillHistogram("3g_tot_true", calculateToT(oPs));
      if(checkVertex(oPsRaw, fTrueVtxValues)) vtxCheck = true;
      for (size_t j = 0; j < promptRawHits.size(); j++) {
        if (checkVertex(promptRawHits[j], {oPsRawHits[i].getMCVtxIndex()})) {
            promptVtxCheck = true;
            true_promptHits.push_back(j); // Store index of true prompt hit
        }
    }

      if(vtxCheck&&promptVtxCheck){
        n_ops++;
        true_oPsHits.push_back(i);
        stats.fillHistogram("ap_3g_tot_true", calculateToT(oPs));
      }
  }

  if (true_oPsHits.size() < 3) return; 

  std::vector<const JPetMCRecoHit*> sortedHits = {
    oPsHits[true_oPsHits[0]],
    oPsHits[true_oPsHits[1]],
    oPsHits[true_oPsHits[2]]
  };

  /*std::sort(sortedHits.begin(), sortedHits.end(), 
  [](const JPetMCRecoHit* a, const JPetMCRecoHit* b) {
      return calculateToT(a) > calculateToT(b);
  });*/

  std::random_device rd;
  std::default_random_engine rng(rd());
  std::shuffle(sortedHits.begin(), sortedHits.end(), rng);
  
  // Assign sorted hits
  const auto* firstHit  = sortedHits[0];
  const auto* secondHit = sortedHits[1];
  const auto* thirdHit  = sortedHits[2];


  double time_hit1 = firstHit->getTime();
  double time_hit2 = secondHit->getTime();
  double time_hit3 = thirdHit->getTime();

  std::vector<double> time_Diffs = {time_hit1 - time_hit2, time_hit2 - time_hit3, time_hit3 - time_hit1};

  vector<double> relativeAngles;
  relativeAngles.push_back(TMath::RadToDeg() * firstHit->getPos().Angle(secondHit->getPos()));
  relativeAngles.push_back(TMath::RadToDeg() * secondHit->getPos().Angle(thirdHit->getPos()));
  relativeAngles.push_back(TMath::RadToDeg() * thirdHit->getPos().Angle(firstHit->getPos()));
  sort(relativeAngles.begin(), relativeAngles.end());

  double transformedX = relativeAngles.at(1) + relativeAngles.at(0);
  double transformedY = relativeAngles.at(1) - relativeAngles.at(0);

  double DOP = calculatePlanePointDistance(firstHit, secondHit, thirdHit, sourcePos);

  for(double timeDiff : time_Diffs) stats.fillHistogram("ap_3g_timeDiff_true", timeDiff);
  stats.fillHistogram("ap_3g_rel_angles_true", transformedX, transformedY);

  if (!oPsHits.empty()) {
      stats.fillHistogram("3g_stats_multi_annihilations_true", oPsHits.size());
      if(n_ops>0) stats.fillHistogram("ap_3g_stats_multi_annihilations_true", n_ops);
      stats.fillHistogram("3g_stats_multi_prompts_true", promptHits.size());
  }
}

bool EventCategorizerTools::checkVertex(const JPetRawMCHit& rawHit, const std::vector<int> fTrueVtxValues){
  return std::find(fTrueVtxValues.begin(), fTrueVtxValues.end(), rawHit.getMCVtxIndex()) != fTrueVtxValues.end();
}


// bool EventCategorizerTools::checkFor1Gamma(const JPetEvent& event, const TVector3& sourcePos, double totCutAnniMin, double totCutAnniMax, double totCutAnniMin_larger, 
//                                           double totCutAnniMax_larger, double totCutDeexMin, double totCutDeexMax, JPetStatistics& stats, bool saveHistos)
// {
//   if (event.getHits().size() != 2)
//   {
//     return false;
//   }
  
//   vector<const JPetPhysRecoHit*> prompts;
//   vector<vector<const JPetPhysRecoHit*>> annihilations;

//   // First check if any of the hits in the event is prompt based on TOT selection
//   for (uint i = 0; i < event.getHits().size(); i++)
//   {
//     bool totCut = false, totCut_larger = false;
//     auto promptHit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(i));
//     if (saveHistos)
//     {
//       double hitTot = promptHit->getToT();
//       stats.fillHistogram("1g_tot", hitTot);
//       if(checkToT(promptHit, totCutAnniMin, totCutAnniMax))
//       {
//         totCut = true;
        
//         stats.fillHistogram("1g_tot_tot", hitTot);
//       }
//       if(checkToT(promptHit, totCutAnniMin_larger, totCutAnniMax_larger))
//       {
//         totCut_larger = true;
//         stats.fillHistogram("1g_tot_larger_tot", hitTot);
//       }
//     }
//     if (checkToT(promptHit, totCutDeexMin, totCutDeexMax))
//     {
//       prompts.push_back(promptHit);
//       if (saveHistos)
//       {
//         stats.fillHistogram("1g_prompt_tot", promptHit->getToT());
//       }
//     }
//   }
//   if (prompts.size() != 1)
//   {
//     return false;
//   }

//   for (uint i = 0; i < event.getHits().size(); i++)
//   {
//     auto hit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(i));
//     if (!hit)
//     {
//       continue;
//     }

//     vector<const JPetPhysRecoHit*> annih_temp = {hit};
//     annihilations.push_back(annih_temp);
//   }

//   if (annihilations.size() == 0)
//   {
//     return false;
//   }

//   double is1Gamma = false;

//   // Iterating over all combinations of found pairs and prompt photons
//   for (auto pair1g : annihilations)
//   {
//     for (auto prompt : prompts)
//     {
//       // Caculate event times - annihilation
//       double annihTime = pair1g[0]->getTime() - calculateTOF(pair1g[0], sourcePos);
//       double promptTime = prompt->getTime() - calculateTOF(prompt, sourcePos);
//       double tot_annih = pair1g[0]->getToT();
//       double tot_prompt = prompt->getToT();
  
//       // Calculate lifetime
//       double lifetime = annihTime - promptTime;

//       if (saveHistos)
//       {
//         stats.fillHistogram("lifetime_1g_prompt", lifetime);
//         stats.fillHistogram("lifetime_1g_prompt_zoom", lifetime);
//         if(checkToT(pair1g[0], totCutAnniMin, totCutAnniMax)){ 
//           is1Gamma = true;
//           stats.fillHistogram("lifetime_tot_1g_prompt", lifetime);
//           stats.fillHistogram("lifetime_tot_1g_prompt_zoom", lifetime);
//         }
//         if(checkToT(pair1g[0], totCutAnniMin_larger, totCutAnniMax_larger)) {
//           stats.fillHistogram("lifetime_tot_larger_1g_prompt", lifetime);
//           stats.fillHistogram("lifetime_tot_larger_1g_prompt_zoom", lifetime);
//         }
//       }
        
//     }
//   }
  
//   return is1Gamma;
// }

/**
 * Method for determining type of event - back to back 2 gamma
 */
bool EventCategorizerTools::checkFor2Gamma(const JPetEvent& event, JPetStatistics& stats, bool saveHistos, double maxThetaDiff, double maxTimeDiff,
                                           double totCutAnniMin, double totCutAnniMax, const TVector3& sourcePos, ScatterTestType testType,
                                           double scatterTestValue, double scatterTimeMin, double scatterTimeMax, double scatterAngleMin,
                                           double scatterAngleMax)
{
  bool isEvent2Gamma = false;
  if (event.getHits().size() != 2)
  {
    return isEvent2Gamma;
  }

  for (uint i = 0; i < event.getHits().size(); i++)
  {
    auto firstHit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(i));
    if (!firstHit)
    {
      continue;
    }

    for (uint j = i + 1; j < event.getHits().size(); j++)
    {
      auto secondHit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(j));
      if (!secondHit)
      {
        continue;
      }

      // Change order or hits, if needed
      if (secondHit->getToT() > firstHit->getToT())
      {
        firstHit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(j));
        secondHit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(i));
      }

      // Skip if scatter
      bool isScatter = checkForScatter(firstHit, secondHit, stats, false, testType, scatterTestValue, scatterTimeMin, scatterTimeMax, scatterAngleMin,
                          scatterAngleMax);

      if (checkFor2Gamma(firstHit, secondHit, stats, saveHistos, maxThetaDiff, maxTimeDiff, totCutAnniMin, totCutAnniMax, sourcePos, isScatter))
      {
        isEvent2Gamma = true;
      }
    }
  }
  return isEvent2Gamma;
}

/**
 * Method for determining type of two hits - back to back 2 gamma
 */
bool EventCategorizerTools::checkFor2Gamma(const JPetPhysRecoHit* firstHit, const JPetPhysRecoHit* secondHit, JPetStatistics& stats, bool saveHistos,
                                           double maxThetaDiff, double maxTimeDiff, double totCutAnniMin, double totCutAnniMax,
                                           const TVector3& sourcePos, bool isScatter)
{
  int scin1ID = firstHit->getScin().getID();
  int scin2ID = secondHit->getScin().getID();

  TVector3 firstVec = firstHit->getPos() - sourcePos;
  TVector3 secondVec = secondHit->getPos() - sourcePos;
  double theta = TMath::RadToDeg() * firstVec.Angle(secondVec);

  // Registration time difference, always positive
  double timeDiff = secondHit->getTime() - firstHit->getTime();
  double dist = calculateDistance(firstHit, secondHit);

  auto tot1 = firstHit->getToT();
  auto tot2 = secondHit->getToT();

  // TOF calculated by convention
  double tof = calculateTOFByConvention(firstHit, secondHit);

  // LOR angle
  TVector3 vechit1_2D(firstHit->getPosX() - sourcePos.X(), 0.0, firstHit->getPosZ() - sourcePos.Z());
  TVector3 vechit2_2D(secondHit->getPosX() - sourcePos.X(), 0.0, secondHit->getPosZ() - sourcePos.Z());
  TVector3 vechit1_1D(firstHit->getPosX() - sourcePos.X(), 0.0, 0.0);
  TVector3 vechit2_1D(secondHit->getPosX() - sourcePos.X(), 0.0, 0.0);


  TVector3 annhilationPoint = calculateAnnihilationPoint(firstHit, secondHit);

  // Pre-cuts histograms
  if (saveHistos)
  {
    stats.fillHistogram("2g_tot", tot1);
    stats.fillHistogram("2g_tot", tot2);
    stats.fillHistogram("2g_theta", theta);
    stats.fillHistogram("2g_timeDiff", timeDiff);
    stats.fillHistogram("2g_scatter_test_time", timeDiff - dist/kLightVelocity_cm_ps);
    stats.fillHistogram("2g_scatter_test_dist", timeDiff*kLightVelocity_cm_ps - dist);

  }

  if(isScatter){
    stats.fillHistogram("scatter_2g_scatter_test_time", timeDiff - dist/kLightVelocity_cm_ps);
    stats.fillHistogram("scatter_2g_scatter_test_dist", timeDiff*kLightVelocity_cm_ps - dist);
  }

  // Checking selection conditions
  bool thetaCut = checkRelativeAngles(firstHit->getPos(), secondHit->getPos(), maxThetaDiff);
  bool tDiffCut = false, totCut = false;

  if (thetaCut)
  {
    if (saveHistos)
    {
      stats.fillHistogram("theta_2g_tot", tot1);
      stats.fillHistogram("theta_2g_tot", tot2);
      stats.fillHistogram("theta_2g_theta", theta);
      stats.fillHistogram("theta_2g_timeDiff", timeDiff);
    }
  }
  // Time difference cut
  if (fabs(timeDiff) < maxTimeDiff)
  {
    tDiffCut = true;
    if (saveHistos)
    {
      stats.fillHistogram("tdiff_2g_tot", tot1);
      stats.fillHistogram("tdiff_2g_tot", tot2);
      stats.fillHistogram("tdiff_2g_theta", theta);
      stats.fillHistogram("tdiff_2g_timeDiff", timeDiff);
    }
  }
  // ToT cut
  if (tot1 > totCutAnniMin && tot1 < totCutAnniMax && tot2 > totCutAnniMin && tot2 < totCutAnniMax)
  {
    totCut = true;
    if (saveHistos)
    {
      stats.fillHistogram("tot_2g_tot", tot1);
      stats.fillHistogram("tot_2g_tot", tot2);
      stats.fillHistogram("tot_2g_theta", theta);
      stats.fillHistogram("tot_2g_timeDiff", timeDiff);
    }
  }

  //if (checkForScatter(firstHit, secondHit, stats, true, testType, scatterTestValue, scatterTimeMin, scatterTimeMax, scatterAngleMin,
  //                    scatterAngleMax)){
  //  if (saveHistos)
  //  {
  //    stats.fillHistogram("scatter_2g_tot", tot1);
  //    stats.fillHistogram("scatter_2g_tot", tot2);
  //    stats.fillHistogram("scatter_2g_theta", theta);
  //    stats.fillHistogram("scatter_2g_timeDiff", timeDiff);  
  //    stats.fillHistogram("scatter_scatter_test_time", timeDiff - dist/kLightVelocity_cm_ps);
  //    stats.fillHistogram("scatter_scatter_test_time", timeDiff*kLightVelocity_cm_ps - dist);
  //  }               
  //}


  // Pair of hits that meet cut conditions are treated as coming from annihilation point
  // Returning event as 2 gamma if meets cut conditions
  if (totCut && tDiffCut && thetaCut)
  {
    if (saveHistos)
    {

      stats.fillHistogram("ap_2g_tot", tot1);
      stats.fillHistogram("ap_2g_tot", tot2);
      stats.fillHistogram("ap_2g_theta", theta);
      stats.fillHistogram("ap_2g_timeDiff", timeDiff);

      stats.fillHistogram("ap_xy", annhilationPoint.X(), annhilationPoint.Y());
      stats.fillHistogram("ap_zx", annhilationPoint.Z(), annhilationPoint.X());
      stats.fillHistogram("ap_zy", annhilationPoint.Z(), annhilationPoint.Y());
      stats.fillHistogram("ap_pos", annhilationPoint.Z(), annhilationPoint.X(), annhilationPoint.Y());
      stats.fillHistogram("ap_xy_zoom", annhilationPoint.X(), annhilationPoint.Y());
      stats.fillHistogram("ap_zx_zoom", annhilationPoint.Z(), annhilationPoint.X());
      stats.fillHistogram("ap_zy_zoom", annhilationPoint.Z(), annhilationPoint.Y());
      stats.fillHistogram("ap_pos_zoom", annhilationPoint.Z(), annhilationPoint.X(), annhilationPoint.Y());
    }
    return true;
  }
  return false;
}

/**
 * Method for determining type of event - 3Gamma
 */
 
bool EventCategorizerTools::checkFor3Gamma(const JPetEvent& event, double minRelAngleCut, double maxTimeDiff, double totCutAnniMin, double totCutAnniMax, JPetStatistics& stats, bool saveHistos)
{
  if (event.getHits().size() != 3)
  {
    return false;
  }
  
  double is3Gamma = false;
  
  // Iteration over the hits in the event
  for (uint i = 0; i < event.getHits().size(); i++)
  {
    for (uint j = i + 1; j < event.getHits().size(); j++)
    {
      for (uint k = j + 1; k < event.getHits().size(); k++)
      {
        auto firstHit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(i));
        auto secondHit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(j));
        auto thirdHit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(k));

        vector<double> relativeAngles;
        relativeAngles.push_back(TMath::RadToDeg() * firstHit->getPos().Angle(secondHit->getPos()));
        relativeAngles.push_back(TMath::RadToDeg() * secondHit->getPos().Angle(thirdHit->getPos()));
        relativeAngles.push_back(TMath::RadToDeg() * thirdHit->getPos().Angle(firstHit->getPos()));
        sort(relativeAngles.begin(), relativeAngles.end());

        double transformedX = relativeAngles.at(1) + relativeAngles.at(0);
        double transformedY = relativeAngles.at(1) - relativeAngles.at(0);

        vector<double> timeDiffs = {firstHit->getTime() - secondHit->getTime(), secondHit->getTime() - thirdHit->getTime(), thirdHit->getTime() - secondHit->getTime()};
        double timeDiff = *max_element(timeDiffs.begin(), timeDiffs.end(), 
                                    [](double a, double b){ return std::fabs(a) < std::fabs(b);});

        if (saveHistos)
        {
          stats.fillHistogram("3g_tot", firstHit->getToT());
          stats.fillHistogram("3g_tot", secondHit->getToT());
          stats.fillHistogram("3g_tot", thirdHit->getToT());
          for(int time_ind = 0; time_ind < 3; time_ind++) stats.fillHistogram("3g_timeDiff", timeDiffs[time_ind]);
	        stats.fillHistogram("3g_rel_angles", transformedX, transformedY);
        }
        bool totCut = false, thetaCut = false, tDiffCut = false;
        if(checkToT(firstHit, totCutAnniMin, totCutAnniMax) && checkToT(secondHit, totCutAnniMin, totCutAnniMax) &&
        checkToT(thirdHit, totCutAnniMin, totCutAnniMax))
        {
          totCut = true;
          if (saveHistos)
          {
            stats.fillHistogram("tot_3g_tot", firstHit->getToT());
            stats.fillHistogram("tot_3g_tot", secondHit->getToT());
            stats.fillHistogram("tot_3g_tot", thirdHit->getToT());
            for(int time_ind = 0; time_ind < 3; time_ind++) stats.fillHistogram("tot_3g_timeDiff", timeDiffs[time_ind]);
            stats.fillHistogram("tot_3g_rel_angles", transformedX, transformedY);
          }
        }
        if(transformedX > minRelAngleCut & transformedY < 80)
        {
          thetaCut = true;
          if (saveHistos)
          {
            stats.fillHistogram("theta_3g_tot", firstHit->getToT());
            stats.fillHistogram("theta_3g_tot", secondHit->getToT());
            stats.fillHistogram("theta_3g_tot", thirdHit->getToT());
            for(int time_ind = 0; time_ind < 3; time_ind++) stats.fillHistogram("theta_3g_timeDiff", timeDiffs[time_ind]);
            stats.fillHistogram("theta_3g_rel_angles", transformedX, transformedY);
          }
        }

        if(fabs(timeDiff) < maxTimeDiff)
        {
          tDiffCut = true;
          if (saveHistos)
          {
            stats.fillHistogram("tdiff_3g_tot", firstHit->getToT());
            stats.fillHistogram("tdiff_3g_tot", secondHit->getToT());
            stats.fillHistogram("tdiff_3g_tot", thirdHit->getToT());
            for(int time_ind = 0; time_ind < 3; time_ind++) stats.fillHistogram("tdiff_3g_timeDiff", timeDiffs[time_ind]);
            stats.fillHistogram("tdiff_3g_rel_angles", transformedX, transformedY);
          }
        }
      if(totCut&&thetaCut&&tDiffCut){
        is3Gamma = true;
        if (saveHistos)
          {
            stats.fillHistogram("ap_3g_tot", firstHit->getToT());
            stats.fillHistogram("ap_3g_tot", secondHit->getToT());
            stats.fillHistogram("ap_3g_tot", thirdHit->getToT());
            for(int time_ind = 0; time_ind < 3; time_ind++) stats.fillHistogram("ap_3g_timeDiff", timeDiffs[time_ind]);
            stats.fillHistogram("ap_3g_rel_angles", transformedX, transformedY);
	        }
      }
      }
    }
  }
  return is3Gamma;
}

bool EventCategorizerTools::checkFor2GammaLifetime(const JPetEvent& event, std::vector<int> bad_ID, JPetStatistics& stats, bool saveHistos, double maxZ, double maxThetaDiff, 
                                                  double maxTimeDiff, double maxDOP, double totCutAnniMin, double totCutAnniMax, double totCutDeexMin, double totCutDeexMax, 
                                                  const TVector3& sourcePos, ScatterTestType testType, double scatterTestValue, double scatterTimeMin, double scatterTimeMax, 
                                                  double scatterAngleMin, double scatterAngleMax) {

  // Step 1: Initial checks and histogram filling for hits
  stats.fillHistogram("multiplicity", event.getHits().size());
  fillHitHistograms(event, stats, saveHistos, "none_2g_tot", "none_2g_z", "none_2g_z_ID");
  if (event.getHits().size() < 2) return false;
  fillHitHistograms(event, stats, saveHistos, "hits_2g_tot", "hits_2g_z", "hits_2g_z_ID");
  // Step 2: Identify prompt and annihilation hits
  std::vector<const JPetPhysRecoHit*> prompts;
  std::vector<std::pair<const JPetPhysRecoHit*, const JPetPhysRecoHit*>> annihilations;

  std::vector<const JPetMCRecoHit*> promptsMC;
  std::vector<std::pair<const JPetMCRecoHit*, const JPetMCRecoHit*>> annihilationsMC;

  identifyAnnihilationHits(event, totCutAnniMin, totCutAnniMax, annihilations, annihilationsMC);

  bool isLifetimeEvent = false;

  vector<pair<double, int>> DOP_values = {}; 
  
  if(!annihilations.empty()) isLifetimeEvent = processHistograms(annihilations, DOP_values, stats, saveHistos, maxZ, maxThetaDiff, maxTimeDiff, maxDOP, 
                           sourcePos, testType, scatterTestValue, scatterTimeMin, scatterTimeMax, scatterAngleMin, scatterAngleMax);
  else if(!annihilationsMC.empty()) isLifetimeEvent = processHistograms(annihilationsMC, DOP_values, stats, saveHistos, maxZ, maxThetaDiff, maxTimeDiff, maxDOP, 
                           sourcePos, testType, scatterTestValue, scatterTimeMin, scatterTimeMax, scatterAngleMin, scatterAngleMax);
  else return false;

  for (uint i = 0; i < event.getHits().size(); i++) {
    auto promptHit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(i));
    auto promptHitMC = dynamic_cast<const JPetMCRecoHit*>(event.getHits().at(i));

      if(promptHit && fabs(promptHit->getPosZ())<maxZ && checkToT(promptHit, totCutDeexMin, totCutDeexMax)) {
        prompts.push_back(promptHit);
      }
      else if(promptHitMC && fabs(promptHitMC->getPosZ())<maxZ && checkToT(promptHitMC, totCutDeexMin, totCutDeexMax)) {
        promptsMC.push_back(promptHitMC);
      }
  }

  if (prompts.size() != 1 && promptsMC.size() != 1) {
      return false;
  }

  if(isLifetimeEvent){
    fillAnnihilationHistograms(annihilations, prompts, DOP_values, stats, sourcePos, totCutAnniMin, totCutAnniMax);
    fillAnnihilationHistograms(annihilationsMC, promptsMC, DOP_values, stats, sourcePos, totCutAnniMin, totCutAnniMax);
  }


  // Step 3: Analyze annihilation pairs and check for lifetime events
  return isLifetimeEvent;
}

bool EventCategorizerTools::checkFor3GammaLifetime( const JPetTimeWindowMC* timeWindowMC, const JPetEvent& event, vector<int> bad_ID, double maxZ, double minRelAngleCut, double minRelPhiCut, double minDistCut, double maxTimeDiff, double maxDOP, JPetStatistics& stats, 
                            bool saveHistos, double totCutAnniMin, double totCutAnniMax, double totCutDeexMin, double totCutDeexMax, const TVector3& sourcePos, 
                            ScatterTestType testType, double scatterTestValue, double scatterTimeMin, double scatterTimeMax, double scatterAngleMin, 
                            double scatterAngleMax)
{

  // Step 1: Initial checks and histogram filling for hits
  fillHitHistograms(event, stats, saveHistos, "none_3g_tot", "none_3g_z", "none_3g_z_ID");
  if (event.getHits().size() < 4) return false;
  fillHitHistograms(event, stats, saveHistos, "hits_3g_tot", "hits_3g_z", "hits_3g_z_ID");

  stats.fillHistogram("efficiency_3g", 1);

  std::vector<const JPetPhysRecoHit*> prompts;
  std::vector<std::vector<const JPetPhysRecoHit*>> annihilations;

  std::vector<const JPetMCRecoHit*> promptsMC;
  std::vector<std::vector<const JPetMCRecoHit*>> annihilationsMC;

  identifyAnnihilationHits(event, totCutAnniMin, totCutAnniMax, annihilations, annihilationsMC);

  vector<pair<double, int>> DOP_values = {};

  bool isLifetimeEvent = false;
  
  if(!annihilations.empty()) isLifetimeEvent = processHistograms(timeWindowMC, annihilations, DOP_values, stats, saveHistos, maxZ, minRelAngleCut, minRelPhiCut, minDistCut, maxTimeDiff, maxDOP, 
                           sourcePos, testType, scatterTestValue, scatterTimeMin, scatterTimeMax, scatterAngleMin, scatterAngleMax);
  else if(!annihilationsMC.empty()) isLifetimeEvent = processHistograms(timeWindowMC, annihilationsMC, DOP_values, stats, saveHistos, maxZ, minRelAngleCut, minRelPhiCut, minDistCut, maxTimeDiff, maxDOP, 
                           sourcePos, testType, scatterTestValue, scatterTimeMin, scatterTimeMax, scatterAngleMin, scatterAngleMax);
  else return false;

  //Check if any of the hits in the event is prompt based on TOT selection
  for (uint i = 0; i < event.getHits().size(); i++) {
    auto promptHit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(i));
    auto promptHitMC = dynamic_cast<const JPetMCRecoHit*>(event.getHits().at(i));

      if(promptHit && fabs(promptHit->getPosZ())<maxZ && checkToT(promptHit, totCutDeexMin, totCutDeexMax)) {
        prompts.push_back(promptHit);
      }
      else if(promptHitMC && fabs(promptHitMC->getPosZ())<maxZ && checkToT(promptHitMC, totCutDeexMin, totCutDeexMax)) {
        promptsMC.push_back(promptHitMC);
      }
  }

  if (prompts.size() != 1 && promptsMC.size() != 1) {
      return false;
  }

  sort(DOP_values.begin(), DOP_values.end());
  

  if(isLifetimeEvent){
    fillAnnihilationHistograms(timeWindowMC, annihilations, prompts, DOP_values, stats, sourcePos, totCutAnniMin, totCutAnniMax);
    fillAnnihilationHistograms(timeWindowMC, annihilationsMC, promptsMC, DOP_values, stats, sourcePos, totCutAnniMin, totCutAnniMax);
  }
  //cout<<"Prompt: "<<prompt_idx<<endl;
  // Then looking for annihilation back to back pairs

  return isLifetimeEvent;
}

template <typename HitType>
bool EventCategorizerTools::processHistograms(const std::vector<std::pair<const HitType*, const HitType*>>& annihilations, std::vector<std::pair<double, int>>& DOP_values,
                                                JPetStatistics& stats, bool saveHistos, double maxZ, double maxThetaDiff, double maxTimeDiff, double maxDOP, const TVector3& sourcePos, ScatterTestType testType, double scatterTestValue, 
                                                double scatterTimeMin, double scatterTimeMax, double scatterAngleMin, double scatterAngleMax) {

  bool isLifetimeEvent = false;
  int annih_ind = 0;

  for (const auto& pair2g : annihilations) {
    // Calculate the annihilation point position
    TVector3 annhilationPoint = calculateAnnihilationPoint(pair2g.first, pair2g.second);

    // Caculate event times - annihilation
    double annihTime1 = pair2g.first->getTime() - calculateTOF(pair2g.first, sourcePos);
    double annihTime2 = pair2g.second->getTime() - calculateTOF(pair2g.second, sourcePos);

    int scinID_annih1 = pair2g.first->getScin().getID();
    int scinID_annih2 = pair2g.second->getScin().getID();   

    auto tot1 = calculateToT(pair2g.first);
    auto tot2 = calculateToT(pair2g.second);

    double dist_annih1_annih2 = calculateDistance(pair2g.second, pair2g.first);

    double timeDiff = pair2g.second->getTime()-pair2g.first->getTime();
    double timeDiff_annih = annihTime2-annihTime1;

    TVector3 firstVec = pair2g.first->getPos();
    TVector3 secondVec = pair2g.second->getPos();
    double theta = TMath::RadToDeg() * firstVec.Angle(secondVec);

    double DOP = (annhilationPoint - sourcePos).Mag();

    bool thetaCut = checkRelativeAngles(pair2g.first->getPos(), pair2g.second->getPos(), maxThetaDiff);
    bool tDiffCut = fabs(timeDiff) < maxTimeDiff;
    bool scatterCut = !checkForScatter(pair2g.first, pair2g.second, stats, false, testType, scatterTestValue, scatterTimeMin, scatterTimeMax, scatterAngleMin, scatterAngleMax);
    bool DOPCut = DOP < maxDOP;
    bool zCut = fabs(pair2g.first->getPosZ()) < maxZ && fabs(pair2g.second->getPosZ()) < maxZ;

  
    if (saveHistos)
    {
      stats.fillHistogram("2g_tot", tot1);
      stats.fillHistogram("2g_tot", tot2);
      stats.fillHistogram("2g_theta", theta);
      stats.fillHistogram("2g_timeDiff_theta", timeDiff, theta);
      stats.fillHistogram("2g_timeDiff", timeDiff_annih);
      stats.fillHistogram("2g_DOP", DOP);
      stats.fillHistogram("2g_scatter_test_time", fabs(timeDiff) - dist_annih1_annih2/kLightVelocity_cm_ps);
      stats.fillHistogram("2g_dist", dist_annih1_annih2);
      stats.fillHistogram("2g_ID_all", scinID_annih1);
      stats.fillHistogram("2g_ID_all", scinID_annih2);
    // Scatter cut
    if (!(checkForScatter(pair2g.first, pair2g.second, stats, false, testType, scatterTestValue, scatterTimeMin, scatterTimeMax, scatterAngleMin,
                          scatterAngleMax))){
        scatterCut = true;
        stats.fillHistogram("scatter_2g_tot", tot1);
        stats.fillHistogram("scatter_2g_tot", tot2);
        stats.fillHistogram("scatter_2g_theta", theta);
        stats.fillHistogram("scatter_2g_timeDiff", timeDiff_annih);
        stats.fillHistogram("scatter_2g_scatter_test_time", timeDiff - dist_annih1_annih2/kLightVelocity_cm_ps);
    }
    if (zCut)
    {
        stats.fillHistogram("z_2g_tot", tot1);
        stats.fillHistogram("z_2g_tot", tot2);
        stats.fillHistogram("z_2g_theta", theta);
        stats.fillHistogram("z_2g_timeDiff", timeDiff_annih);
        stats.fillHistogram("z_2g_scatter_test_time", timeDiff - dist_annih1_annih2/kLightVelocity_cm_ps);
    }
    if (tDiffCut)
    {
        stats.fillHistogram("tdiff_2g_tot", tot1);
        stats.fillHistogram("tdiff_2g_tot", tot2);
        stats.fillHistogram("tdiff_2g_theta", theta);
        stats.fillHistogram("tdiff_2g_timeDiff", timeDiff_annih);
        stats.fillHistogram("tdiff_2g_scatter_test_time", fabs(timeDiff) - dist_annih1_annih2/kLightVelocity_cm_ps);
        stats.fillHistogram("tdiff_2g_dist", dist_annih1_annih2);
    }
      // Time difference cut
    if(thetaCut)
    {  
        stats.fillHistogram("theta_2g_tot", tot1);
        stats.fillHistogram("theta_2g_tot", tot2);
        stats.fillHistogram("theta_2g_theta", theta);
        stats.fillHistogram("theta_2g_timeDiff_theta", timeDiff, theta);
        stats.fillHistogram("theta_2g_timeDiff", timeDiff_annih);
        stats.fillHistogram("theta_2g_ID_all", scinID_annih1);
        stats.fillHistogram("theta_2g_ID_all", scinID_annih2);
        stats.fillHistogram("theta_2g_time_annih_12", pair2g.first->getTime(), pair2g.second->getTime());
        stats.fillHistogram("theta_2g_scatter_test_time", fabs(timeDiff) - dist_annih1_annih2/kLightVelocity_cm_ps);
        stats.fillHistogram("theta_2g_dist", dist_annih1_annih2);

    }  
    if(scatterCut && tDiffCut && thetaCut && DOPCut && zCut)
    {
      isLifetimeEvent = true;
      stats.fillHistogram("ap_2g_tot", tot1);
      stats.fillHistogram("ap_2g_tot", tot2);
      stats.fillHistogram("ap_2g_theta", theta);
      stats.fillHistogram("ap_2g_timeDiff", timeDiff_annih);
      stats.fillHistogram("ap_2g_scatter_test_time", timeDiff - dist_annih1_annih2/kLightVelocity_cm_ps);
      stats.fillHistogram("ap_xy", annhilationPoint.X(), annhilationPoint.Y());
      stats.fillHistogram("ap_zx", annhilationPoint.Z(), annhilationPoint.X());
      stats.fillHistogram("ap_zy", annhilationPoint.Z(), annhilationPoint.Y());
      stats.fillHistogram("ap_pos", annhilationPoint.Z(), annhilationPoint.X(), annhilationPoint.Y());
      stats.fillHistogram("ap_xy_zoom", annhilationPoint.X(), annhilationPoint.Y());
      stats.fillHistogram("ap_zx_zoom", annhilationPoint.Z(), annhilationPoint.X());
      stats.fillHistogram("ap_zy_zoom", annhilationPoint.Z(), annhilationPoint.Y());
      stats.fillHistogram("ap_pos_zoom", annhilationPoint.Z(), annhilationPoint.X(), annhilationPoint.Y());
      DOP_values.push_back(make_pair(DOP, annih_ind));
    }
    }
    annih_ind++;
  }

  return isLifetimeEvent;
}

template <typename HitType>
bool EventCategorizerTools::processHistograms(const JPetTimeWindowMC* timeWindowMC, const std::vector<std::vector<const HitType*>>& annihilations, std::vector<std::pair<double, int>>& DOP_values,
                                                JPetStatistics& stats, bool saveHistos, double maxZ, double minRelAngleCut, double minRelPhiCut, double minDistCut, double maxTimeDiff,  double maxDOP, const TVector3& sourcePos, ScatterTestType testType, double scatterTestValue, 
                                                double scatterTimeMin, double scatterTimeMax, double scatterAngleMin, double scatterAngleMax) {

  bool isLifetimeEvent = false;

  bool zCut_event = false, thetaCut_event = false, tDiffCut_event = false, DOPCut_event = false, vtxCut_event = false;
  bool zCut_signal = false, thetaCut_signal = false, tDiffCut_signal = false, DOPCut_signal = false, vtxCut_signal = false;

  int annih_ind = 0;

  for (auto pair3g : annihilations)
  {
    // Caculate event times - annihilation
    double annihTime1 = pair3g[0]->getTime()-calculateTOF(pair3g[0], sourcePos);
    double annihTime2 = pair3g[1]->getTime()-calculateTOF(pair3g[1], sourcePos);
    double annihTime3 = pair3g[2]->getTime()-calculateTOF(pair3g[2], sourcePos);

    auto tot1 = calculateToT(pair3g[0]);
    auto tot2 = calculateToT(pair3g[1]);
    auto tot3 = calculateToT(pair3g[2]);

    double DOP = calculatePlanePointDistance(pair3g[0], pair3g[1], pair3g[2], sourcePos);

    TVector3 ap = calculateAnnihilationPoint(*pair3g[0], *pair3g[1], *pair3g[2]);

    bool zCut = false, thetaCut = false, tDiffCut = false, DOPCut = false, distCut = false, vtxCut = false, phiCut = false;

    vector<double> relativeAngles;
    relativeAngles.push_back(TMath::RadToDeg() * pair3g[0]->getPos().Angle(pair3g[1]->getPos()));
    relativeAngles.push_back(TMath::RadToDeg() * pair3g[1]->getPos().Angle(pair3g[2]->getPos()));
    relativeAngles.push_back(TMath::RadToDeg() * pair3g[2]->getPos().Angle(pair3g[0]->getPos()));
    sort(relativeAngles.begin(), relativeAngles.end());

    vector<double> relative2DAngles;
    relative2DAngles.push_back(TMath::RadToDeg() * pair3g[0]->getPos().DeltaPhi(pair3g[1]->getPos()));
    relative2DAngles.push_back(TMath::RadToDeg() * pair3g[1]->getPos().DeltaPhi(pair3g[2]->getPos()));
    relative2DAngles.push_back(TMath::RadToDeg() * pair3g[2]->getPos().DeltaPhi(pair3g[0]->getPos()));

    double transformedX = relativeAngles.at(1) + relativeAngles.at(0);
    double transformedY = relativeAngles.at(1) - relativeAngles.at(0);

    vector<double> timeDiffs = {pair3g[1]->getTime() - pair3g[0]->getTime(), pair3g[2]->getTime() - pair3g[0]->getTime(), pair3g[2]->getTime() - pair3g[1]->getTime()};
    double timeDiff = *max_element(timeDiffs.begin(), timeDiffs.end(), 
                                  [](double a, double b){ return std::fabs(a) < std::fabs(b);});

    vector<double> timeDiffs_annih = {annihTime1 - annihTime2, annihTime2 - annihTime3, annihTime3 - annihTime1};
    double timeDiff_annih = *max_element(timeDiffs_annih.begin(), timeDiffs_annih.end(), 
                                  [](double a, double b){ return std::fabs(a) < std::fabs(b);});

    vector<double> dist_annih = {calculateDistance(pair3g[0], pair3g[1]), calculateDistance(pair3g[2], pair3g[0]), calculateDistance(pair3g[1], pair3g[2])};
    vector<double> dist_annih2D = {calculateDistance2D(pair3g[0], pair3g[1]), calculateDistance2D(pair3g[2], pair3g[1]), calculateDistance2D(pair3g[2], pair3g[0])};

    std::vector<int> gammaTagsTriplet = {-1000, -1000, -1000};
    int promptGammaTag = -1000;
    std::vector<int> gammaVtxTriplet = {-1000, -2000, -3000};
    int promptGammaVtx = -1000;


    for (size_t i = 0; i < pair3g.size(); ++i) {
      if (const auto* mcHitReco = dynamic_cast<const JPetMCRecoHit*>(pair3g[i])) {
        const auto mcHit = timeWindowMC->getMCHit<JPetRawMCHit>(mcHitReco->getMCindex());
        gammaTagsTriplet[i] = mcHit.getGammaTag();
        gammaVtxTriplet[i] = mcHit.getMCVtxIndex();
      }
    }
  

    bool allVerticesSame = true;
    if (!gammaVtxTriplet.empty()) {
      const auto refVtx = gammaVtxTriplet[0];
      for (const auto& vtx : gammaVtxTriplet) {
        if (vtx != refVtx) {
          allVerticesSame = false;
          break;
        }
      }
    }

    bool allTripletTagsAre3 = std::all_of(gammaTagsTriplet.begin(), gammaTagsTriplet.end(),
    [](int tag) { return tag == 3; });

    if (saveHistos)
    {
      stats.fillHistogram("3g_tot", tot1);
      stats.fillHistogram("3g_tot", tot2);
      stats.fillHistogram("3g_tot", tot3);
      stats.fillHistogram("3g_DOP", DOP);
      stats.fillHistogram("3g_rel_angles", transformedX, transformedY);
      for(int tdiff_ind = 0; tdiff_ind < timeDiffs.size(); tdiff_ind++) {
        stats.fillHistogram("3g_timeDiff", timeDiffs_annih[tdiff_ind]);
        stats.fillHistogram("3g_dist", dist_annih[tdiff_ind]);
        stats.fillHistogram("3g_phi", relative2DAngles[tdiff_ind]);
        stats.fillHistogram("3g_dist2D", dist_annih2D[tdiff_ind]);
        stats.fillHistogram("3g_scatter_test_time", timeDiffs[tdiff_ind] - dist_annih[tdiff_ind]/kLightVelocity_cm_ps);
      }
      stats.fillHistogram("3g_annihilation_point_xy", ap.X(), ap.Y());
      stats.fillHistogram("3g_annihilation_point_xz", ap.X(), ap.Z());

      if(fabs(pair3g[0]->getPosZ()) < maxZ && fabs(pair3g[1]->getPosZ()) < maxZ && fabs(pair3g[2]->getPosZ()) < maxZ ){
        zCut = true;
        zCut_event = true;
        if (allTripletTagsAre3 && allVerticesSame) zCut_signal = true;
        stats.fillHistogram("z_3g_DOP", DOP);
        stats.fillHistogram("z_3g_tot", tot1);
        stats.fillHistogram("z_3g_tot", tot2);
        stats.fillHistogram("z_3g_tot", tot3);
        stats.fillHistogram("z_3g_rel_angles", transformedX, transformedY);
        for(int tdiff_ind = 0; tdiff_ind < timeDiffs.size(); tdiff_ind++) {
          stats.fillHistogram("z_3g_timeDiff", timeDiffs_annih[tdiff_ind]);
          stats.fillHistogram("z_3g_phi", relative2DAngles[tdiff_ind]);
          stats.fillHistogram("z_3g_dist", dist_annih[tdiff_ind]);
          stats.fillHistogram("z_3g_dist2D", dist_annih2D[tdiff_ind]);
          stats.fillHistogram("z_3g_scatter_test_time", timeDiffs[tdiff_ind] - dist_annih[tdiff_ind]/kLightVelocity_cm_ps);
        }
        stats.fillHistogram("z_3g_annihilation_point_xy", ap.X(), ap.Y());
        stats.fillHistogram("z_3g_annihilation_point_xz", ap.X(), ap.Z());
      }

      if(DOP < maxDOP){
        DOPCut = true;
        DOPCut_event = true;
        if (allTripletTagsAre3 && allVerticesSame) DOPCut_signal = true;
        stats.fillHistogram("DOP_3g_DOP", DOP);
        stats.fillHistogram("DOP_3g_tot", tot1);
        stats.fillHistogram("DOP_3g_tot", tot2);
        stats.fillHistogram("DOP_3g_tot", tot3);
        stats.fillHistogram("DOP_3g_rel_angles", transformedX, transformedY);
        for(int tdiff_ind = 0; tdiff_ind < timeDiffs.size(); tdiff_ind++) {
          stats.fillHistogram("DOP_3g_timeDiff", timeDiffs_annih[tdiff_ind]);
          stats.fillHistogram("DOP_3g_phi", relative2DAngles[tdiff_ind]);
          stats.fillHistogram("DOP_3g_dist", dist_annih[tdiff_ind]);
          stats.fillHistogram("DOP_3g_dist2D", dist_annih2D[tdiff_ind]);
          stats.fillHistogram("DOP_3g_scatter_test_time", timeDiffs[tdiff_ind] - dist_annih[tdiff_ind]/kLightVelocity_cm_ps);
        }
        stats.fillHistogram("DOP_3g_annihilation_point_xy", ap.X(), ap.Y());
        stats.fillHistogram("DOP_3g_annihilation_point_xz", ap.X(), ap.Z());
      }
      if(fabs(timeDiff_annih) < maxTimeDiff){
        tDiffCut = true;
        tDiffCut_event = true;
        if (allTripletTagsAre3 && allVerticesSame) tDiffCut_signal = true;
        stats.fillHistogram("tdiff_3g_DOP", DOP);
        stats.fillHistogram("tdiff_3g_tot", tot1);
        stats.fillHistogram("tdiff_3g_tot", tot2);
        stats.fillHistogram("tdiff_3g_tot", tot3);
        stats.fillHistogram("tdiff_3g_rel_angles", transformedX, transformedY);
        for(int tdiff_ind = 0; tdiff_ind < timeDiffs.size(); tdiff_ind++) {
          stats.fillHistogram("tdiff_3g_timeDiff", timeDiffs_annih[tdiff_ind]);
          stats.fillHistogram("tdiff_3g_phi", relative2DAngles[tdiff_ind]);
          stats.fillHistogram("tdiff_3g_dist", dist_annih[tdiff_ind]);
          stats.fillHistogram("tdiff_3g_dist2D", dist_annih2D[tdiff_ind]);
          stats.fillHistogram("tdiff_3g_scatter_test_time", timeDiffs[tdiff_ind] - dist_annih[tdiff_ind]/kLightVelocity_cm_ps);
        }
        stats.fillHistogram("tdiff_3g_annihilation_point_xy", ap.X(), ap.Y());
        stats.fillHistogram("tdiff_3g_annihilation_point_xz", ap.X(), ap.Z());
      }
      if(transformedX > minRelAngleCut)
      {
        thetaCut=true;
        thetaCut_event=true;
        if (allTripletTagsAre3 && allVerticesSame) thetaCut_signal = true;
        stats.fillHistogram("theta_3g_DOP", DOP);
        stats.fillHistogram("theta_3g_tot", tot1);
        stats.fillHistogram("theta_3g_tot", tot2);
        stats.fillHistogram("theta_3g_tot", tot3);
        stats.fillHistogram("theta_3g_rel_angles", transformedX, transformedY);
        for(int tdiff_ind = 0; tdiff_ind < timeDiffs.size(); tdiff_ind++) {
          stats.fillHistogram("theta_3g_timeDiff", timeDiffs_annih[tdiff_ind]);
          stats.fillHistogram("theta_3g_phi", relative2DAngles[tdiff_ind]);
          stats.fillHistogram("theta_3g_dist", dist_annih[tdiff_ind]);
          stats.fillHistogram("theta_3g_dist2D", dist_annih2D[tdiff_ind]);
          stats.fillHistogram("theta_3g_scatter_test_time", timeDiffs[tdiff_ind] - dist_annih[tdiff_ind]/kLightVelocity_cm_ps);
        }          
        stats.fillHistogram("theta_3g_annihilation_point_xy", ap.X(), ap.Y());
        stats.fillHistogram("theta_3g_annihilation_point_xz", ap.X(), ap.Z());
      }
      if(fabs(relative2DAngles[0]) > minRelPhiCut && fabs(relative2DAngles[1]) > minRelPhiCut && fabs(relative2DAngles[2]) > minRelPhiCut)
      {
        phiCut=true;
        stats.fillHistogram("phi_3g_DOP", DOP);
        stats.fillHistogram("phi_3g_tot", tot1);
        stats.fillHistogram("phi_3g_tot", tot2);
        stats.fillHistogram("phi_3g_tot", tot3);
        stats.fillHistogram("phi_3g_rel_angles", transformedX, transformedY);
        for(int tdiff_ind = 0; tdiff_ind < timeDiffs.size(); tdiff_ind++) {
          stats.fillHistogram("phi_3g_timeDiff", timeDiffs_annih[tdiff_ind]);
          stats.fillHistogram("phi_3g_phi", relative2DAngles[tdiff_ind]);
          stats.fillHistogram("phi_3g_dist", dist_annih[tdiff_ind]);
          stats.fillHistogram("phi_3g_dist2D", dist_annih2D[tdiff_ind]);
          stats.fillHistogram("phi_3g_scatter_test_time", timeDiffs[tdiff_ind] - dist_annih[tdiff_ind]/kLightVelocity_cm_ps);
        }          
        stats.fillHistogram("phi_3g_annihilation_point_xy", ap.X(), ap.Y());
        stats.fillHistogram("phi_3g_annihilation_point_xz", ap.X(), ap.Z());
      }
      if(fabs(relative2DAngles[0]) <= minRelPhiCut || fabs(relative2DAngles[1]) <= minRelPhiCut || fabs(relative2DAngles[2]) <= minRelPhiCut)
      {
        stats.fillHistogram("phi0_3g_DOP", DOP);
        stats.fillHistogram("phi0_3g_tot", tot1);
        stats.fillHistogram("phi0_3g_tot", tot2);
        stats.fillHistogram("phi0_3g_tot", tot3);
        stats.fillHistogram("phi0_3g_rel_angles", transformedX, transformedY);
        for(int tdiff_ind = 0; tdiff_ind < timeDiffs.size(); tdiff_ind++) {
          stats.fillHistogram("phi0_3g_timeDiff", timeDiffs_annih[tdiff_ind]);
          stats.fillHistogram("phi0_3g_phi", relative2DAngles[tdiff_ind]);
          stats.fillHistogram("phi0_3g_dist", dist_annih[tdiff_ind]);
          stats.fillHistogram("phi0_3g_dist2D", dist_annih2D[tdiff_ind]);
          stats.fillHistogram("phi0_3g_scatter_test_time", timeDiffs[tdiff_ind] - dist_annih[tdiff_ind]/kLightVelocity_cm_ps);
        }          
        stats.fillHistogram("phi0_3g_annihilation_point_xy", ap.X(), ap.Y());
        stats.fillHistogram("phi0_3g_annihilation_point_xz", ap.X(), ap.Z());
      }


      if(*min_element(dist_annih.begin(), dist_annih.end()) > minDistCut)
      {
        distCut=true;
        stats.fillHistogram("dist_3g_DOP", DOP);
        stats.fillHistogram("dist_3g_tot", tot1);
        stats.fillHistogram("dist_3g_tot", tot2);
        stats.fillHistogram("dist_3g_tot", tot3);
        stats.fillHistogram("dist_3g_rel_angles", transformedX, transformedY);
        for(int tdiff_ind = 0; tdiff_ind < timeDiffs.size(); tdiff_ind++) {
          stats.fillHistogram("dist_3g_timeDiff", timeDiffs_annih[tdiff_ind]);
          stats.fillHistogram("dist_3g_phi", relative2DAngles[tdiff_ind]);
          stats.fillHistogram("dist_3g_dist", dist_annih[tdiff_ind]);
          stats.fillHistogram("dist_3g_dist2D", dist_annih2D[tdiff_ind]);
          stats.fillHistogram("dist_3g_scatter_test_time", timeDiffs[tdiff_ind] - dist_annih[tdiff_ind]/kLightVelocity_cm_ps);
        }          
        stats.fillHistogram("dist_3g_annihilation_point_xy", ap.X(), ap.Y());
        stats.fillHistogram("dist_3g_annihilation_point_xz", ap.X(), ap.Z());
      }
      
      if(ap.X()*ap.X()+ap.Y()*ap.Y() < 10000000*20.*20. && fabs(ap.Z())<10000000*10.)
      {
        vtxCut=true;
        vtxCut_event=true;
        if (allTripletTagsAre3 && allVerticesSame) vtxCut_signal = true;
        stats.fillHistogram("vtx_3g_DOP", DOP);
        stats.fillHistogram("vtx_3g_tot", tot1);
        stats.fillHistogram("vtx_3g_tot", tot2);
        stats.fillHistogram("vtx_3g_tot", tot3);
        stats.fillHistogram("vtx_3g_rel_angles", transformedX, transformedY);
        for(int tdiff_ind = 0; tdiff_ind < timeDiffs.size(); tdiff_ind++) {
          stats.fillHistogram("vtx_3g_timeDiff", timeDiffs_annih[tdiff_ind]);
          stats.fillHistogram("vtx_3g_phi", relative2DAngles[tdiff_ind]);
          stats.fillHistogram("vtx_3g_dist", dist_annih[tdiff_ind]);
          stats.fillHistogram("vtx_3g_dist2D", dist_annih2D[tdiff_ind]);
          stats.fillHistogram("vtx_3g_scatter_test_time", timeDiffs[tdiff_ind] - dist_annih[tdiff_ind]/kLightVelocity_cm_ps);
        }          
        stats.fillHistogram("vtx_3g_annihilation_point_xy", ap.X(), ap.Y());
        stats.fillHistogram("vtx_3g_annihilation_point_xz", ap.X(), ap.Z());
      }
      if(zCut&&thetaCut&&tDiffCut&&DOPCut&&vtxCut){
        DOP_values.push_back(make_pair(DOP, annih_ind));
        isLifetimeEvent = true;
      }

      annih_ind++;  

    }
  }

  if(zCut_event) stats.fillHistogram("efficiency_3g", 2);
  if(thetaCut_event) stats.fillHistogram("efficiency_3g", 3);
  if(tDiffCut_event) stats.fillHistogram("efficiency_3g", 4);
  if(DOPCut_event) stats.fillHistogram("efficiency_3g", 5);
  if(vtxCut_event) stats.fillHistogram("efficiency_3g", 6);

  if(zCut_signal) stats.fillHistogram("efficiency_3g_signal", 2);
  if(thetaCut_signal) stats.fillHistogram("efficiency_3g_signal", 3);
  if(tDiffCut_signal) stats.fillHistogram("efficiency_3g_signal", 4);
  if(DOPCut_signal) stats.fillHistogram("efficiency_3g_signal", 5);
  if(vtxCut_signal) stats.fillHistogram("efficiency_3g_signal", 6);

  sort(DOP_values.begin(), DOP_values.end());

  if(isLifetimeEvent){
    auto pair3g = annihilations[DOP_values[0].second];

    double annihTime1 = pair3g[0]->getTime()-calculateTOF(pair3g[0], sourcePos);
    double annihTime2 = pair3g[1]->getTime()-calculateTOF(pair3g[1], sourcePos);
    double annihTime3 = pair3g[2]->getTime()-calculateTOF(pair3g[2], sourcePos);

    auto tot1 = calculateToT(pair3g[0]);
    auto tot2 = calculateToT(pair3g[1]);
    auto tot3 = calculateToT(pair3g[2]);

    TVector3 ap = calculateAnnihilationPoint(*pair3g[0], *pair3g[1], *pair3g[2]);

    vector<double> relativeAngles;
    relativeAngles.push_back(TMath::RadToDeg() * pair3g[0]->getPos().Angle(pair3g[1]->getPos()));
    relativeAngles.push_back(TMath::RadToDeg() * pair3g[1]->getPos().Angle(pair3g[2]->getPos()));
    relativeAngles.push_back(TMath::RadToDeg() * pair3g[2]->getPos().Angle(pair3g[0]->getPos()));
    sort(relativeAngles.begin(), relativeAngles.end());
    
    double transformedX = relativeAngles.at(1) + relativeAngles.at(0);
    double transformedY = relativeAngles.at(1) - relativeAngles.at(0);

    vector<double> relative2DAngles;
    relative2DAngles.push_back(TMath::RadToDeg() * pair3g[0]->getPos().DeltaPhi(pair3g[1]->getPos()));
    relative2DAngles.push_back(TMath::RadToDeg() * pair3g[1]->getPos().DeltaPhi(pair3g[2]->getPos()));
    relative2DAngles.push_back(TMath::RadToDeg() * pair3g[2]->getPos().DeltaPhi(pair3g[0]->getPos()));

    vector<double> timeDiffs = {pair3g[1]->getTime() - pair3g[0]->getTime(), pair3g[2]->getTime() - pair3g[1]->getTime(), pair3g[2]->getTime() - pair3g[0]->getTime()};
    double timeDiff = *max_element(timeDiffs.begin(), timeDiffs.end(), 
                                    [](double a, double b){ return std::fabs(a) < std::fabs(b);});

    vector<double> timeDiffs_annih = {annihTime1 - annihTime2, annihTime2 - annihTime3, annihTime3 - annihTime1};
    double timeDiff_annih = *max_element(timeDiffs_annih.begin(), timeDiffs_annih.end(), 
                                  [](double a, double b){ return std::fabs(a) < std::fabs(b);});
                                    
    vector<double> dist_annih = {calculateDistance(pair3g[0], pair3g[1]), calculateDistance(pair3g[2], pair3g[1]), calculateDistance(pair3g[2], pair3g[0])};
    vector<double> dist_annih2D = {calculateDistance2D(pair3g[0], pair3g[1]), calculateDistance2D(pair3g[2], pair3g[1]), calculateDistance2D(pair3g[2], pair3g[0])};

    std::vector<int> gammaTagsTriplet = {-1000, -1000, -1000};
    int promptGammaTag = -1000;
    std::vector<int> gammaVtxTriplet = {-1000, -2000, -3000};
    int promptGammaVtx = -1000;

    for (size_t i = 0; i < pair3g.size(); ++i) {
      if (const auto* mcHitReco = dynamic_cast<const JPetMCRecoHit*>(pair3g[i])) {
        const auto mcHit = timeWindowMC->getMCHit<JPetRawMCHit>(mcHitReco->getMCindex());
        gammaTagsTriplet[i] = mcHit.getGammaTag();
        gammaVtxTriplet[i] = mcHit.getMCVtxIndex();
      }
    }
  

    bool allVerticesSame = true;
    if (!gammaVtxTriplet.empty()) {
      const auto refVtx = gammaVtxTriplet[0];
      for (const auto& vtx : gammaVtxTriplet) {
        if (vtx != refVtx) {
          allVerticesSame = false;
          break;
        }
      }
    }

    bool allTripletTagsAre3 = std::all_of(gammaTagsTriplet.begin(), gammaTagsTriplet.end(),
    [](int tag) { return tag == 3; });

    stats.fillHistogram("efficiency_3g", 7);
    if(allVerticesSame && allTripletTagsAre3) stats.fillHistogram("efficiency_3g_signal", 7);

    stats.fillHistogram("ap_3g_tot", tot1);
    stats.fillHistogram("ap_3g_tot", tot2);
    stats.fillHistogram("ap_3g_tot", tot3);
    stats.fillHistogram("ap_3g_DOP", DOP_values[0].first);
    stats.fillHistogram("ap_3g_rel_angles", transformedX, transformedY);
    for(int tdiff_ind = 0; tdiff_ind < timeDiffs.size(); tdiff_ind++) {
      stats.fillHistogram("ap_3g_timeDiff", timeDiffs[tdiff_ind]);
      stats.fillHistogram("ap_3g_timeDiff_annih", timeDiffs_annih[tdiff_ind]);
      stats.fillHistogram("ap_3g_phi", relative2DAngles[tdiff_ind]);
      stats.fillHistogram("ap_3g_dist", dist_annih[tdiff_ind]);
      stats.fillHistogram("ap_3g_dist2D", dist_annih2D[tdiff_ind]);
      stats.fillHistogram("ap_3g_scatter_test_time", timeDiffs[tdiff_ind] - dist_annih[tdiff_ind]/kLightVelocity_cm_ps);
      stats.fillHistogram("ap_3g_z", pair3g[tdiff_ind]->getPosZ());
    }            
    stats.fillHistogram("ap_3g_annihilation_point_xy", ap.X(), ap.Y()); 
    stats.fillHistogram("ap_3g_annihilation_point_xz", ap.X(), ap.Z()); 
  }

  return isLifetimeEvent;
}

template <typename HitType>
void EventCategorizerTools::fillAnnihilationHistograms(const std::vector<std::pair<const HitType*, const HitType*>>& annihilations,
                                                       const std::vector<const HitType*>& prompts, std::vector<std::pair<double, int>>& DOP_values, JPetStatistics& stats, 
                                                       const TVector3& sourcePos, double totCutAnniMin, double totCutAnniMax) {
  if (!annihilations.empty() && prompts.size() == 1) {
    const auto& pair2g = annihilations[DOP_values[0].second];
    double annihTime1 = pair2g.first->getTime() - calculateTOF(pair2g.first, sourcePos);
    double annihTime2 = pair2g.second->getTime() - calculateTOF(pair2g.second, sourcePos);
    double promptTime = prompts[0]->getTime() - calculateTOF(prompts[0], sourcePos);

    TVector3 annhilationPoint = calculateAnnihilationPoint(pair2g.first, pair2g.second);\

    double lifetime = (annihTime1 + annihTime2) / 2.0 - promptTime;
    double dist_annih1_annih2 = calculateDistance(pair2g.second, pair2g.first);

      // Fill histograms
    stats.fillHistogram("lifetime_ap_2g_prompt", lifetime);
    stats.fillHistogram("lifetime_ap_2g_prompt_zoom", lifetime);
    stats.fillHistogram("ap_2g_tot_lifetime", calculateToT(pair2g.first));
    stats.fillHistogram("ap_2g_tot_lifetime", calculateToT(pair2g.second));
    stats.fillHistogram("ap_2g_timeDiff_lifetime", pair2g.first->getTime() - pair2g.second->getTime());
    stats.fillHistogram("ap_xy_lifetime", annhilationPoint.X(), annhilationPoint.Y());
    stats.fillHistogram("ap_zx_lifetime", annhilationPoint.Z(), annhilationPoint.X());
    stats.fillHistogram("ap_zy_lifetime", annhilationPoint.Z(), annhilationPoint.Y());
    stats.fillHistogram("ap_pos_lifetime", annhilationPoint.Z(), annhilationPoint.X(), annhilationPoint.Y());
    
  }
}


template <typename HitType>
void EventCategorizerTools::fillAnnihilationHistograms(const JPetTimeWindowMC* timeWindowMC, const std::vector<std::vector<const HitType*>>& annihilations,
                                                       const std::vector<const HitType*>& prompts, std::vector<pair<double, int>>& DOP_values, JPetStatistics& stats, 
                                                       const TVector3& sourcePos, double totCutAnniMin, double totCutAnniMax) {
  if (!annihilations.empty() && prompts.size() == 1) {
    stats.fillHistogram("3g_stats_multi_annihilations", annihilations.size());
    stats.fillHistogram("3g_stats_multi_prompts", prompts.size());

 {     //stats.fillHistogram("3g_cut_stats", 6);  

      auto pair3g = annihilations[DOP_values[0].second];

      std::vector<int> gammaTagsTriplet = {-1000, -1000, -1000};
      int promptGammaTag = -1000;
      std::vector<int> gammaVtxTriplet = {-1000, -2000, -3000};
      int promptGammaVtx = -1000;

      if(timeWindowMC){
        for (size_t i = 0; i < pair3g.size(); ++i) {
          if (const auto* mcHitReco = dynamic_cast<const JPetMCRecoHit*>(pair3g[i])) {
            const auto mcHit = timeWindowMC->getMCHit<JPetRawMCHit>(mcHitReco->getMCindex());
            gammaTagsTriplet[i] = mcHit.getGammaTag();
            gammaVtxTriplet[i] = mcHit.getMCVtxIndex();
          }
        }
        
        if (!prompts.empty()) {
          if (const auto* promptMcHitReco = dynamic_cast<const JPetMCRecoHit*>(prompts[0])) {
            const auto promptMcHit = timeWindowMC->getMCHit<JPetRawMCHit>(promptMcHitReco->getMCindex());
            promptGammaTag = promptMcHit.getGammaTag();
            promptGammaVtx = promptMcHit.getMCVtxIndex();
          }
        }
      }

      bool allVerticesSame = true;
      if (!gammaVtxTriplet.empty()) {
        const auto refVtx = gammaVtxTriplet[0];
        for (const auto& vtx : gammaVtxTriplet) {
          if (vtx != refVtx) {
            allVerticesSame = false;
            break;
          }
        }
        if (allVerticesSame && !prompts.empty()) {
          if (refVtx != promptGammaVtx) {
            allVerticesSame = false;
          }
        }
      }

      double annihTime1 = pair3g[0]->getTime()-calculateTOF(pair3g[0], sourcePos);
      double annihTime2 = pair3g[1]->getTime()-calculateTOF(pair3g[1], sourcePos);
      double annihTime3 = pair3g[2]->getTime()-calculateTOF(pair3g[2], sourcePos);
      double promptTime = prompts[0]->getTime()-calculateTOF(prompts[0], sourcePos);

      auto tot1 = calculateToT(pair3g[0]);
      auto tot2 = calculateToT(pair3g[1]);
      auto tot3 = calculateToT(pair3g[2]);

      TVector3 ap = calculateAnnihilationPoint(*pair3g[0], *pair3g[1], *pair3g[2]);

      // Calculate lifetime
      double lifetime = (annihTime1 + annihTime2 + annihTime3) / 3.0 - promptTime;

      vector<double> relativeAngles;
      relativeAngles.push_back(TMath::RadToDeg() * pair3g[0]->getPos().Angle(pair3g[1]->getPos()));
      relativeAngles.push_back(TMath::RadToDeg() * pair3g[1]->getPos().Angle(pair3g[2]->getPos()));
      relativeAngles.push_back(TMath::RadToDeg() * pair3g[2]->getPos().Angle(pair3g[0]->getPos()));
      sort(relativeAngles.begin(), relativeAngles.end());
      
      double transformedX = relativeAngles.at(1) + relativeAngles.at(0);
      double transformedY = relativeAngles.at(1) - relativeAngles.at(0);

      vector<double> relative2DAngles;
      relative2DAngles.push_back(TMath::RadToDeg() * pair3g[0]->getPos().DeltaPhi(pair3g[1]->getPos()));
      relative2DAngles.push_back(TMath::RadToDeg() * pair3g[1]->getPos().DeltaPhi(pair3g[2]->getPos()));
      relative2DAngles.push_back(TMath::RadToDeg() * pair3g[2]->getPos().DeltaPhi(pair3g[0]->getPos()));

      vector<double> timeDiffs = {pair3g[1]->getTime() - pair3g[0]->getTime(), pair3g[2]->getTime() - pair3g[1]->getTime(), pair3g[2]->getTime() - pair3g[0]->getTime()};
      double timeDiff = *max_element(timeDiffs.begin(), timeDiffs.end(), 
                                      [](double a, double b){ return std::fabs(a) < std::fabs(b);});

    vector<double> timeDiffs_annih = {annihTime1 - annihTime2, annihTime2 - annihTime3, annihTime3 - annihTime1};
    double timeDiff_annih = *max_element(timeDiffs_annih.begin(), timeDiffs_annih.end(), 
                                  [](double a, double b){ return std::fabs(a) < std::fabs(b);});

      vector<double> dist_annih = {calculateDistance(pair3g[0], pair3g[1]), calculateDistance(pair3g[2], pair3g[1]), calculateDistance(pair3g[2], pair3g[0])};
      vector<double> dist_annih2D = {calculateDistance2D(pair3g[0], pair3g[1]), calculateDistance2D(pair3g[2], pair3g[1]), calculateDistance2D(pair3g[2], pair3g[0])};

      bool allTripletTagsAre3 = std::all_of(gammaTagsTriplet.begin(), gammaTagsTriplet.end(),
                                      [](int tag) { return tag == 3; });

      if (allTripletTagsAre3 && promptGammaTag == 1 && allVerticesSame) {
        stats.fillHistogram("efficiency_3g_signal", 8);
        stats.fillHistogram("ap_3g_tot_signal", tot1);
        stats.fillHistogram("ap_3g_tot_signal", tot2);
        stats.fillHistogram("ap_3g_tot_signal", tot3);
        stats.fillHistogram("ap_3g_DOP_signal", DOP_values[0].first);
        stats.fillHistogram("ap_3g_rel_angles_signal", transformedX, transformedY);
        for(int tdiff_ind = 0; tdiff_ind < timeDiffs.size(); tdiff_ind++) {
          stats.fillHistogram("ap_3g_timeDiff_signal", timeDiffs_annih[tdiff_ind]);
        }
      }
      else{
        stats.fillHistogram("ap_3g_tot_bkg", tot1);
        stats.fillHistogram("ap_3g_tot_bkg", tot2);
        stats.fillHistogram("ap_3g_tot_bkg", tot3); 
        stats.fillHistogram("ap_3g_DOP_bkg", DOP_values[0].first);
        stats.fillHistogram("ap_3g_rel_angles_bkg", transformedX, transformedY);
        for(int tdiff_ind = 0; tdiff_ind < timeDiffs.size(); tdiff_ind++) {
          stats.fillHistogram("ap_3g_timeDiff_bkg", timeDiffs_annih[tdiff_ind]);
        }
      }
      
      stats.fillHistogram("efficiency_3g", 8);
      stats.fillHistogram("lifetime_ap_3g_prompt", lifetime);
      stats.fillHistogram("lifetime_ap_3g_prompt_zoom", lifetime);
      stats.fillHistogram("ap_3g_tot_lifetime", tot1);
      stats.fillHistogram("ap_3g_tot_lifetime", tot2);
      stats.fillHistogram("ap_3g_tot_lifetime", tot3);
      stats.fillHistogram("ap_3g_DOP_lifetime", DOP_values[0].first);
      stats.fillHistogram("ap_3g_rel_angles_lifetime", transformedX, transformedY);
      for(int tdiff_ind = 0; tdiff_ind < timeDiffs.size(); tdiff_ind++) {
        stats.fillHistogram("ap_3g_timeDiff_lifetime", timeDiffs_annih[tdiff_ind]);
        stats.fillHistogram("ap_3g_phi_lifetime", relative2DAngles[tdiff_ind]);
        stats.fillHistogram("ap_3g_dist_lifetime", dist_annih[tdiff_ind]);
        stats.fillHistogram("ap_3g_dist2D_lifetime", dist_annih2D[tdiff_ind]);
        stats.fillHistogram("ap_3g_scatter_test_time_lifetime", timeDiffs[tdiff_ind] - dist_annih[tdiff_ind]/kLightVelocity_cm_ps);
      }            
      stats.fillHistogram("ap_3g_annihilation_point_xy_lifetime", ap.X(), ap.Y()); 
      stats.fillHistogram("ap_3g_annihilation_point_xz_lifetime", ap.X(), ap.Z()); 
      stats.fillHistogram("ap_3g_stats_prompts", prompts.size());}
  }
}


void EventCategorizerTools::fillHitHistograms(const JPetEvent& event, JPetStatistics& stats, bool saveHistos, const std::string& histogramName, const std::string& histogramName_z, const std::string& histogramName_z_ID) {
  if(saveHistos) {
    for (uint i = 0; i < event.getHits().size(); i++) {
      auto hit = event.getHits().at(i);
      if (const auto* firstHit = dynamic_cast<const JPetPhysRecoHit*>(hit)) {
          stats.fillHistogram(histogramName.c_str(), calculateToT(firstHit));
          stats.fillHistogram(histogramName_z.c_str(), static_cast<double>(firstHit->getPosZ()));
          stats.fillHistogram(histogramName_z_ID.c_str(), firstHit->getScin().getID(), static_cast<double>(firstHit->getPosZ()));
      } else if (const auto* firstHitMC = dynamic_cast<const JPetMCRecoHit*>(hit)) {
          stats.fillHistogram(histogramName.c_str(), calculateToT(firstHitMC));
          stats.fillHistogram(histogramName_z.c_str(), static_cast<double>(firstHitMC->getPosZ()));
      }
    }
  }
}

void EventCategorizerTools::identifyAnnihilationHits(const JPetEvent& event, double totCutAnniMin, double totCutAnniMax,
  std::vector<std::pair<const JPetPhysRecoHit*, const JPetPhysRecoHit*>>& annihilations, std::vector<std::pair<const JPetMCRecoHit*, const JPetMCRecoHit*>>& annihilationsMC
) {
  std::vector<const JPetPhysRecoHit*> temp;
  std::vector<const JPetMCRecoHit*> tempMC;

  std::random_device rd;
  std::default_random_engine rng(rd());

  const auto& hits = event.getHits();
  const size_t nHits = hits.size();

  for (size_t i = 0; i < nHits; ++i) {
    auto firstHit = dynamic_cast<const JPetPhysRecoHit*>(hits[i]);
    auto firstHitMC = dynamic_cast<const JPetMCRecoHit*>(hits[i]);

    for (size_t j = i + 1; j < nHits; ++j) {
      auto secondHit = dynamic_cast<const JPetPhysRecoHit*>(hits[j]);
      auto secondHitMC = dynamic_cast<const JPetMCRecoHit*>(hits[j]);

      temp.clear();
      tempMC.clear();

      if (firstHit && secondHit &&
          checkToT(firstHit, totCutAnniMin, totCutAnniMax) &&
          checkToT(secondHit, totCutAnniMin, totCutAnniMax)) {
        
        temp = {firstHit, secondHit};
        std::shuffle(temp.begin(), temp.end(), rng);
        annihilations.emplace_back(temp[0], temp[1]);
      }
      else if (firstHitMC && secondHitMC &&
               checkToT(firstHitMC, totCutAnniMin, totCutAnniMax) &&
               checkToT(secondHitMC, totCutAnniMin, totCutAnniMax)) {
        
        tempMC = {firstHitMC, secondHitMC};
        std::shuffle(tempMC.begin(), tempMC.end(), rng);
        annihilationsMC.emplace_back(tempMC[0], tempMC[1]);
      }
    }
  }
}


void EventCategorizerTools::identifyAnnihilationHits(const JPetEvent& event, double totCutAnniMin, double totCutAnniMax, 
                               std::vector<std::vector<const JPetPhysRecoHit*>>& annihilations,
                               std::vector<std::vector<const JPetMCRecoHit*>>& annihilationsMC) {

    std::vector<const JPetPhysRecoHit*> temp = {};
    std::vector<const JPetMCRecoHit*> tempMC = {};
    for (uint i = 0; i < event.getHits().size(); i++) {
        auto firstHit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(i));
        auto firstHitMC = dynamic_cast<const JPetMCRecoHit*>(event.getHits().at(i));
        for (uint j = i + 1; j < event.getHits().size(); j++) {
            auto secondHit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(j));
            auto secondHitMC = dynamic_cast<const JPetMCRecoHit*>(event.getHits().at(j));
            for (uint k = j + 1; k < event.getHits().size(); k++) {
              auto thirdHit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(k));
              auto thirdHitMC = dynamic_cast<const JPetMCRecoHit*>(event.getHits().at(k));
              temp = {};
              tempMC = {};
              if(firstHit && secondHit && thirdHit && checkToT(firstHit, totCutAnniMin, totCutAnniMax) && checkToT(secondHit, totCutAnniMin, totCutAnniMax) && checkToT(thirdHit, totCutAnniMin, totCutAnniMax)) {
                  temp={firstHit, secondHit, thirdHit};
                  std::random_device rd;
                  std::default_random_engine rng(rd());
                  std::shuffle(temp.begin(), temp.end(), rng);
                  //std::sort(temp.begin(), temp.end(), [](const JPetPhysRecoHit* a, const JPetPhysRecoHit* b) {return calculateToT(a) > calculateToT(b);});
                  annihilations.emplace_back(temp);
              }
              else if(firstHitMC && secondHitMC && thirdHitMC && checkToT(firstHitMC, totCutAnniMin, totCutAnniMax) && checkToT(secondHitMC, totCutAnniMin, totCutAnniMax) && checkToT(thirdHitMC, totCutAnniMin, totCutAnniMax)) {
                  tempMC={firstHitMC, secondHitMC, thirdHitMC};
                  std::random_device rd;
                  std::default_random_engine rng(rd());
                  std::shuffle(tempMC.begin(), tempMC.end(), rng);
                  //std::sort(tempMC.begin(), tempMC.end(), [](const JPetMCRecoHit* a, const JPetMCRecoHit* b) {return calculateToT(a) > calculateToT(b);});
                  annihilationsMC.emplace_back(tempMC);
              }
            }
        }
    }
}


bool EventCategorizerTools::checkToT(const JPetPhysRecoHit* hit, double minToT, double maxToT)
{
  return (hit->getToT() > minToT && hit->getToT() < maxToT);
}

bool EventCategorizerTools::checkToT(const JPetMCRecoHit* hit, double minToT, double maxToT)
{
  return (hit->getEnergy() > minToT && hit->getEnergy() < maxToT);
}

double EventCategorizerTools::calculateToT(const JPetPhysRecoHit* hit)
{
  return ToTToEnergy(hit->getToT()/pow(10,6));
}

double EventCategorizerTools::calculateToT(const JPetMCRecoHit* hit)
{
  return EnergyToToT(hit->getEnergy())*pow(10,6);
}

double EventCategorizerTools::EnergyToToT(double hit_energy)
{
  double A0 = 5.6;
  double A1 = 7.0;
  double A2 = 0.9944;
  return A0 - A1 * pow(A2, hit_energy);
}

double EventCategorizerTools::ToTToEnergy(double hit_tot)
{
  double A0 = 5.6;
  double A1 = 7.0;
  double A2 = 0.9944;
  double x = (A0 - hit_tot)/A1;
  return log10(x)/log10(A2);
}


bool EventCategorizerTools::checkRelativeAngles(const TVector3& pos1, const TVector3& pos2, double maxThetaDiff)
{
  double theta = TMath::RadToDeg() * pos1.Angle(pos2);
  return (180.0 - theta < maxThetaDiff);
}

/**
 * Checking if pair of hits meet scattering condition
 */
bool EventCategorizerTools::checkForScatter(const JPetBaseHit* primaryHit, const JPetBaseHit* scatterHit, JPetStatistics& stats, bool saveHistos,
                                            ScatterTestType testType, double scatterTestValue, double scatterTimeMin, double scatterTimeMax,
                                            double scatterAngleMin, double scatterAngleMax)
{
  bool isScatter = false;

  double dist = calculateDistance(primaryHit, scatterHit);
  double timeDiff = scatterHit->getTime() - primaryHit->getTime();
  double testTimeRel = timeDiff - dist / kLightVelocity_cm_ps;
  double testTimeAbs = fabs(testTimeRel);
  double testDistRel = timeDiff * kLightVelocity_cm_ps - dist;
  double testDistAbs = fabs(testDistRel);
  double scatterAngle = calculateScatteringAngle(primaryHit, scatterHit);

  if (saveHistos)
  {
    stats.fillHistogram("scatter_test_time_rel", testTimeRel);
    stats.fillHistogram("scatter_test_time_abs", testTimeAbs);
    stats.fillHistogram("scatter_test_dist_rel", testDistRel);
    stats.fillHistogram("scatter_test_dist_abs", testDistAbs);

    stats.fillHistogram("scatter_angle_time", testTimeRel, scatterAngle);
    stats.fillHistogram("scatter_angle_time_small", testTimeRel, scatterAngle);
  }

  if (testType == EventCategorizerTools::kSimpleParam)
  {
    isScatter = testTimeAbs < scatterTestValue;
  }

  if (testType == EventCategorizerTools::kMinMaxParams)
  {
    isScatter = !(scatterTimeMin < testTimeRel && testTimeRel < scatterTimeMax && scatterAngleMin < scatterAngle && scatterAngle < scatterAngleMax);
  }

  if (saveHistos)
  {
    if (isScatter)
    {
      stats.fillHistogram("scatter_test_rel_pass", testTimeRel);
      stats.fillHistogram("scatter_test_abs_pass", testTimeAbs);
      stats.fillHistogram("scatter_angle_time_pass", testTimeRel, scatterAngle);
    }
    else
    {
      stats.fillHistogram("scatter_test_rel_fail", testTimeRel);
      stats.fillHistogram("scatter_test_abs_fail", testTimeAbs);
      stats.fillHistogram("scatter_angle_time_fail", testTimeRel, scatterAngle);
    }
  }

  return isScatter;
}

/**
 * Calculation of distance between two hits
 */
double EventCategorizerTools::calculateDistance(const JPetBaseHit* hit1, const JPetBaseHit* hit2) { return (hit1->getPos() - hit2->getPos()).Mag(); }

double EventCategorizerTools::calculateDistance2D(const JPetBaseHit* hit1, const JPetBaseHit* hit2) { 
  return sqrt((hit1->getPosX()-hit2->getPosX())*(hit1->getPosX()-hit2->getPosX()) + (hit1->getPosY()-hit2->getPosY())*(hit1->getPosY()-hit2->getPosY())); }


/**
 * Calculation of time that light needs to travel the distance between primary gamma
 * and scattered gamma. Return value in picoseconds.
 */
double EventCategorizerTools::calculateScatteringTime(const JPetBaseHit* hit1, const JPetBaseHit* hit2)
{
  return calculateDistance(hit1, hit2) / kLightVelocity_cm_ps;
}

/**
 * Calculation of scatter angle between primary hit and scattered hit.
 * This function assumes that source of first gamma was in (0,0,0).
 * Angle is calculated from scalar product, return value in degrees.
 */
double EventCategorizerTools::calculateScatteringAngle(const JPetBaseHit* hit1, const JPetBaseHit* hit2)
{
  return TMath::RadToDeg() * hit1->getPos().Angle(hit2->getPos() - hit1->getPos());
}

double EventCategorizerTools::calculateTOFByConvention(const JPetBaseHit* hitA, const JPetBaseHit* hitB)
{
  if (hitA->getScin().getSlot().getTheta() < hitB->getScin().getSlot().getTheta())
  {
    return calculateTOF(hitA, hitB);
  }
  else
  {
    return calculateTOF(hitB, hitA);
  }
}

double EventCategorizerTools::calculateTOF(const JPetBaseHit* hitA, const JPetBaseHit* hitB)
{
  return EventCategorizerTools::calculateTOF(hitA->getTime(), hitB->getTime());
}

double EventCategorizerTools::calculateTOF(double time1, double time2) { return (time1 - time2); }

double EventCategorizerTools::calculateTOF(const JPetBaseHit* hit, const TVector3& sourcePos)  { return (hit->getPos()-sourcePos).Mag()/kLightVelocity_cm_ps;} 

/**
 * @brief Calculation of an annihilation point based on LOR and TOFof two hits.
 *
 * Line of Response between two hits is used to estimate annihilation point by calculating its
 * middle and shifting it along LOR based on Time of Flight value toward one hit or the other.
 */
TVector3 EventCategorizerTools::calculateAnnihilationPoint(const JPetBaseHit* hit1, const JPetBaseHit* hit2)
{
  TVector3 middleOfLOR = 0.5 * (hit1->getPos() + hit2->getPos());
  TVector3 versorOnLOR = (hit2->getPos() - hit1->getPos()).Unit();

  double tof = EventCategorizerTools::calculateTOF(hit1, hit2);
  double shift = 0.5 * tof * kLightVelocity_cm_ps;

  TVector3 annihilationPoint;
  annihilationPoint.SetX(middleOfLOR.X() + shift * versorOnLOR.X());
  annihilationPoint.SetY(middleOfLOR.Y() + shift * versorOnLOR.Y());
  annihilationPoint.SetZ(middleOfLOR.Z() + shift * versorOnLOR.Z());
  return annihilationPoint;
}

/**
 * Calculating distance from the center of the decay plane
 */
double EventCategorizerTools::calculatePlaneCenterDistance(const JPetBaseHit& firstHit, const JPetBaseHit& secondHit, const JPetBaseHit& thirdHit)
{
  TVector3 crossProd = (secondHit.getPos() - firstHit.getPos()).Cross(thirdHit.getPos() - secondHit.getPos());
  double distCoef = -crossProd.X() * secondHit.getPosX() - crossProd.Y() * secondHit.getPosY() - crossProd.Z() * secondHit.getPosZ();
  if (crossProd.Mag() != 0)
  {
    return fabs(distCoef) / crossProd.Mag();
  }
  else
  {
    ERROR("One of the hit has zero position vector - unable to calculate distance from the center of the surface");
    return -1.;
  }
}

double EventCategorizerTools::calculatePlanePointDistance(
  const JPetBaseHit* firstHit, const JPetBaseHit* secondHit, const JPetBaseHit* thirdHit, const TVector3& decayPoint)
{
  TVector3 n_vector = ((secondHit->getPos() - firstHit->getPos()).Cross(thirdHit->getPos() - secondHit->getPos())).Unit();
  TVector3 v = firstHit->getPos() - decayPoint;
  return fabs(n_vector*v);
}

/**
 * @brief Calculation of an annihilation point based on positions of three hits.
 */
TVector3 EventCategorizerTools::calculateAnnihilationPoint(const JPetBaseHit& hit1, const JPetBaseHit& hit2, const JPetBaseHit& hit3)
{
  // Calculating norm vector for a surface created by 3 hits (vectors of their positions)
  TVector3 surfaceVec;
  surfaceVec.SetX((hit2.getPosY() - hit1.getPosY()) * (hit3.getPosZ() - hit1.getPosZ()) -
                  (hit2.getPosZ() - hit1.getPosZ()) * (hit3.getPosY() - hit1.getPosY()));
  surfaceVec.SetY((hit2.getPosZ() - hit1.getPosZ()) * (hit3.getPosX() - hit1.getPosX()) -
                  (hit2.getPosX() - hit1.getPosX()) * (hit3.getPosZ() - hit1.getPosZ()));
  surfaceVec.SetZ((hit2.getPosX() - hit1.getPosX()) * (hit3.getPosY() - hit1.getPosY()) -
                  (hit2.getPosY() - hit1.getPosY()) * (hit3.getPosX() - hit1.getPosX()));

  // Unitary perpendicular vector
  TVector3 perpVec(-surfaceVec.Y(), surfaceVec.X(), 0);
  perpVec = perpVec.Unit();

  double theta = -acos(surfaceVec.Z() / surfaceVec.Mag());

  // Defining rotation transformation to 2D plane and its reverse
  TVector3 rotX, rotY, rotZ, rotXr, rotYr, rotZr;
  rotX.SetX(cos(theta) + perpVec.X() * perpVec.X() * (1 - cos(theta)));
  rotX.SetY(perpVec.X() * perpVec.Y() * (1 - cos(theta)));
  rotX.SetZ(perpVec.Y() * sin(theta));
  rotY.SetX(perpVec.X() * perpVec.Y() * (1 - cos(theta)));
  rotY.SetY(cos(theta) + perpVec.Y() * perpVec.Y() * (1 - cos(theta)));
  rotY.SetZ(-perpVec.X() * sin(theta));
  rotZ.SetX(-perpVec.Y() * sin(theta));
  rotZ.SetY(perpVec.X() * sin(theta));
  rotZ.SetZ(cos(theta));
  rotXr.SetX(cos(-theta) + perpVec.X() * perpVec.X() * (1 - cos(-theta)));
  rotXr.SetY(perpVec.X() * perpVec.Y() * (1 - cos(-theta)));
  rotXr.SetZ(perpVec.Y() * sin(-theta));
  rotYr.SetX(perpVec.X() * perpVec.Y() * (1 - cos(-theta)));
  rotYr.SetY(cos(-theta) + perpVec.Y() * perpVec.Y() * (1 - cos(-theta)));
  rotYr.SetZ(-perpVec.X() * sin(-theta));
  rotZr.SetX(-perpVec.Y() * sin(-theta));
  rotZr.SetY(perpVec.X() * sin(-theta));
  rotZr.SetZ(cos(-theta));

  // Centers of circles transformed
  TVector3 p1(rotX * hit1.getPos(), rotY * hit1.getPos(), rotZ * hit1.getPos());
  TVector3 p2(rotX * hit2.getPos(), rotY * hit2.getPos(), rotZ * hit2.getPos());
  TVector3 p3(rotX * hit3.getPos(), rotY * hit3.getPos(), rotZ * hit3.getPos());

  // Time differences of hits registration
  double tdiff21 = hit2.getTime() - hit1.getTime();
  double tdiff31 = hit3.getTime() - hit1.getTime();

  TVector3 intersection = findIntersection(p1, p2, p3, tdiff21, tdiff31);

  // Transforming back found intersection by reverse rotation
  TVector3 annihilationPoint(rotXr * intersection, rotYr * intersection, rotZr * intersection);
  return annihilationPoint;
}

/**
 * Method for determining type of event for streaming - 3 gamma annihilation
 */
/*bool EventCategorizerTools::stream3Gamma(const JPetEvent& event, JPetStatistics& stats, bool saveHistos, double d3SlotthetaMin, double d3TimeDiff,
                                         double d3PlaneCenterDist, double maxScatter)
{
  if (event.getHits().size() < 3)
  {
    return false;
  }
  for (uint i = 0; i < event.getHits().size(); i++)
  {
    JPetBaseHit firstHit = event.getHits().at(i);

    for (uint j = i + 1; j < event.getHits().size(); j++)
    {
      JPetBaseHit secondHit = event.getHits().at(j);

      if (checkForScatter(firstHit, secondHit, stats, saveHistos, maxScatter))
      {
        continue;
      }

      for (uint k = j + 1; k < event.getHits().size(); k++)
      {
        JPetBaseHit thirdHit = event.getHits().at(k);

        if (checkForScatter(firstHit, thirdHit, stats, saveHistos, maxScatter))
        {
          continue;
        }

        if (checkForScatter(secondHit, thirdHit, stats, saveHistos, maxScatter))
        {
          continue;
        }

        vector<double> relativeAngles;
        relativeAngles.push_back(TMath::RadToDeg() * firstHit.getPos().Angle(secondHit.getPos()));
        relativeAngles.push_back(TMath::RadToDeg() * secondHit.getPos().Angle(thirdHit.getPos()));
        relativeAngles.push_back(TMath::RadToDeg() * thirdHit.getPos().Angle(firstHit.getPos()));
        sort(relativeAngles.begin(), relativeAngles.end());

        double transformedX = relativeAngles.at(1) + relativeAngles.at(0);
        double transformedY = relativeAngles.at(1) - relativeAngles.at(0);
        double timeDiff = fabs(thirdHit.getTime() - firstHit.getTime());
        double planeCenterDist = calculatePlaneCenterDistance(firstHit, secondHit, thirdHit);

        if (saveHistos)
        {
          stats.fillHistogram("stream3g_thetas", transformedX, transformedY);
          stats.fillHistogram("stream3g_plane_dist", planeCenterDist);
          stats.fillHistogram("stream3g_tdiff", timeDiff);
        }
        if (transformedX > d3SlotthetaMin && timeDiff < d3TimeDiff && planeCenterDist < d3PlaneCenterDist)
        {
          if (saveHistos)
          {
            TVector3 ap = calculateAnnihilationPoint(firstHit, secondHit, thirdHit);
            stats.fillHistogram("ap_yx", ap.Y(), ap.X());
            stats.fillHistogram("ap_zx", ap.Z(), ap.X());
            stats.fillHistogram("ap_zy", ap.Z(), ap.Y());
            stats.fillHistogram("ap_yx_zoom", ap.Y(), ap.X());
            stats.fillHistogram("ap_zx_zoom", ap.Z(), ap.X());
            stats.fillHistogram("ap_zy_zoom", ap.Z(), ap.Y());
          }
          return true;
        }
      }
    }
  }
  return false;
}*/

/**
 * Helper method for estimating anihilation point
 */
TVector3 EventCategorizerTools::findIntersection(TVector3 hit1Pos, TVector3 hit2Pos, TVector3 hit3Pos, double t21, double t31)
{
  double R21 = sqrt(pow(hit2Pos(0) - hit1Pos(0), 2) + pow(hit2Pos(1) - hit1Pos(1), 2));
  double R32 = sqrt(pow(hit3Pos(0) - hit2Pos(0), 2) + pow(hit3Pos(1) - hit2Pos(1), 2));
  double R13 = sqrt(pow(hit1Pos(0) - hit3Pos(0), 2) + pow(hit1Pos(1) - hit3Pos(1), 2));

  double TDiffTOR1 = 0.0;
  double TDiffTOR2 = t21;
  double TDiffTOR3 = t31;

  TDiffTOR2 = kLightVelocity_cm_ps * TDiffTOR2;
  TDiffTOR3 = kLightVelocity_cm_ps * TDiffTOR3;

  double R0 = 0.0;

  if (R0 < (R21 - TDiffTOR2) / 2.0)
  {
    R0 = (R21 - TDiffTOR2) / 2.0;
  }
  if (R0 < (R32 - TDiffTOR2 - TDiffTOR3) / 2.0)
  {
    R0 = (R32 - TDiffTOR2 - TDiffTOR3) / 2.0;
  }
  if (R0 < (R13 - TDiffTOR3) / 2.0)
  {
    R0 = (R13 - TDiffTOR3) / 2.0;
  }

  double R1 = 0.;
  double R2 = 0.;
  double R3 = 0.;

  vector<double> temp, temp2;
  vector<vector<double>> points;
  temp.push_back(0.0);
  temp.push_back(0.0);

  points.push_back(temp);
  points.push_back(temp);
  points.push_back(temp);
  points.push_back(temp);
  points.push_back(temp);
  points.push_back(temp);
  temp.clear();

  double Distance = 0.0;
  double MinDistance = 0.0;
  double PreviousDistance = 10000000.0;

  int test = 1;
  while (test)
  {
    R1 = TDiffTOR1 + R0 + 1;
    R2 = TDiffTOR2 + R0 + 1;
    R3 = TDiffTOR2 + R0 + 1;
    points = findIntersectiosOfCircles(hit1Pos, hit2Pos, hit3Pos, R1, R2, R3, R13, R21, R32);

    MinDistance = 1000000.0;
    for (unsigned i = 0; i < 2; i++)
    {
      for (unsigned j = 0; j < 2; j++)
      {
        for (unsigned k = 0; k < 2; k++)
        {
          Distance = sqrt(pow(points[i][0] - points[j + 2][0], 2) + pow(points[i][1] - points[j + 2][1], 2)) +
                     sqrt(pow(points[i][0] - points[k + 4][0], 2) + pow(points[i][1] - points[k + 4][1], 2)) +
                     sqrt(pow(points[k + 4][0] - points[j + 2][0], 2) + pow(points[k + 4][1] - points[j + 2][1], 2));
          if (Distance < MinDistance)
          {
            MinDistance = Distance;
            temp.clear();
            temp.push_back(points[i][0]);
            temp.push_back(points[i][1]);
            temp.push_back(points[2 + j][0]);
            temp.push_back(points[2 + j][1]);
            temp.push_back(points[4 + k][0]);
            temp.push_back(points[4 + k][1]);
          }
        }
      }
    }
    test++;
    if (test % 50 == 0)
    {
      if (MinDistance == 1000000.0)
      {
        temp.clear();
        temp.push_back(100.0);
        temp.push_back(100.0);
        temp.push_back(100.0);
        temp.push_back(100.0);
        temp.push_back(100.0);
        temp.push_back(100.0);
        break;
      }
    }
    if (MinDistance > PreviousDistance)
      test = 0;
    else
    {
      PreviousDistance = MinDistance;
      temp2 = temp;
    }
    R0 += 1;
  }
  vector<double> R0s, Distances;
  if (MinDistance != 1000000.0)
    test = 1;

  double MinnDistance = 1000000.0;
  while (test)
  {
    R1 = TDiffTOR1 + R0 + 1;
    R2 = TDiffTOR2 + R0 + 1;
    R3 = TDiffTOR2 + R0 + 1;
    points = findIntersectiosOfCircles(hit1Pos, hit2Pos, hit3Pos, R1, R2, R3, R13, R21, R32);

    MinDistance = 1000000.;
    for (unsigned i = 0; i < 2; i++)
    {
      for (unsigned j = 0; j < 2; j++)
      {
        for (unsigned k = 0; k < 2; k++)
        {
          Distance = sqrt(pow(points[i][0] - points[j + 2][0], 2) + pow(points[i][1] - points[j + 2][1], 2)) +
                     sqrt(pow(points[i][0] - points[k + 4][0], 2) + pow(points[i][1] - points[k + 4][1], 2)) +
                     sqrt(pow(points[k + 4][0] - points[j + 2][0], 2) + pow(points[k + 4][1] - points[j + 2][1], 2));
          if (Distance < MinDistance)
          {
            MinDistance = Distance;
            temp.clear();
            temp.push_back(points[i][0]);
            temp.push_back(points[i][1]);
            temp.push_back(points[2 + j][0]);
            temp.push_back(points[2 + j][1]);
            temp.push_back(points[4 + k][0]);
            temp.push_back(points[4 + k][1]);
          }
        }
      }
    }
    if (MinDistance != 1000000.0)
    {
      R0s.push_back(R0);
      Distances.push_back(MinDistance);
    }

    test++;
    if (test % 50 == 0)
    {
      if (MinDistance == 1000000.0)
      {
        temp.clear();
        temp.push_back(100.0);
        temp.push_back(100.0);
        temp.push_back(100.0);
        temp.push_back(100.0);
        temp.push_back(100.0);
        temp.push_back(100.0);
        break;
      }
      test = 0;
    }
    if (MinDistance < MinnDistance)
    {
      MinnDistance = MinDistance;
    }
    PreviousDistance = MinDistance;
    temp2 = temp;
    R0 -= 0.1;
  }

  if (MinnDistance != 1000000.0)
  {
    double R0Min;
    double minEle = *min_element(begin(Distances), end(Distances));
    if (minEle == Distances[0])
    {
      R0Min = R0s[0];
    }
    else if (minEle == Distances[Distances.size() - 1])
    {
      R0Min = R0s[R0s.size() - 1];
    }
    else
    {
      R0Min = findMinimumFromDerivative(R0s, Distances);
    }
    R1 = TDiffTOR1 + R0Min + 1;
    R2 = TDiffTOR2 + R0Min + 1;
    R3 = TDiffTOR2 + R0Min + 1;
    points = findIntersectiosOfCircles(hit1Pos, hit2Pos, hit3Pos, R1, R2, R3, R13, R21, R32);

    MinDistance = 1000000.0;
    for (unsigned i = 0; i < 2; i++)
    {
      for (unsigned j = 0; j < 2; j++)
      {
        for (unsigned k = 0; k < 2; k++)
        {
          Distance = sqrt(pow(points[i][0] - points[j + 2][0], 2) + pow(points[i][1] - points[j + 2][1], 2)) +
                     sqrt(pow(points[i][0] - points[k + 4][0], 2) + pow(points[i][1] - points[k + 4][1], 2)) +
                     sqrt(pow(points[k + 4][0] - points[j + 2][0], 2) + pow(points[k + 4][1] - points[j + 2][1], 2));
          if (Distance < MinDistance)
          {
            MinDistance = Distance;
            temp.clear();
            temp.push_back(points[i][0]);
            temp.push_back(points[i][1]);
            temp.push_back(points[2 + j][0]);
            temp.push_back(points[2 + j][1]);
            temp.push_back(points[4 + k][0]);
            temp.push_back(points[4 + k][1]);
          }
        }
      }
    }
  }

  TVector3 recoPoint((temp[0] + temp[2] + temp[4]) / 3, (temp[1] + temp[3] + temp[5]) / 3, hit1Pos(2));
  return recoPoint;
}

vector<vector<double>> EventCategorizerTools::findIntersectiosOfCircles(TVector3 hit1Pos, TVector3 hit2Pos, TVector3 hit3Pos, double R1, double R2,
                                                                        double R3, double R13, double R21, double R32)
{
  vector<vector<double>> points;
  vector<double> temp;
  temp.push_back(0);
  temp.push_back(0);
  points.push_back(temp);
  points.push_back(temp);
  points.push_back(temp);
  points.push_back(temp);
  points.push_back(temp);
  points.push_back(temp);

  points[0][0] =
      (hit1Pos(0) + hit2Pos(0)) / 2 + (pow(R1, 2) - pow(R2, 2)) * (hit2Pos(0) - hit1Pos(0)) / 2 / pow(R21, 2) +
      0.5 * (hit2Pos(1) - hit1Pos(1)) * sqrt(2 * (pow(R1, 2) + pow(R2, 2)) / pow(R21, 2) - pow(pow(R1, 2) - pow(R2, 2), 2) / pow(R21, 4) - 1);
  points[0][1] =
      (hit1Pos(1) + hit2Pos(1)) / 2 + (pow(R1, 2) - pow(R2, 2)) * (hit2Pos(1) - hit1Pos(1)) / 2 / pow(R21, 2) +
      0.5 * (hit1Pos(0) - hit2Pos(0)) * sqrt(2 * (pow(R1, 2) + pow(R2, 2)) / pow(R21, 2) - pow(pow(R1, 2) - pow(R2, 2), 2) / pow(R21, 4) - 1);
  points[1][0] =
      (hit1Pos(0) + hit2Pos(0)) / 2 + (pow(R1, 2) - pow(R2, 2)) * (hit2Pos(0) - hit1Pos(0)) / 2 / pow(R21, 2) -
      0.5 * (hit2Pos(1) - hit1Pos(1)) * sqrt(2 * (pow(R1, 2) + pow(R2, 2)) / pow(R21, 2) - pow(pow(R1, 2) - pow(R2, 2), 2) / pow(R21, 4) - 1);
  points[1][1] =
      (hit1Pos(1) + hit2Pos(1)) / 2 + (pow(R1, 2) - pow(R2, 2)) * (hit2Pos(1) - hit1Pos(1)) / 2 / pow(R21, 2) -
      0.5 * (hit1Pos(0) - hit2Pos(0)) * sqrt(2 * (pow(R1, 2) + pow(R2, 2)) / pow(R21, 2) - pow(pow(R1, 2) - pow(R2, 2), 2) / pow(R21, 4) - 1);

  points[2][0] =
      (hit2Pos(0) + hit3Pos(0)) / 2 + (pow(R2, 2) - pow(R3, 2)) * (hit3Pos(0) - hit2Pos(0)) / 2 / pow(R32, 2) +
      0.5 * (hit3Pos(1) - hit2Pos(1)) * sqrt(2 * (pow(R2, 2) + pow(R3, 2)) / pow(R32, 2) - pow(pow(R2, 2) - pow(R3, 2), 2) / pow(R32, 4) - 1);
  points[2][1] =
      (hit2Pos(1) + hit3Pos(1)) / 2 + (pow(R2, 2) - pow(R3, 2)) * (hit3Pos(1) - hit2Pos(1)) / 2 / pow(R32, 2) +
      0.5 * (hit2Pos(0) - hit3Pos(0)) * sqrt(2 * (pow(R2, 2) + pow(R3, 2)) / pow(R32, 2) - pow(pow(R2, 2) - pow(R3, 2), 2) / pow(R32, 4) - 1);
  points[3][0] =
      (hit2Pos(0) + hit3Pos(0)) / 2 + (pow(R2, 2) - pow(R3, 2)) * (hit3Pos(0) - hit2Pos(0)) / 2 / pow(R32, 2) -
      0.5 * (hit3Pos(1) - hit2Pos(1)) * sqrt(2 * (pow(R2, 2) + pow(R3, 2)) / pow(R32, 2) - pow(pow(R2, 2) - pow(R3, 2), 2) / pow(R32, 4) - 1);
  points[3][1] =
      (hit2Pos(1) + hit3Pos(1)) / 2 + (pow(R2, 2) - pow(R3, 2)) * (hit3Pos(1) - hit2Pos(1)) / 2 / pow(R32, 2) -
      0.5 * (hit2Pos(0) - hit3Pos(0)) * sqrt(2 * (pow(R2, 2) + pow(R3, 2)) / pow(R32, 2) - pow(pow(R2, 2) - pow(R3, 2), 2) / pow(R32, 4) - 1);

  points[4][0] =
      (hit1Pos(0) + hit3Pos(0)) / 2 + (pow(R3, 2) - pow(R1, 2)) * (hit1Pos(0) - hit3Pos(0)) / 2 / pow(R13, 2) +
      0.5 * (hit1Pos(1) - hit3Pos(1)) * sqrt(2 * (pow(R3, 2) + pow(R1, 2)) / pow(R13, 2) - pow(pow(R3, 2) - pow(R1, 2), 2) / pow(R13, 4) - 1);
  points[4][1] =
      (hit1Pos(1) + hit3Pos(1)) / 2 + (pow(R3, 2) - pow(R1, 2)) * (hit1Pos(1) - hit3Pos(1)) / 2 / pow(R13, 2) +
      0.5 * (hit3Pos(0) - hit1Pos(0)) * sqrt(2 * (pow(R3, 2) + pow(R1, 2)) / pow(R13, 2) - pow(pow(R3, 2) - pow(R1, 2), 2) / pow(R13, 4) - 1);
  points[5][0] =
      (hit1Pos(0) + hit3Pos(0)) / 2 + (pow(R3, 2) - pow(R1, 2)) * (hit1Pos(0) - hit3Pos(0)) / 2 / pow(R13, 2) -
      0.5 * (hit1Pos(1) - hit3Pos(1)) * sqrt(2 * (pow(R3, 2) + pow(R1, 2)) / pow(R13, 2) - pow(pow(R3, 2) - pow(R1, 2), 2) / pow(R13, 4) - 1);
  points[5][1] =
      (hit1Pos(1) + hit3Pos(1)) / 2 + (pow(R3, 2) - pow(R1, 2)) * (hit1Pos(1) - hit3Pos(1)) / 2 / pow(R13, 2) -
      0.5 * (hit3Pos(0) - hit1Pos(0)) * sqrt(2 * (pow(R3, 2) + pow(R1, 2)) / pow(R13, 2) - pow(pow(R3, 2) - pow(R1, 2), 2) / pow(R13, 4) - 1);

  return points;
}

double EventCategorizerTools::findMinimumFromDerivative(std::vector<double> x_vec, std::vector<double> y_vec)
{
  // Checking which element i of y values vecotr is a minimum, smaller than elements i-1 and i+1
  unsigned minIndex = 0;
  for (unsigned i = 1; i < y_vec.size() - 1; ++i)
  {
    if (y_vec.at(i) < y_vec.at(i - 1) && y_vec.at(i) < y_vec.at(i + 1))
    {
      minIndex = i;
      break;
    }
  }

  double a = (y_vec[minIndex + 1] - y_vec[minIndex] - (y_vec[minIndex] - y_vec[minIndex - 1])) /
             ((x_vec[minIndex + 1] + x_vec[minIndex]) / 2.0 - (x_vec[minIndex] + x_vec[minIndex - 1]) / 2.0);

  double b = y_vec[minIndex + 1] - y_vec[minIndex] - a * (x_vec[minIndex + 1] + x_vec[minIndex]) / 2.0;

  if (a > 0.0)
  {
    return -b / a;
  }
  else
  {
    return 0.0;
  }
}
