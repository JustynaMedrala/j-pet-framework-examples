/**
 *  @copyright Copyright 2024 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  @file EventCategorizerTools.h
 */

#ifndef EVENTCATEGORIZERTOOLS_H
#define EVENTCATEGORIZERTOOLS_H

#include <Hits/JPetPhysRecoHit/JPetPhysRecoHit.h>
#include <Hits/JPetMCRecoHit/JPetMCRecoHit.h>
#include <JPetEvent/JPetEvent.h>
#include <JPetStatistics/JPetStatistics.h>
#include <TVector3.h>
#include <boost/property_tree/ptree.hpp>

static const double kLightVelocity_cm_ps = 0.0299792458;
static const double kUndefinedValue = 999.0;



/**
 * @brief Tools for Event Categorization
 *
 * Lots of tools in constatnt developement.
 */
class EventCategorizerTools
{
public:
  enum ScatterTestType
  {
    kSimpleParam,
    kMinMaxParams
  };

  // Categorizing methods
  static bool checkFor1Gamma(const JPetEvent& event, double totCutAnniMin, double totCutAnniMax, double totCutAnniMin_larger, 
                                          double totCutAnniMax_larger, double totCutDeexMin, double totCutDeexMax, JPetStatistics& stats, bool saveHistos);
  static bool checkFor2Gamma(const JPetEvent& event, JPetStatistics& stats, bool saveHistos, double maxThetaDiff, double maxTimeDiff,
                             double totCutAnniMin, double totCutAnniMax, const TVector3& sourcePos, ScatterTestType testType, double scatterTestValue,
                             double scatterTimeMin, double scatterTimeMax, double scatterAngleMin, double scatterAngleMax);

  static bool checkFor2Gamma(const JPetPhysRecoHit* firstHit, const JPetPhysRecoHit* secondHit, JPetStatistics& stats, bool saveHistos,
                             double maxThetaDiff, double maxTimeDiff, double totCutAnniMin, double totCutAnniMax, const TVector3& sourcePos, bool isScatter);

  static bool checkFor3Gamma(const JPetEvent& event, double minRelAngleCut, double maxTimeDiff, double totCutAnniMin, double totCutAnniMax, JPetStatistics& stats, bool saveHistos);

  static bool checkFor2GammaLifetime(const JPetEvent& event, std::vector<int> bad_ID, JPetStatistics& stats, bool saveHistos, double maxThetaDiff, double maxTimeDiff,
                                     double maxDOP, double totCutAnniMin, double totCutAnniMax, double totCutDeexMin, double totCutDeexMax,
                                     const TVector3& sourcePos, ScatterTestType testType, double scatterTestValue, double scatterTimeMin,
                                     double scatterTimeMax, double scatterAngleMin, double scatterAngleMax);

  template <typename HitType>
  static bool processHistograms(const std::vector<std::pair<const HitType*, const HitType*>>& annihilations, 
                                   JPetStatistics& stats, bool saveHistos, double maxThetaDiff, double maxTimeDiff, double maxDOP, const TVector3& sourcePos, ScatterTestType testType, double scatterTestValue, 
                                   double scatterTimeMin, double scatterTimeMax, double scatterAngleMin, double scatterAngleMax);

  template <typename HitType>
  static bool processHistograms(const std::vector<std::vector<const HitType*>>& annihilations, std::vector<std::pair<double, int>>& DOP_values,
                          JPetStatistics& stats, bool saveHistos, double minRelAngleCut, double minRelPhiCut, double minDistCut, double maxTimeDiff,  double maxDOP, const TVector3& sourcePos, ScatterTestType testType, double scatterTestValue, 
                          double scatterTimeMin, double scatterTimeMax, double scatterAngleMin, double scatterAngleMax);

  template <typename HitType>
  static void fillAnnihilationHistograms(const std::vector<std::pair<const HitType*, const HitType*>>& annihilations,
                                          const std::vector<const HitType*>& prompts, JPetStatistics& stats, 
                                          const TVector3& sourcePos, double totCutAnniMin, double totCutAnniMax);

  template <typename HitType>
  static void fillAnnihilationHistograms(const std::vector<std::vector<const HitType*>>& annihilations,
                                          const std::vector<const HitType*>& prompts, std::vector<std::pair<double, int>>& DOP_values, JPetStatistics& stats, 
                                          const TVector3& sourcePos, double totCutAnniMin, double totCutAnniMax);

  static bool checkFor3GammaLifetime(const JPetEvent& event, std::vector<int> bad_ID, double minRelAngleCut, double minRelPhiCut, double minDistCut, double maxTimeDiff, double maxDOP, JPetStatistics& stats, 
                            bool saveHistos, double totCutAnniMin, double totCutAnniMax, double totCutDeexMin, double totCutDeexMax, const TVector3& sourcePos, 
                            ScatterTestType testType, double scatterTestValue, double scatterTimeMin, double scatterTimeMax, double scatterAngleMin, 
                            double scatterAngleMax);

  // Single observable tests
  static void fillHitHistograms(const JPetEvent& event, JPetStatistics& stats, bool saveHistos, const std::string& histogramName);
  
  static void identifyAnnihilationHits(const JPetEvent& event, double totCutAnniMin, double totCutAnniMax, 
                               std::vector<std::pair<const JPetPhysRecoHit*, const JPetPhysRecoHit*>>& annihilations,
                               std::vector<std::pair<const JPetMCRecoHit*, const JPetMCRecoHit*>>& annihilationsMC);

  static void identifyAnnihilationHits(const JPetEvent& event, double totCutAnniMin, double totCutAnniMax, 
                               std::vector<std::vector<const JPetPhysRecoHit*>>& annihilations,
                               std::vector<std::vector<const JPetMCRecoHit*>>& annihilationsMC);

  static bool checkToT(const JPetPhysRecoHit* hit, double minToT, double maxToT);

  static bool checkToT(const JPetMCRecoHit* hit, double minToT, double maxToT);

  static double calculateToT(const JPetPhysRecoHit* hit);

  static double calculateToT(const JPetMCRecoHit* hit);

  static bool checkRelativeAngles(const TVector3& pos1, const TVector3& pos2, double maxThetaDiff);

  static bool checkForScatter(const JPetBaseHit* primaryHit, const JPetBaseHit* scatterHit, JPetStatistics& stats, bool saveHistos,
                              ScatterTestType testType, double scatterTestValue, double scatterTimeMin, double scatterTimeMax, double scatterAngleMin,
                              double scatterAngleMax);

  // Utilities
  static double calculateDistance(const JPetBaseHit* hit1, const JPetBaseHit* hit2);

  static double calculateDistance2D(const JPetBaseHit* hit1, const JPetBaseHit* hit2);

  static double calculateScatteringTime(const JPetBaseHit* hit1, const JPetBaseHit* hit2);

  static double calculateScatteringAngle(const JPetBaseHit* hit1, const JPetBaseHit* hit2);

  /// ToF is calculated as time1-time2.
  static double calculateTOF(const JPetBaseHit* hitA, const JPetBaseHit* hitB);

  static double calculateTOF(double time1, double time2);

  static double calculateTOF(const JPetBaseHit* hit);

  /// Tof calculated with the ordered hits with respect to scintillator number.
  /// The first one will be hit with smaller theta angle.
  /// See also: http://koza.if.uj.edu.pl/petwiki/index.php/Coordinate_system_in_Big_Barrel
  // cppcheck-suppress unusedFunction
  static double calculateTOFByConvention(const JPetBaseHit* hitA, const JPetBaseHit* hitB);

  // Annijilation point estimations
  static TVector3 calculateAnnihilationPoint(const JPetBaseHit* hit1, const JPetBaseHit* hit2);

  static double calculatePlaneCenterDistance(const JPetBaseHit& firstHit, const JPetBaseHit& secondHit, const JPetBaseHit& thirdHit);

  static double calculatePlanePointDistance(const JPetBaseHit* firstHit, const JPetBaseHit* secondHit, const JPetBaseHit* thirdHit, const TVector3& decayPoint);

  static TVector3 calculateAnnihilationPoint(const JPetBaseHit& hit1, const JPetBaseHit& hit2, const JPetBaseHit& hit3);

  // static bool stream2Gamma(const JPetEvent& event, JPetStatistics& stats, bool saveHistos, double b2bSlotThetaDiff, double b2bTimeDiff,
  //                          double maxScatter);

  // static bool stream3Gamma(const JPetEvent& event, JPetStatistics& stats, bool saveHistos, double d3SlotThetaMin, double d3TimeDiff,
  //                          double d3DistanceFromCenter, double maxScatter);

  static TVector3 findIntersection(TVector3 hit1Pos, TVector3 hit2Pos, TVector3 hit3Pos, double t21, double t31);

  static std::vector<std::vector<double>> findIntersectiosOfCircles(TVector3 Hit1Pos, TVector3 Hit2Pos, TVector3 Hit3Pos, double R1, double R2,
                                                                     double R3, double R13, double R21, double R32);

  static double findMinimumFromDerivative(std::vector<double> x_vec, std::vector<double> y_vec);
};

#endif /* !EVENTCATEGORIZERTOOLS_H */
