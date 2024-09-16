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
 *  @file EventCategorizerTools.cpp
 */

#include "EventCategorizerTools.h"
#include "HitFinderTools.h"
#include <TMath.h>
#include <vector>

using namespace std;

/**
* Method for determining type of event - back to back 2 gamma
*/
bool EventCategorizerTools::checkFor2Gamma(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos,
  double b2bSlotThetaDiff, double b2bTimeDiff
)
{
  if (event.getHits().size() < 2) {
    return false;
  }
  for (uint i = 0; i < event.getHits().size(); i++) {
    for (uint j = i + 1; j < event.getHits().size(); j++) {
      JPetHit firstHit, secondHit;
      if (event.getHits().at(i).getTime() < event.getHits().at(j).getTime()) {
        firstHit = event.getHits().at(i);
        secondHit = event.getHits().at(j);
      } else {
        firstHit = event.getHits().at(j);
        secondHit = event.getHits().at(i);
      }
      // Checking for back to back
      double timeDiff = fabs(firstHit.getTime() - secondHit.getTime());
      double deltaLor = (secondHit.getTime() - firstHit.getTime()) * kLightVelocity_cm_ps / 2.;
      double theta1 = min(firstHit.getBarrelSlot().getTheta(), secondHit.getBarrelSlot().getTheta());
      double theta2 = max(firstHit.getBarrelSlot().getTheta(), secondHit.getBarrelSlot().getTheta());
      double thetaDiff = min(theta2 - theta1, 360.0 - theta2 + theta1);
      if (saveHistos) {
        stats.fillHistogram("2Gamma_Zpos", firstHit.getPosZ());
        stats.fillHistogram("2Gamma_Zpos", secondHit.getPosZ());
        stats.fillHistogram("2Gamma_TimeDiff", timeDiff / 1000.0);
        stats.fillHistogram("2Gamma_DLOR", deltaLor);
        stats.fillHistogram("2Gamma_ThetaDiff", thetaDiff);
        stats.fillHistogram("2Gamma_Dist", calculateDistance(firstHit, secondHit));
        stats.fillHistogram("2Gamma_Energy", firstHit.getEnergy());
        stats.fillHistogram("2Gamma_Energy", secondHit.getEnergy());
        stats.fillHistogram("2Gamma_XYpos", firstHit.getPosX(), firstHit.getPosY());
        stats.fillHistogram("2Gamma_XYpos", secondHit.getPosX(), secondHit.getPosY());
        stats.fillHistogram("2Gamma_Z_Energy", firstHit.getPosZ(),  firstHit.getEnergy());
	stats.fillHistogram("2Gamma_Z_Energy", secondHit.getPosZ(),  secondHit.getEnergy());
      }
      if (fabs(thetaDiff - 180.0) < b2bSlotThetaDiff && timeDiff < b2bTimeDiff) {
        if (saveHistos) {
          TVector3 annhilationPoint = calculateAnnihilationPoint(firstHit, secondHit);
          stats.fillHistogram("Annih_TOF", calculateTOFByConvention(firstHit, secondHit));
          stats.fillHistogram("AnnihPoint_XY", annhilationPoint.X(), annhilationPoint.Y());
          stats.fillHistogram("AnnihPoint_ZX", annhilationPoint.Z(), annhilationPoint.X());
          stats.fillHistogram("AnnihPoint_ZY", annhilationPoint.Z(), annhilationPoint.Y());
          stats.fillHistogram("Annih_DLOR", deltaLor);
        }
        return true;
      }
    }
  }
  return false;
}

/**
* Method for determining type of event - 3Gamma
*/
bool EventCategorizerTools::checkFor3Gamma(const JPetEvent& event, JPetStatistics& stats, bool saveHistos)
{
if (event.getHits().size() != 4) return false;
  int count = 0;
  JPetHit promptHit;
  for (uint i = 0; i < event.getHits().size(); i++) {
    //Prompts
    if(event.getHits().at(i).getEnergy() > 360. && event.getHits().at(i).getEnergy() < 122000. ){ 
      promptHit = event.getHits().at(i);
      count++;}
  }
  auto hits = event.getHits();
  for (uint i = 0; i < event.getHits().size(); i++) {
    for (uint j = i + 1; j < event.getHits().size(); j++) {
      for (uint k = j + 1; k < event.getHits().size(); k++) {
        JPetHit firstHit = hits.at(i);
        JPetHit secondHit = hits.at(j);
        JPetHit thirdHit = hits.at(k);
	      if(firstHit.getEnergy() > 360.||secondHit.getEnergy() > 360.||thirdHit.getEnergy() > 360.) continue;
        if(firstHit.getEnergy() < 10.||secondHit.getEnergy() < 10.||thirdHit.getEnergy() < 10.) continue;

  vector<double> thetaAngles;
  thetaAngles.push_back(firstHit.getBarrelSlot().getTheta());
  thetaAngles.push_back(secondHit.getBarrelSlot().getTheta());
  thetaAngles.push_back(thirdHit.getBarrelSlot().getTheta());
  sort(thetaAngles.begin(), thetaAngles.end());

  vector<double> relativeAngles;
  relativeAngles.push_back(thetaAngles.at(1) - thetaAngles.at(0));
  relativeAngles.push_back(thetaAngles.at(2) - thetaAngles.at(1));
  relativeAngles.push_back(360.0 - thetaAngles.at(2) + thetaAngles.at(0));
  sort(relativeAngles.begin(), relativeAngles.end());
  double transformedX = relativeAngles.at(1) + relativeAngles.at(0);
  double transformedY = relativeAngles.at(1) - relativeAngles.at(0);

  TVector3 anniFirst(firstHit.getPos());
  TVector3 anniSecond(secondHit.getPos());
  TVector3 anniThird(thirdHit.getPos());
  Double_t delta12 = TMath::RadToDeg()*anniFirst.Angle(anniSecond);
  Double_t delta23 = TMath::RadToDeg()*anniSecond.Angle(anniThird);
  vector<Double_t> relativeAngles3D;
  relativeAngles3D.push_back(delta12);
  relativeAngles3D.push_back(delta23);
  relativeAngles3D.push_back(360. - (delta23+delta12));
  sort(relativeAngles3D.begin(), relativeAngles3D.end());

  double transformedX3D = relativeAngles3D.at(1) + relativeAngles3D.at(0);
  double transformedY3D = relativeAngles3D.at(1) - relativeAngles3D.at(0);

  if (saveHistos) {
    if(count == 1){
	stats.fillHistogram("3Gamma_Angles", transformedX, transformedY);
	stats.fillHistogram("3Gamma_Angles3D", transformedX3D, transformedY3D);
        if(transformedX3D>=190){
        stats.fillHistogram("3Gamma_Angles_after_cut", transformedX, transformedY);
        stats.fillHistogram("3Gamma_Angles3D_after_cut", transformedX3D, transformedY3D);
	      stats.fillHistogram("3Gamma_ToT", firstHit.getEnergy());
        stats.fillHistogram("3Gamma_ToT", secondHit.getEnergy());
        stats.fillHistogram("3Gamma_ToT", thirdHit.getEnergy());
        stats.fillHistogram("3Gamma_XYpos", firstHit.getPosX(), firstHit.getPosY());
        stats.fillHistogram("3Gamma_XYpos", secondHit.getPosX(), secondHit.getPosY());
        stats.fillHistogram("3Gamma_XYpos", thirdHit.getPosX(), thirdHit.getPosY());
        stats.fillHistogram("3Gamma_Z_ToT", firstHit.getEnergy(), firstHit.getPosZ());
        stats.fillHistogram("3Gamma_Z_ToT", secondHit.getEnergy(), secondHit.getPosZ());
        stats.fillHistogram("3Gamma_Z_ToT", thirdHit.getEnergy(), thirdHit.getPosZ());
	      double average_time = firstHit.getTime() + secondHit.getTime() + thirdHit.getTime();
        double lifetime = average_time/3. - promptHit.getTime(); 
        stats.fillHistogram("Lifetime_3Gamma", lifetime);
      bool isVertex = findVertex(firstHit, secondHit, thirdHit, stats);
      if(isVertex) {
        stats.fillHistogram("Lifetime_3Gamma_vertex_cut", lifetime);
        stats.fillHistogram("3Gamma_ToT_vertex_cut", firstHit.getEnergy());
        stats.fillHistogram("3Gamma_ToT_vertex_cut", secondHit.getEnergy());
        stats.fillHistogram("3Gamma_ToT_vertex_cut", thirdHit.getEnergy());}
      }}
      }
      }
    }
  }
  return true;
}

/**
* Method for determining type of event - prompt
*/
bool EventCategorizerTools::checkForPrompt(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos,
  double deexTOTCutMin, double deexTOTCutMax, std::string fTOTCalculationType)
{
  for (unsigned i = 0; i < event.getHits().size(); i++) {
    double tot = event.getHits().at(i).getEnergy();
    if (tot > deexTOTCutMin && tot < deexTOTCutMax) {
      if (saveHistos) {
        stats.fillHistogram("Deex_TOT_cut", tot);
      }
      return true;
    }
  }
  return false;
}

/**
* Method for determining type of event - scatter
*/
bool EventCategorizerTools::checkForScatter(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos, double scatterTOFTimeDiff, 
  std::string fTOTCalculationType, std::vector<bool> fHitInfo)
{
//HitInfo: fContainsPrompt, fIsOPs, fIsPickOff, fIsScattered, fIsSecondary
  
  stats.fillHistogram("Hit_multiplicity_scatt0", event.getHits().size()); // no cuts
  if(event.getHits().size() == 4){ 
  int scatters_count = 0;
  for (auto i = 0; i < event.getHits().size(); i++) {
    for (auto j = i + 1; j < event.getHits().size(); j++) {
     JPetHit primaryHit, scatterHit;
      if (event.getHits().at(i).getTime() < event.getHits().at(j).getTime()) {
        primaryHit = event.getHits().at(i);
        scatterHit = event.getHits().at(j);
      } else {
        primaryHit = event.getHits().at(j);
        scatterHit = event.getHits().at(i);
      }
      if(primaryHit.getEnergy() < 650. && scatterHit.getEnergy() < 650.){
      double scattAngle = calculateScatteringAngle(primaryHit, scatterHit);
      double scattTOF = calculateScatteringTime(primaryHit, scatterHit);
      double timeDiff = scatterHit.getTime() - primaryHit.getTime();
      double dist = calculateDistance(primaryHit,scatterHit);
      double M = (scattTOF - timeDiff);
      double dist2D = calculateDistance2D(primaryHit,scatterHit);
 

      if (saveHistos) {
        stats.fillHistogram("ScatterTOF_TimeDiff_before", scattTOF - timeDiff);
        stats.fillHistogram("ScatterAngle_PrimaryTOT_before", scattAngle, primaryHit.getEnergy());
        stats.fillHistogram("ScatterAngle_ScatterTOT_before", scattAngle, scatterHit.getEnergy());
        stats.fillHistogram("Scatt_Dist_vs_M_before", dist,M);
        stats.fillHistogram("Scatt_TOF_vs_TimeDiff_before", scattTOF,timeDiff);
              stats.fillHistogram("Scatt_Dist_2D_vs_M_before", dist2D,M);
        if(fHitInfo.at(0)&&fHitInfo.at(1)&&!fHitInfo.at(2)&&!fHitInfo.at(3)&&!fHitInfo.at(4)){
          stats.fillHistogram("Scatt_Dist_vs_M_before_signal", dist,M);
          stats.fillHistogram("Scatt_TOF_vs_TimeDiff_before_signal", scattTOF,timeDiff);}
        if(fHitInfo.at(3)){
          stats.fillHistogram("Scatt_Dist_vs_M_before_scatter", dist,M);
          stats.fillHistogram("Scatt_Dist_2D_vs_M_before_scatter", dist2D,M);
          stats.fillHistogram("Scatt_TOF_vs_TimeDiff_before_scatter", scattTOF,timeDiff);}
        if(fHitInfo.at(0)){
          stats.fillHistogram("Scatt_Dist_vs_M_before_prompt", dist,M);
          stats.fillHistogram("Scatt_TOF_vs_TimeDiff_before_prompt", scattTOF,timeDiff);}
      }

      if (fabs(scattTOF - timeDiff) < scatterTOFTimeDiff) {
        if (saveHistos) {
          stats.fillHistogram("Scattered_Energy", scatterHit.getEnergy());
          stats.fillHistogram("Scattered_Z_Energy", scatterHit.getEnergy(), scatterHit.getPosZ());
          stats.fillHistogram("ScatterAngle_PrimaryTOT", scattAngle, primaryHit.getEnergy());
          stats.fillHistogram("ScatterAngle_ScatterTOT", scattAngle, scatterHit.getEnergy());
          stats.fillHistogram("Scatt_Dist_vs_M", dist,M);
          stats.fillHistogram("Scatt_TOF_vs_TimeDiff", scattTOF,timeDiff);
          stats.fillHistogram("Scatt_Dist_2D_vs_M", dist2D,M);
	        stats.fillHistogram("ScatterTOF_TimeDiff",M);
	      }
	    scatters_count++;
      }
    }
    }
  }
  if(scatters_count) return true;
  }
  stats.fillHistogram("Hit_multiplicity_scatt1", event.getHits().size());
  return false;
}

bool EventCategorizerTools::checkForScatter_minimalDifference(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos, double scatterTOFTimeDiff, 
  std::string fTOTCalculationType, std::vector<bool> fHitInfo)
{
//HitInfo: fContainsPrompt, fIsOPs, fIsPickOff, fIsScattered, fIsSecondary
JPetHit primaryHit, scatterHit;
  vector<double> TOF_TimeDiff = {};
  vector<vector<int>> ind_Primary_Scatter = {}; 
  if (event.getHits().size() != 4) return false;
  for (auto i = 0; i < event.getHits().size(); i++) {
    for (auto j = i + 1; j < event.getHits().size(); j++) {
      if (event.getHits().at(i).getTime() < event.getHits().at(j).getTime()) {
        primaryHit = event.getHits().at(i);
        scatterHit = event.getHits().at(j);
      } else {
        primaryHit = event.getHits().at(j);
        scatterHit = event.getHits().at(i);
      }
      if(fabs(primaryHit.getPosZ()) < 23 && fabs(scatterHit.getPosZ()) < 23 && primaryHit.getEnergy() < 650. && scatterHit.getEnergy() < 650.){
        
        double timeDiff = scatterHit.getTime() - primaryHit.getTime();
        double scattTOF = calculateScatteringTime(primaryHit, scatterHit);
        double M = fabs(scattTOF - timeDiff);
        TOF_TimeDiff.push_back(M);
        ind_Primary_Scatter.push_back({static_cast<int>(i), static_cast<int>(j)});
      }
    }
  }

  if(!TOF_TimeDiff.empty()){
    std::vector<double>::iterator it = std::min_element(std::begin(TOF_TimeDiff), std::end(TOF_TimeDiff));
    int primary_ind = ind_Primary_Scatter.at(std::distance(std::begin(TOF_TimeDiff), it)).at(0);
    int scatter_ind = ind_Primary_Scatter.at(std::distance(std::begin(TOF_TimeDiff), it)).at(1);
    primaryHit = event.getHits().at(primary_ind);
    scatterHit = event.getHits().at(scatter_ind);
    double scattAngle = calculateScatteringAngle(primaryHit, scatterHit);
    double scattTOF = calculateScatteringTime(primaryHit, scatterHit);
    double timeDiff = scatterHit.getTime() - primaryHit.getTime();
    double dist = calculateDistance(primaryHit,scatterHit);
    double M = (scattTOF - timeDiff);
    if (fabs(scattTOF - timeDiff) < scatterTOFTimeDiff){
      stats.fillHistogram("smallest_difference_Scattered_Energy", scatterHit.getEnergy());
      stats.fillHistogram("smallest_difference_Scattered_Z_Energy", scatterHit.getEnergy(), scatterHit.getPosZ());
      stats.fillHistogram("smallest_difference_ScatterAngle_PrimaryTOT", scattAngle, primaryHit.getEnergy());
      stats.fillHistogram("smallest_difference_ScatterAngle_ScatterTOT", scattAngle, scatterHit.getEnergy());
      stats.fillHistogram("smallest_difference_Scatt_Dist_vs_M", dist,M);
      stats.fillHistogram("smallest_difference_Scatt_TOF_vs_TimeDiff", scattTOF,timeDiff);
      stats.fillHistogram("smallest_difference_ScatterTOF_TimeDiff",M);
      stats.fillHistogram("smallest_difference_Scattered_XYpos", scatterHit.getPosX(), scatterHit.getPosY());  
      return true;
    }
  }
  return false;
}

bool EventCategorizerTools::findVertex(const JPetHit& hit1, const JPetHit& hit2, const JPetHit& hit3, JPetStatistics& stats){
  Reconstructor * fReconstructor = Reconstructor::GetInstance();

  double scin_length = 500.;
  scin_length *= 0.1; // in cm
  fReconstructor->setBarrelLength( scin_length );
  fReconstructor->setChamberRadius( 10.0 );

  fReconstructor->setGamma(0, hit1);
  fReconstructor->setGamma(1, hit2);
  fReconstructor->setGamma(2, hit3);

  int error = 0;
  TVector3 sol[2];
  double t[2];
  error = fReconstructor->getSolution(sol[0], t[0], 0);
  error = fReconstructor->getSolution(sol[1], t[1], 1);
  
  int errFlag = 0;
  if(error == 1) errFlag = 1;
  if( t[1] < -20000. || t[1] > 20000.  ){
	  errFlag = 2;
  }else if( sol[1].Perp() > 50. || fabs( sol[1].Z() ) > scin_length/2. ){
	  errFlag = 3;
  }else if( TMath::IsNaN( t[1] ) ||
		TMath::IsNaN( sol[1].X() ) ||
		TMath::IsNaN( sol[1].Y() ) ||
		TMath::IsNaN( sol[1].Z() )){
	errFlag = 4;
  }

  stats.fillHistogram("decay point error all", error);
  stats.fillHistogram("decay point XY all", sol[1].X(), sol[1].Y());
  stats.fillHistogram("decay point XZ all", sol[1].Z(), sol[1].X());
  stats.fillHistogram("decay point time all", t[1]);

  if(errFlag == 0){
    stats.fillHistogram("decay point XY within detector flag", sol[1].X(), sol[1].Y());
    stats.fillHistogram("decay point XZ within detector flag", sol[1].Z(), sol[1].X());
    stats.fillHistogram("decay point time within detector flag", t[1]);
    stats.fillHistogram("decay point error within detector flag", error);}

  if( sol[1].Perp()<50 && fabs(sol[1].Z()) < 25. && fabs(t[1]) < 20000. && error != 1 ){
    stats.fillHistogram("decay point time within detector", t[1]);
    stats.fillHistogram("decay point XY within detector", sol[1].X(), sol[1].Y());
    stats.fillHistogram("decay point XZ within detector", sol[1].Z(), sol[1].X());
    stats.fillHistogram("decay point error within detector", error);
  }
  if(error == 0){
    stats.fillHistogram("decay point XY", sol[1].X(), sol[1].Y());
    stats.fillHistogram("decay point XZ", sol[1].Z(), sol[1].X());
    stats.fillHistogram("decay point time", t[1]);
    stats.fillHistogram("decay point error", error);
  return true; }

  return false;

}


bool EventCategorizerTools::checkForAccidental(const JPetEvent& event, JPetStatistics& stats, bool saveHistos, std::vector<unsigned int> fVtxIndex)
{
      if(std::equal(fVtxIndex.begin() + 1, fVtxIndex.end(), fVtxIndex.begin())) return false;
      else{
         for (uint i = 0; i < event.getHits().size(); i++) {
            JPetHit hit = event.getHits().at(i);
	    stats.fillHistogram("Acc_Energy", hit.getEnergy());
            stats.fillHistogram("Acc_XYpos", hit.getPosX(), hit.getPosY());
            stats.fillHistogram("Acc_Z_Energy", hit.getPosZ(),  hit.getEnergy());
          }

        return true;
      }
}


/**
* Calculation of distance between two hits
*/
double EventCategorizerTools::calculateDistance(const JPetHit& hit1, const JPetHit& hit2)
{
  return (hit1.getPos() - hit2.getPos()).Mag();
}

double EventCategorizerTools::calculateDistance2D(const JPetHit& hit1, const JPetHit& hit2)
{
  return sqrt((hit1.getPosX() - hit2.getPosX())*(hit1.getPosX() - hit2.getPosX())+(hit1.getPosY() - hit2.getPosY())*(hit1.getPosY() - hit2.getPosY()));
}

/**
* Calculation of time that light needs to travel the distance between primary gamma
* and scattered gamma. Return value in picoseconds.
*/
double EventCategorizerTools::calculateScatteringTime(const JPetHit& hit1, const JPetHit& hit2)
{
  return calculateDistance(hit1, hit2) / kLightVelocity_cm_ps;
}

/**
* Calculation of scatter angle between primary hit and scattered hit.
* This function assumes that source of first gamma was in (0,0,0).
* Angle is calculated from scalar product, return value in degrees.
*/
double EventCategorizerTools::calculateScatteringAngle(const JPetHit& hit1, const JPetHit& hit2)
{
  return TMath::RadToDeg() * hit1.getPos().Angle(hit2.getPos() - hit1.getPos());
}

/**
* Calculation point in 3D, where annihilation occured
*/
TVector3 EventCategorizerTools::calculateAnnihilationPoint(const JPetHit& hitA, const JPetHit& hitB)
{
  double tof = EventCategorizerTools::calculateTOF(hitA, hitB);
  return calculateAnnihilationPoint(hitA.getPos(), hitB.getPos(), tof);
}

TVector3 EventCategorizerTools::calculateAnnihilationPoint(const TVector3& hitA, const TVector3& hitB, double tof)
{
  TVector3 middleOfLOR = 0.5 * (hitA + hitB);
  TVector3 versorOnLOR = (hitB - hitA).Unit()  ;

  double shift = 0.5 * tof  * kLightVelocity_cm_ps;
  TVector3 annihilationPoint(middleOfLOR.X() + shift * versorOnLOR.X(),
                             middleOfLOR.Y() + shift * versorOnLOR.Y(),
                             middleOfLOR.Z() + shift * versorOnLOR.Z());
  return annihilationPoint;
}

double EventCategorizerTools::calculateTOFByConvention(const JPetHit& hitA, const JPetHit& hitB)
{
  if (hitA.getBarrelSlot().getTheta() < hitB.getBarrelSlot().getTheta()) {
    return calculateTOF(hitA, hitB);
  } else {
    return calculateTOF(hitB, hitA);
  }
}

double EventCategorizerTools::calculateTOF(const JPetHit& hitA, const JPetHit& hitB)
{
  return EventCategorizerTools::calculateTOF(hitA.getTime(), hitB.getTime());
}

double EventCategorizerTools::calculateTOF(double time1, double time2)
{
  return (time1 - time2);
}

/**
* Calculating distance from the center of the decay plane
*/
double EventCategorizerTools::calculatePlaneCenterDistance(
  const JPetHit& firstHit, const JPetHit& secondHit, const JPetHit& thirdHit)
{
  TVector3 crossProd = (secondHit.getPos() - firstHit.getPos()).Cross(thirdHit.getPos() - secondHit.getPos());
  double distCoef = -crossProd.X() * secondHit.getPosX() - crossProd.Y() * secondHit.getPosY() - crossProd.Z() * secondHit.getPosZ();
  if (crossProd.Mag() != 0) {
    return fabs(distCoef) / crossProd.Mag();
  } else {
    ERROR("One of the hit has zero position vector - unable to calculate distance from the center of the surface");
    return -1.;
  }
}
