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
#include <algorithm>

using namespace std;

/**
* Method for determining type of event - back to back 2 gamma
*/
bool EventCategorizerTools::checkFor2Gamma(const JPetEvent& event, JPetStatistics& stats, bool saveHistos,
                          double b2bSlotThetaDiff, double b2bTimeDiff, double DeexTOTCutMin, double DeexTOTCutMax, 
                          double ToTCut2AnniMin, double ToTCut2AnniMax)
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
      if(firstHit.getEnergy() > ToTCut2AnniMax||secondHit.getEnergy() > ToTCut2AnniMax) continue;
      if(firstHit.getEnergy() < ToTCut2AnniMin||secondHit.getEnergy() < ToTCut2AnniMin) continue;
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
        stats.fillHistogram("2Gamma_ToT", firstHit.getEnergy());
        stats.fillHistogram("2Gamma_ToT", secondHit.getEnergy());
        stats.fillHistogram("2Gamma_XYpos", firstHit.getPosX(), firstHit.getPosY());
        stats.fillHistogram("2Gamma_XYpos", secondHit.getPosX(), secondHit.getPosY());
        stats.fillHistogram("2Gamma_Z_ToT", firstHit.getPosZ(),  firstHit.getEnergy());
	stats.fillHistogram("2Gamma_Z_ToT", secondHit.getPosZ(),  secondHit.getEnergy());
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
bool EventCategorizerTools::checkFor3Gamma(const JPetEvent& event, JPetStatistics& stats, bool saveHistos, double DeexTOTCutMin, 
                            double DeexTOTCutMax, double ToTCut3AnniMin, double ToTCut3AnniMax)
{
  if (event.getHits().size() != 4) return false;
  int count = 0;
  JPetHit promptHit;
  for (uint i = 0; i < event.getHits().size(); i++) {
    //Prompts
    if(event.getHits().at(i).getEnergy() > DeexTOTCutMin && event.getHits().at(i).getEnergy() < DeexTOTCutMax ){ 
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
	      if(firstHit.getEnergy() > ToTCut3AnniMax||secondHit.getEnergy() > ToTCut3AnniMax||thirdHit.getEnergy() > ToTCut3AnniMax) continue;
        if(firstHit.getEnergy() < ToTCut3AnniMin||secondHit.getEnergy() < ToTCut3AnniMin||thirdHit.getEnergy() < ToTCut3AnniMin) continue;

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
        stats.fillHistogram("Lifetime_3Gamma_zoom", lifetime);
      TVector3 annhilationPoint = calculateAnnihilationPoint(firstHit, secondHit);
      stats.fillHistogram("3Gamma_AnnihPoint_XY", annhilationPoint.X(), annhilationPoint.Y());
      stats.fillHistogram("3Gamma_AnnihPoint_ZX", annhilationPoint.Z(), annhilationPoint.X());
      stats.fillHistogram("3Gamma_AnnihPoint_ZY", annhilationPoint.Z(), annhilationPoint.Y());
      int isVertex = findVertex(firstHit, secondHit, thirdHit, stats);
      if(isVertex == 3){
        stats.fillHistogram("Lifetime_3Gamma_vertex_reconstruction", lifetime);
        stats.fillHistogram("Lifetime_3Gamma_vertex_reconstruction_zoom", lifetime);
      }
      if(isVertex == 2){
        stats.fillHistogram("Lifetime_3Gamma_vertex_cut_larger", lifetime);
        stats.fillHistogram("Lifetime_3Gamma_vertex_cut_larger_zoom", lifetime);
      }
      if(isVertex == 1) {
        stats.fillHistogram("Lifetime_3Gamma_vertex_cut", lifetime);
        stats.fillHistogram("Lifetime_3Gamma_vertex_cut_zoom", lifetime);
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
    double tot = HitFinderTools::calculateTOT(event.getHits().at(i), 
                                              HitFinderTools::getTOTCalculationType(fTOTCalculationType));
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
  std::string fTOTCalculationType, double ToTCutScattMax)
{
  stats.fillHistogram("Hit_multiplicity_scatt0", event.getHits().size()); // no cuts
  int scatters_count = 0;
  if(event.getHits().size() == 4){
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
      double scattTOF = calculateScatteringTime(primaryHit, scatterHit);
      double timeDiff = scatterHit.getTime() - primaryHit.getTime();
      double dist = calculateDistance(primaryHit,scatterHit);
      double M = (scattTOF - timeDiff);
      if(saveHistos){
        if (fabs(scattTOF - timeDiff) < scatterTOFTimeDiff && abs(primaryHit.getScintillator().getID() - scatterHit.getScintillator().getID()) <= 2){
          stats.fillHistogram("Scatt_Dist_vs_M_all", dist,M);      
          stats.fillHistogram("ScatterTOF_TimeDiff_all",M);
          stats.fillHistogram("Scatter_ToT_all", scatterHit.getEnergy());}
      }
      if(primaryHit.getEnergy() < ToTCutScattMax && scatterHit.getEnergy() < ToTCutScattMax){
      double scattAngle = calculateScatteringAngle(primaryHit, scatterHit);
      double dist_2D = calculateDistance2D(primaryHit,scatterHit);
      stats.fillHistogram("Scatt_Dist_vs_M_before_without_ID_cut", dist,M);
      if (abs(primaryHit.getScintillator().getID() - scatterHit.getScintillator().getID()) <= 2)
      {
        if (saveHistos) {
          stats.fillHistogram("ScatterTOF_TimeDiff_before", M);
	  stats.fillHistogram("ScatterAngle_PrimaryTOT_before", scattAngle, primaryHit.getEnergy());
  	  stats.fillHistogram("ScatterAngle_ScatterTOT_before", scattAngle, scatterHit.getEnergy());
	  stats.fillHistogram("Scatt_Dist_vs_M_before", dist,M);
          stats.fillHistogram("Scatt_TOF_vs_TimeDiff_before", scattTOF,timeDiff);
          stats.fillHistogram("Scatt_Dist_2D_vs_M_before", dist_2D,M);
        }

        if (fabs(scattTOF - timeDiff) < scatterTOFTimeDiff) {
          if (saveHistos) {
	    stats.fillHistogram("Scatter_ToT", scatterHit.getEnergy());
            stats.fillHistogram("Scatter_Z_ToT", scatterHit.getEnergy(), scatterHit.getPosZ());
            stats.fillHistogram("ScatterAngle_PrimaryTOT", scattAngle, primaryHit.getEnergy());
            stats.fillHistogram("ScatterAngle_ScatterTOT", scattAngle, scatterHit.getEnergy());
	    stats.fillHistogram("Scatt_Dist_vs_M", dist,M);
	    stats.fillHistogram("Scatt_TOF_vs_TimeDiff", scattTOF,timeDiff);
            stats.fillHistogram("Scatt_Dist_2D_vs_M", dist_2D,M);
            stats.fillHistogram("ScatterTOF_TimeDiff",M);
            stats.fillHistogram("Scatter_XYpos", scatterHit.getPosX(), scatterHit.getPosY());
	  }
          scatters_count++;
        }
      }
    }
    }
  }
  if(scatters_count != 0){
    double energy = 0;
    double energy_syn = 0;
    for (auto i = 0; i < event.getHits().size(); i++) {
      energy_syn+=event.getHits().at(i).getEnergy();
      energy+=HitFinderTools::calculateTOT(event.getHits().at(i), HitFinderTools::getTOTCalculationType(fTOTCalculationType));
    }
    stats.fillHistogram("Scatter_ToT_sum", energy);
    stats.fillHistogram("Scatter_ToT_sum_syn", energy);
    return true;
  } 
  }
  stats.fillHistogram("Hit_multiplicity_scatt1", event.getHits().size());
  return false;
}

bool EventCategorizerTools::checkForScatter_minimalDifference(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos, double scatterTOFTimeDiff, 
  std::string fTOTCalculationType, double ToTCutScattMax)
{
  JPetHit primaryHit, scatterHit;
  vector<double> TOF_TimeDiff = {};
  vector<vector<int>> ind_Primary_Scatter = {}; 
  if(event.getHits().size() == 4){
  double energy = 0;
  double energy_syn = 0;
  for (auto i = 0; i < event.getHits().size(); i++) {
  energy_syn+=event.getHits().at(i).getEnergy();
  energy+=HitFinderTools::calculateTOT(event.getHits().at(i), HitFinderTools::getTOTCalculationType(fTOTCalculationType));
    for (auto j = i + 1; j < event.getHits().size(); j++) {
      if (event.getHits().at(i).getTime() < event.getHits().at(j).getTime()) {
        primaryHit = event.getHits().at(i);
        scatterHit = event.getHits().at(j);
      } else {
        primaryHit = event.getHits().at(j);
        scatterHit = event.getHits().at(i);
      }
      if(fabs(primaryHit.getPosZ()) < 23 && fabs(scatterHit.getPosZ()) < 23 && primaryHit.getEnergy() < ToTCutScattMax && scatterHit.getEnergy() < ToTCutScattMax){
      if (abs(primaryHit.getScintillator().getID() - scatterHit.getScintillator().getID()) <= 2){
        double timeDiff = scatterHit.getTime() - primaryHit.getTime();
        double scattTOF = calculateScatteringTime(primaryHit, scatterHit);
        double M = fabs(scattTOF - timeDiff);
        TOF_TimeDiff.push_back(M);
        vector<int> indices = {i, j};
        ind_Primary_Scatter.push_back(indices);
      }
      }
    }
  }

  if(!TOF_TimeDiff.empty()){
    stats.fillHistogram("smallest_difference_Scatter_ToT_sum", energy);
    stats.fillHistogram("smallest_difference_Scatter_ToT_sum_syn", energy_syn);
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
    double dist_2D = calculateDistance2D(primaryHit,scatterHit);
    if (fabs(scattTOF - timeDiff) < scatterTOFTimeDiff){
      stats.fillHistogram("smallest_difference_Scatter_ToT", scatterHit.getEnergy());
      stats.fillHistogram("smallest_difference_Scatter_Z_ToT", scatterHit.getEnergy(), scatterHit.getPosZ());
      stats.fillHistogram("smallest_difference_ScatterAngle_PrimaryTOT", scattAngle, primaryHit.getEnergy());
      stats.fillHistogram("smallest_difference_ScatterAngle_ScatterTOT", scattAngle, scatterHit.getEnergy());
      stats.fillHistogram("smallest_difference_Scatt_Dist_vs_M", dist,M);
      stats.fillHistogram("smallest_difference_Scatt_TOF_vs_TimeDiff", scattTOF,timeDiff);
      stats.fillHistogram("smallest_difference_ScatterTOF_TimeDiff",M);
      stats.fillHistogram("smallest_difference_Scatter_XYpos", scatterHit.getPosX(), scatterHit.getPosY());  
      return true;
    }
  }
  }
  return false;
}


void EventCategorizerTools::HitTypeMC(std::vector<unsigned int> fHitType, bool& fIsScattered, bool& fContainsPrompt, bool& fIsPickOff, bool& fIsOPs, 
      bool& fIsSecondary,  std::vector<int>& fScatterIndex, std::vector<int>& fPromptIndex, std::vector<int>& fSignalIndex){
  int ind = 0;
      
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
}

int EventCategorizerTools::findVertex(const JPetHit& hit1, const JPetHit& hit2, const JPetHit& hit3, JPetStatistics& stats){
  Reconstructor * fReconstructor = Reconstructor::GetInstance();

  double scin_length = 500;
  scin_length *= 0.1; // in cm
  fReconstructor->setBarrelLength( scin_length );
  fReconstructor->setChamberRadius( 10.0 );

  fReconstructor->setGamma(0, hit1);
  fReconstructor->setGamma(1, hit2);
  fReconstructor->setGamma(2, hit3);

  int error = 0;
  int isVertex = 0;
  TVector3 sol[2];
  double t[2];
  error = fReconstructor->getSolution(sol[0], t[0], 0);
  error = fReconstructor->getSolution(sol[1], t[1], 1);
 
  double DOP = calculatePlanePointDistance(hit1, hit2, hit3, sol[1]);

  stats.fillHistogram("decay point XY all", sol[1].X(), sol[1].Y());
  stats.fillHistogram("decay point XZ all", sol[1].Z(), sol[1].X());
  stats.fillHistogram("decay point error all", error);
  stats.fillHistogram("decay point time all", t[1]);

  stats.fillHistogram("DOP all", DOP);

  if(error == 0){
  	stats.fillHistogram("decay point XY", sol[1].X(), sol[1].Y());
	  stats.fillHistogram("decay point XZ", sol[1].Z(), sol[1].X());
    stats.fillHistogram("decay point time", t[1]);
    stats.fillHistogram("decay point error", error);
    stats.fillHistogram("DOP", DOP);
    isVertex = 3;
    if(sol[1].Perp()<20. && fabs(sol[1].Z()) < 10.){
      stats.fillHistogram("decay point XY vertex cut larger", sol[1].X(), sol[1].Y());
      stats.fillHistogram("decay point XZ vertex cut larger", sol[1].Z(), sol[1].X());
      stats.fillHistogram("decay point time vertex cut larger", t[1]);
      stats.fillHistogram("DOP vertex cut", DOP);
      isVertex = 2;
    }
    if(sol[1].Perp()<10. && fabs(sol[1].Z()) < 5.){
      stats.fillHistogram("decay point XY vertex cut", sol[1].X(), sol[1].Y());
      stats.fillHistogram("decay point XZ vertex cut", sol[1].Z(), sol[1].X());
      stats.fillHistogram("decay point time vertex cut", t[1]);
      stats.fillHistogram("DOP vertex cut", DOP);
      isVertex = 1;
    }}

  return isVertex;

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

/**
 * * Calculating distance between point and plane for 3 hit
 * */
double EventCategorizerTools::calculatePlanePointDistance(
  const JPetHit& firstHit, const JPetHit& secondHit, const JPetHit& thirdHit, const TVector3& decayPoint)
{
  TVector3 n_vector = ((secondHit.getPos() - firstHit.getPos()).Cross(thirdHit.getPos() - secondHit.getPos())).Unit();
  TVector3 v = firstHit.getPos() - decayPoint;
  return fabs(n_vector*v);
}

