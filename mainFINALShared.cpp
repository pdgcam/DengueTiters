//
//  main.cpp
//  DengueTiters
//
//  Created by HSALJE on 7/15/16.
//  Copyright (c) 2016 hsalje. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <limits>
#include <random>
//#include <omp.h>
#include "main.h"
using namespace std;

//========== data ==========
int _numberOfSubject;
int _numberOfSymptomatic;
int *_numberOfTests;            //[subject]
int *_numberOfHITests;            //[subject]
int _totalTests;
int _totalHITests;
int *_vaccineSubject;                 //[subject]
int *_numberVaccines;                 //[subject]
int *_maxDayStudy;              //[subject]
int *_minDayStudy;              //[subject]
int **_dayOfTesting;            //[subject][test #]
int **_dayVacc;
int **_dayOfHITesting;            //[subject][test #]
double ***_resultOfTesting;        //[#ind][#tests][5]
double ***_resultOfHITesting;        //[#ind][#tests][5]
double ***_meanTiterindividual; //[#ind][#tests][5]
double ***_meanHITiterindividual; //[#ind][#tests][5]
double **_probabInfecHistoryIndividual; //[#ind][5]
int **_daySymptomaticSerKnown;
int **_daySymptomaticSerUnknown;
int **_serotypeSymptomatic;
int *_freqSymptomaticSerKnown;  //[subject]
int *_freqSymptomaticSerUnknown;//[subject]
int *_serostatusSubject;
int _minDay;
int _maxDay;
int _nSerAssay=1;
int _nSer=4;
int _sympOnlyEpiCurve;
int _sympOnly;
int _maxDelayRise;
int _linearDecay;
int _simThenEstimate;
double **_minDelaySinceInfection; //[subject][testno]
double **_minTimPreInfection; //[subject][testno]
int _numberOfSubjectSim;
double _normConstant;
int _addPreviousAugmentedInfs;
int _temporaryTitersAdd;
double _incubationPeriod; //Time from infection to time to rise
double _delayRise; //Time from symptom onset to rise
int _iterationPerSim;
double _minSlope;
double _maxErrorTiter;

//========== augmented data ==========
int *_numberAsymptomatic;       //[subject]
int **_serotypeAsymptomatic;    //[subject][# infections]
int **_dayAsymptomatic;         //[subject][# infections]
int _numberMissingSymptomaticSerotype;
int _numberNotMissingSymptomaticSerotype;
int **_augmentedSerotypeSymptomatic;    //[subject][# infections]
int *_idsMissingSymptomaticSerotype;
int *_ranksMissingSymptomaticSerotype;
double **_baselineTiter;         //[subject][5]
double **_IndividualImpactSlowDecay;         //[subject][5]
double **_baselineHITiter;         //[subject][5]
int *_baselineNaive;         //[subject]
double *_IndividualImpactTiter;         //[subject]
double ***_IndividualImpactTiterByInfec;         //[subject][# infections][ser]
double ***_IndividualImpactTempTiterByInfec;         //[subject][# infections][ser]
double ***_IndividualImpactSlopeByInfec;         //[subject][# infections][ser]
double ***_IndividualImpactTiterByInfecVac;         //[subject][3][ser]
double ***_IndividualImpactTempTiterByInfecVac;         //[subject][3][ser]
double ***_IndividualImpactSlopeByInfecVac;         //[subject][3][ser]
double *_IndividualImpactTempTiter;         //[subject]
double *_IndividualImpactSlope;         //[subject]
double *_IndividualImpactDelay;         //[subject]
int _nbIterRJMCMC;
int _nbIterRJMCMCDelay;
int _nbIterRJMCMCpostDelay;
int _nAsymptomatic;
double *_normFactor;    //[subject]
int _estNaiveBaseline;
int _minDayInfected;

//=========== parameters ============
int _numberOfParameters;
int _minTimeInfecPostVaccine;
int _seperateVaccineParameters;
int _individualLevelEffectsT;
int _individualLevelEffectsR;
int _individualLevelEffectsS;
int _individualLevelEffectsD;
int _individualLevelEffectsTByInfec;
int _individualLevelEffectsRByInfec;
int _individualLevelEffectsSByInfec;
int _individualLevelEffectsAllByInfec;
int _individualLevelEffectsAllByInfecVac;
int _individualLevelEffectsSlowDecay;
int _vaccineParOnly;
double *_parameter;
double *_initialParameter;
double *_rateForRandomWalk;
int _stopEstimatingPars;
int *_idOfSelectedParameter;
int _numberOfSelectedParameter;
int *_selectedParameter;
double *_probaTiterIndividual;
double _newLogLikGlobal;
double _logLikGlobal;
int maxTimeBetweenInfections;
int _maxDelaySympCalc;
int _maxDelayPreSympCalc;
int _calcMeanTiterImpact;
double _minTempTiterRise;
int _daysInMonth[12];
double _epiPDF[365*9];
double _epiCDF[365*9];
int _discreteTits;
double **_oldTiters;  //holders for changing these values
double **_oldTempTiters;
double **_oldSlopes;
double **_oldTitersVac;  //holders for changing these values
double **_oldTempTitersVac;
double **_oldSlopesVac;
double *_tmpTvec;
double *_tmpRvec;
double *_tmpSvec;
double *_tmpTvecVac;
double *_tmpRvecVac;
double *_tmpSvecVac;
double *_oldBaselines;
double *_oldSlowDecay;
double **_newIndEffectsT;
double **_newIndEffectsR;
double **_newIndEffectsS;


//=========== timing and serotype of asymptomatic infections ============

double drawBaselineTiter(int seronaive){
    if(seronaive==1){
        return 0;
    }
  double titer=0+runif()*10;
  return titer;
}

int drawPrimary(){
  int prim=int(runif()*2);
  return prim;
}

void updateEpiPDF(){
  int counter=0;
  int year, month, nday, day;
  
  double epiYear[365];
  for (day=0;day<365;day++){
    epiYear[day]=(1+_parameter[31]*(cos(_parameter[32]+2*3.14159265358979323846*day/365)))/365;
  }
  double tot=0;
  for (year=0; year<9; year++){tot+=_parameter[21+year];}
  for (year=0; year<9; year++){
    for (day=0;day<365;day++){
      _epiPDF[counter]=epiYear[day]*_parameter[21+year]/tot;
      if(counter==0){
        _epiCDF[counter]=_epiPDF[counter];}else{_epiCDF[counter]=_epiCDF[counter-1]+_epiPDF[counter];}
        counter++;
    }
  }
  
  tot=0;
  counter=0;
  for (year=0; year<9; year++){
    for (day=0;day<365;day++){
      tot+=_epiPDF[counter];
      counter++;
    }
  }
  
}






// Sort idx based on values in v
template <typename T>
void sort_indices_by_vector(const T v[], int size, int idx[]) {
  sort(idx, idx + size, [&v](int i1, int i2) { return v[i1] < v[i2]; });
}

double meanTiter(int symp, double delay, int infecting, int subject, int previous, int ser, int primary, int parity, int returnType=0, double indTeffect=1, double indReffect=1,double indSeffect=1,int vaccineResponse=0)
{
  double S,R,T,titer;
  T=indTeffect;
  R=indReffect;
  S=indSeffect;
  
  if(vaccineResponse==1){
      if(_individualLevelEffectsSlowDecay==1){
          T=T*exp(-_IndividualImpactSlowDecay[subject][ser]*delay);
          if(_linearDecay==1){T=max(0.0,T+(-_IndividualImpactSlowDecay[subject][ser]*delay));}
      }else{
          T=T*exp(-_parameter[3]*delay);
      }
    titer=R*exp(-S*delay)+T;
      if(_linearDecay==1){titer=max(0.0,R+(-S*delay))+T;}
    if (returnType==1) return T;
    if (returnType==2) return R;
    if (returnType==3) return S;
    return(titer);
  }
  
  if(infecting==1){
    T=T+_parameter[4];
    R=R+_parameter[4];
  }
  if(infecting==1 & primary==0 & _vaccineSubject[subject]==0){
    T=T+_parameter[5];
    R=R+_parameter[5];
  }
    
//Long term decay of T
    if(_individualLevelEffectsSlowDecay==1){
        T=T*exp(-_IndividualImpactSlowDecay[subject][ser]*delay);
        if(_linearDecay==1){T=max(0.0,T+(-_IndividualImpactSlowDecay[subject][ser]*delay));}
    }else{
        T=T*exp(-_parameter[3]*delay);
    }
    
  if (returnType==1) return T;
  if (returnType==2) return R;
  if (returnType==3) return S;
  
  titer=R*exp(-S*delay)+T;
    if(_linearDecay==1){titer=max(0.0,R+(-S*delay))+T;}
  
  return(titer);
}

//========== titers ==========
double computeTiterIndividualAtDay(int subject, int day, int ser)
{
  int infecting, noVaccines,pastSymp, infecNo, infecType,noInfectionsSerKnown,noInfectionsSerUnknown,noInfectionsAsymp,previousInfec,primary,bb;
  double delay, titerOut, titEff, titEffR,titEffS;
  noInfectionsSerKnown=_freqSymptomaticSerKnown[subject];
  noInfectionsSerUnknown=_freqSymptomaticSerUnknown[subject];
  noInfectionsAsymp=_numberAsymptomatic[subject];
  noVaccines=_numberVaccines[subject];
  int totInfVacc=noInfectionsSerKnown+noInfectionsSerUnknown+noInfectionsAsymp+noVaccines;
  int totInf=noInfectionsSerKnown+noInfectionsSerUnknown+noInfectionsAsymp;
  
  //Add in long-term decay of titers
    double aaa,bbb;
    delay=double(day)-double(_minDayStudy[subject]);
      aaa=_baselineTiter[subject][ser];
      bbb=_IndividualImpactSlowDecay[subject][ser];
      if(_individualLevelEffectsSlowDecay==1){titerOut=_baselineTiter[subject][ser]*exp(-_IndividualImpactSlowDecay[subject][ser]*delay);}else{
          titerOut=_baselineTiter[subject][ser]*exp(-_parameter[3]*delay);
      }
      if(_linearDecay==1&_individualLevelEffectsSlowDecay==1){titerOut=_baselineTiter[subject][ser]+(-_IndividualImpactSlowDecay[subject][ser]*delay);}
    if(titerOut<0){titerOut=0;}
  
  if (totInfVacc==0) return titerOut;
      
  double dates[totInf];
  double datesVacc[totInfVacc];
  int type[totInfVacc];
  int sers[totInfVacc];
  //    int infecNoByType[totInfVacc];
  int counter=0;
  double minDate=9999;
  for (infecNo=0;infecNo<noInfectionsSerKnown; infecNo++)
  {   dates[counter]=int(double(_daySymptomaticSerKnown[subject][infecNo])+_parameter[1]+_IndividualImpactDelay[subject]);
    datesVacc[counter]=dates[counter];
    type[counter]=0;
    sers[counter]=_serotypeSymptomatic[subject][infecNo];
    if(dates[counter]<minDate){minDate=dates[counter];}
    counter++;
  }
  for (infecNo=0;infecNo<noInfectionsSerUnknown; infecNo++)
  {   dates[counter]=int(double(_daySymptomaticSerUnknown[subject][infecNo])+_parameter[1]+_IndividualImpactDelay[subject]);
    datesVacc[counter]=dates[counter];
    type[counter]=1;
    sers[counter]=_augmentedSerotypeSymptomatic[subject][infecNo];
    if(dates[counter]<minDate){minDate=dates[counter];}
    counter++;
  }
  for (infecNo=0;infecNo<noInfectionsAsymp; infecNo++)
  {   dates[counter]=int(double(_dayAsymptomatic[subject][infecNo])+_parameter[1]);
    datesVacc[counter]=dates[counter];
    type[counter]=2;
    sers[counter]=_serotypeAsymptomatic[subject][infecNo];
    if(dates[counter]<minDate){minDate=dates[counter];}
    counter++;
  }
  int idx[totInf];
  int aa;
  for (aa=0;aa<totInf;aa++){idx[aa]=aa;}
  sort_indices_by_vector(dates,totInf, idx);
  for (infecNo=0;infecNo<noVaccines; infecNo++)
  {   datesVacc[counter]=int(double(_dayVacc[subject][infecNo]));
    type[counter]=3;
    sers[counter]=999;
    bb=datesVacc[counter];
    if(datesVacc[counter]<minDate){minDate=datesVacc[counter];}
    counter++;
  }
  int idxVacc[totInfVacc];
  for (aa=0;aa<totInfVacc;aa++){idxVacc[aa]=aa;}
  sort_indices_by_vector(datesVacc,totInfVacc, idxVacc);
  
  counter=0;
  for (aa=0;aa<totInfVacc;aa++){
    if(datesVacc[idxVacc[aa]]>day) break;
    counter++;
  }
  totInfVacc=counter;

    if (day<=minDate) return titerOut;

  double permTiter=titerOut;
  double day2;
  int vaccInd;
  int VaccIndicator=0;
  int InfectionIndicator=0;
  for (infecNo=0;infecNo<totInfVacc;infecNo++){
    infecType=type[idxVacc[infecNo]];
    day2=day;
    delay=double(day2)-double(datesVacc[idxVacc[infecNo]]);
    if(delay<=0)continue;
    if(_temporaryTitersAdd==0){titerOut=permTiter;}
    pastSymp=1;
    if(infecType==2){pastSymp=0;}
    primary=0;
    if(infecNo==0&_baselineNaive[subject]==1& _vaccineSubject[subject]==0){primary=1;}
    previousInfec=0;
    if(ser==sers[idxVacc[0]]){previousInfec=1;}
    infecting=0;
    if(ser==sers[idxVacc[infecNo]]){infecting=1;}
    vaccInd=0;
    if(infecType==3){
      vaccInd=1;
        if(_individualLevelEffectsAllByInfecVac==1){
          titEff=_IndividualImpactTiterByInfecVac[subject][VaccIndicator][ser];
          titEffR=_IndividualImpactTempTiterByInfecVac[subject][VaccIndicator][ser];
          titEffS=_IndividualImpactSlopeByInfecVac[subject][VaccIndicator][ser];
        }else{
            if(_baselineNaive[subject]==0){
                titEff=_parameter[34];
                titEffR=_parameter[36];
                titEffS=_parameter[38];}
            if(_baselineNaive[subject]==1){
                titEff=_parameter[44];
                titEffR=_parameter[46];
                titEffS=_parameter[48];}
        }
      VaccIndicator++;
    }
    if(infecType!=3){
        if(_individualLevelEffectsAllByInfec==1){
          titEff=_IndividualImpactTiterByInfec[subject][idx[InfectionIndicator]][ser];
          titEffR=_IndividualImpactTempTiterByInfec[subject][idx[InfectionIndicator]][ser];
          titEffS=_IndividualImpactSlopeByInfec[subject][idx[InfectionIndicator]][ser];
        }else{
            if(_vaccineSubject[subject]==0){
                titEff=_parameter[11];
                titEffR=_parameter[0];
                titEffS=_parameter[2];
            }
            if(_vaccineSubject[subject]==1&_baselineNaive[subject]==1){
                titEff=_parameter[54];
                titEffR=_parameter[56];
                titEffS=_parameter[2];
            }
            if(_vaccineSubject[subject]==1&_baselineNaive[subject]==0){
                titEff=_parameter[64];
                titEffR=_parameter[66];
                titEffS=_parameter[2];
            }
        }
      InfectionIndicator++;
    }
    
    titerOut+=meanTiter(pastSymp,delay,infecting,subject,previousInfec,ser,primary,infecNo+1,0,titEff, titEffR,titEffS,vaccInd);
    permTiter+=meanTiter(pastSymp,delay,infecting,subject,previousInfec,ser,primary,infecNo+1,1,titEff, titEffR,titEffS,vaccInd);
  }
  
  return(titerOut);
}


int checkMinRise(int subject,int includeMaxError=1){
    double titpre, titpost,serImpact,finalTiter,obsLogTiter, diffTiter;
    int noInfec, ser, i, minRise, day, noTests, testNo;
    minRise=0;
    int noInfection=0;
    double tmpRise;
    
    if(_vaccineSubject[subject]==1){
        titpre=computeTiterIndividualAtDay(subject,_dayVacc[subject][0]-1,0);
        titpost=computeTiterIndividualAtDay(subject,_dayVacc[subject][0]+365,0);
        if(titpost-titpre<_minTempTiterRise){minRise=1;}
    }
    
    
    noInfec=_freqSymptomaticSerKnown[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++){
          day=int(double(_daySymptomaticSerKnown[subject][i])+_parameter[1]+_IndividualImpactDelay[subject]);
          for (ser=0; ser<_nSerAssay;ser++){
              titpre=computeTiterIndividualAtDay(subject,day-1,ser);
              titpost=computeTiterIndividualAtDay(subject,day+365,ser);
              if(titpost-titpre<_minTempTiterRise){minRise=1;}
              
              tmpRise=_IndividualImpactTempTiterByInfec[subject][noInfection][ser]*exp(-_IndividualImpactSlopeByInfec[subject][noInfection][ser]*365);
              if(tmpRise<_minTempTiterRise){minRise=1;}
          }
          noInfection=noInfection+1;
      }
    }
    noInfec=_freqSymptomaticSerUnknown[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++){
          day=int(double(_daySymptomaticSerUnknown[subject][i])+_parameter[1]+_IndividualImpactDelay[subject]);
          for (ser=0; ser<_nSerAssay;ser++){
              titpre=computeTiterIndividualAtDay(subject,day-1,ser);
              titpost=computeTiterIndividualAtDay(subject,day+365,ser);
              if(titpost-titpre<_minTempTiterRise){minRise=1;}
              tmpRise=_IndividualImpactTempTiterByInfec[subject][noInfection][ser]*exp(-_IndividualImpactSlopeByInfec[subject][noInfection][ser]*365);
              if(tmpRise<_minTempTiterRise){minRise=1;}
          }
          noInfection=noInfection+1;
      }
    }
    noInfec=_numberAsymptomatic[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++){
          day=int(double(_dayAsymptomatic[subject][i])+_parameter[1]);
          for (ser=0; ser<_nSerAssay;ser++){
              titpre=computeTiterIndividualAtDay(subject,day-1,ser);
              titpost=computeTiterIndividualAtDay(subject,day+365,ser);
              if(titpost-titpre<_minTempTiterRise){minRise=1;}
              tmpRise=_IndividualImpactTempTiterByInfec[subject][noInfection][ser]*exp(-_IndividualImpactSlopeByInfec[subject][noInfection][ser]*365);
              if(tmpRise<_minTempTiterRise){minRise=1;}
          }
          noInfection=noInfection+1;
      }
    }
    
    //    If test within 2 days of 'titer rise day' don't include as we're not modelling the rate of rise
    int inf;
    int symp=_freqSymptomaticSerKnown[subject];
    int symp2=_freqSymptomaticSerUnknown[subject];
    int asymp=_numberAsymptomatic[subject];
    noTests=_numberOfTests[subject];
    int indicator,needInfec, interveninginfecs,currDaytest,prevDaytest;
    double prevTit,currTit;
    
    if(includeMaxError==1){
        int noBigErrors=0;
        noTests=_numberOfTests[subject];
        for (testNo=0;testNo<noTests;testNo++){
            indicator=0;
            needInfec=0;
            interveninginfecs=0;
            if(testNo>0){
                currTit=_resultOfTesting[subject][testNo][0];
                currDaytest=_dayOfTesting[subject][testNo];
                prevDaytest=_dayOfTesting[subject][testNo-1];
                if(currTit-prevTit>1.5&currDaytest-prevDaytest>300){
                    needInfec=1;
                }
            }
            prevTit=_resultOfTesting[subject][testNo][0];
            for (inf=0;inf<symp; inf++){
              day=_daySymptomaticSerKnown[subject][inf]+_parameter[1]+_IndividualImpactDelay[subject];
              if(_dayOfTesting[subject][testNo]>=day&_dayOfTesting[subject][testNo]<2+day){
                indicator=1;
                  break;
              };
                if(needInfec==1&testNo>0){
                    if(day-_parameter[1]-_IndividualImpactDelay[subject]>=_dayOfTesting[subject][testNo-1]&day-_parameter[1]-_IndividualImpactDelay[subject]<_dayOfTesting[subject][testNo]){
                        interveninginfecs=1;
                    }
                }
            }
            for (inf=0;inf<symp2;inf++){
              day=_daySymptomaticSerUnknown[subject][inf]+_parameter[1]+_IndividualImpactDelay[subject];
              if(_dayOfTesting[subject][testNo]>=day&_dayOfTesting[subject][testNo]<2+day){
                indicator=1;
                  break;
              }
                if(needInfec==1&testNo>0){
                    if(day-_parameter[1]-_IndividualImpactDelay[subject]>=_dayOfTesting[subject][testNo-1]&day-_parameter[1]-_IndividualImpactDelay[subject]<_dayOfTesting[subject][testNo]){
                        interveninginfecs=1;
                    }
                }
            }
            for (inf=0;inf<asymp;inf++){
                day=_dayAsymptomatic[subject][inf]+_parameter[1];
                if(needInfec==1&testNo>0){
                    if(day-_parameter[1]>=_dayOfTesting[subject][testNo-1]&day-_parameter[1]<_dayOfTesting[subject][testNo]){
                        interveninginfecs=1;
                    }
                }
            }
            if(needInfec==1&interveninginfecs==0&indicator==0){
                minRise=1;
                break;
            }
            if (indicator==1)continue;
            for (ser=0; ser<_nSerAssay; ser++)
            { finalTiter=_meanTiterindividual[subject][testNo][ser];
              serImpact=0;
              if(ser==1){serImpact=_parameter[16];}
              if(ser==2){serImpact=_parameter[17];}
              if(ser==3){serImpact=_parameter[18];}
              if(finalTiter==0){serImpact=0;}
              finalTiter=max(0.0,finalTiter+(serImpact));
              obsLogTiter=_resultOfTesting[subject][testNo][ser];
                diffTiter=abs(obsLogTiter-finalTiter);
                if(diffTiter>_maxErrorTiter){
                    noBigErrors+=1;
                }
            }
        }
        if(noBigErrors>1){minRise=1;}
    }
            
    
    return minRise;
}


void addAsymptomaticInfection(int subject,int newDate,int newSerotype)
{
  int oldNoInfections,i,counter,ser;
  oldNoInfections=_numberAsymptomatic[subject];
  int noSympInfec1=_freqSymptomaticSerKnown[subject];
  int noSympInfec2=_freqSymptomaticSerUnknown[subject];
  int newDates[oldNoInfections+1];
  int newSerotypes[oldNoInfections+1];
  
  counter=0;
  for (i=0;i<oldNoInfections+1;i++){
    if(i==0){
      newDates[i]=newDate;
      newSerotypes[i]=newSerotype;
    }else{
      newDates[i]=_dayAsymptomatic[subject][counter];               //vector with old dates
      newSerotypes[i]=_serotypeAsymptomatic[subject][counter];      //vector with old serotypes
      counter++;
    }
  }
    
  for (i=0;i<(oldNoInfections+1);i++){
    _dayAsymptomatic[subject][i]=newDates[i];
    _serotypeAsymptomatic[subject][i]=newSerotypes[i];
  }
  
  counter=0;

    if(_individualLevelEffectsAllByInfec==1){
  for (i=0;i<oldNoInfections+noSympInfec1+noSympInfec2+1;i++){
    for (ser=0; ser<_nSerAssay; ser++){
      if(i==noSympInfec1+noSympInfec2){
        _newIndEffectsT[i][ser]=_tmpTvec[ser];
        _newIndEffectsR[i][ser]=_tmpRvec[ser];
        _newIndEffectsS[i][ser]=_tmpSvec[ser];
      }else{
        _newIndEffectsT[i][ser]=_IndividualImpactTiterByInfec[subject][counter][ser];
        _newIndEffectsR[i][ser]=_IndividualImpactTempTiterByInfec[subject][counter][ser];
        _newIndEffectsS[i][ser]=_IndividualImpactSlopeByInfec[subject][counter][ser];
        if(ser==(_nSerAssay-1)){counter++;};
      }
    }
  }

  double aa;
  for (i=0;i<oldNoInfections+noSympInfec1+noSympInfec2+1;i++){
    for (ser=0;ser<_nSerAssay;ser++){
        aa=_newIndEffectsT[i][ser];
        aa=_IndividualImpactTiterByInfec[subject][i][ser];
        _IndividualImpactTiterByInfec[subject][i][ser]=_newIndEffectsT[i][ser];
      _IndividualImpactTempTiterByInfec[subject][i][ser]=_newIndEffectsR[i][ser];
      _IndividualImpactSlopeByInfec[subject][i][ser]=_newIndEffectsS[i][ser];
      aa=_IndividualImpactTempTiterByInfec[subject][i][ser];
    }
  }
    }
  
  _numberAsymptomatic[subject]++;
  _nAsymptomatic++;
}


void removeAsymptomaticInfection(int subject,int infectionNumber)
{   int oldNoInfections,i,ser;
  oldNoInfections=_numberAsymptomatic[subject];
  int noSympInfec1=_freqSymptomaticSerKnown[subject];
  int noSympInfec2=_freqSymptomaticSerUnknown[subject];
  int newDates[oldNoInfections-1];
  int newSerotypes[oldNoInfections-1];

  int counter=0;
  for (i=0;i<oldNoInfections;i++){
    if(i!=infectionNumber){
      newDates[counter]=_dayAsymptomatic[subject][i];               //vector with old dates
      newSerotypes[counter]=_serotypeAsymptomatic[subject][i];      //vector with old serotype
      counter++;
    }
  }
  counter=0;
    if(_individualLevelEffectsAllByInfec==1){
  for (i=0;i<oldNoInfections+noSympInfec1+noSympInfec2;i++){
    for (ser=0;ser<_nSerAssay;ser++){
      if(i!=infectionNumber+noSympInfec1+noSympInfec2){
        _newIndEffectsT[counter][ser]=_IndividualImpactTiterByInfec[subject][i][ser];
        _newIndEffectsR[counter][ser]=_IndividualImpactTempTiterByInfec[subject][i][ser];
        _newIndEffectsS[counter][ser]=_IndividualImpactSlopeByInfec[subject][i][ser];
        if(ser==_nSerAssay-1){counter++;};
      }
    }
  }

  for (i=0;i<(oldNoInfections+noSympInfec1+noSympInfec2-1);i++){
    for (ser=0;ser<_nSerAssay;ser++){
      _IndividualImpactTiterByInfec[subject][i][ser]=_newIndEffectsT[i][ser];
      _IndividualImpactTempTiterByInfec[subject][i][ser]=_newIndEffectsR[i][ser];
      _IndividualImpactSlopeByInfec[subject][i][ser]=_newIndEffectsS[i][ser];
    }
  }
    }
  
  if(oldNoInfections+noSympInfec1+noSympInfec2-1>0){
    for (i=0;i<(oldNoInfections-1);i++){
      _dayAsymptomatic[subject][i]=newDates[i];
      _serotypeAsymptomatic[subject][i]=newSerotypes[i];
    }
  }
  
  _numberAsymptomatic[subject]--;
  _nAsymptomatic--;
}

void addSymptomaticInfectionSerKnown(int subject,int newDate,int newSerotype)
{   int counter, oldNoInfections,i, ser;
  double aa,bb,cc;
  oldNoInfections=_freqSymptomaticSerKnown[subject];
  int noAsympInfec=_numberAsymptomatic[subject];
  int noSympInfec1=_freqSymptomaticSerKnown[subject];
  int noSympInfec2=_freqSymptomaticSerUnknown[subject];
  int newDates[oldNoInfections+1];
  int newSerotypes[oldNoInfections+1];
  
  counter=0;
  for (i=0;i<oldNoInfections+1;i++){
    if(i==0){
      newDates[i]=newDate;
      newSerotypes[i]=newSerotype;
    }else{
      newDates[i]=_daySymptomaticSerKnown[subject][counter];               //vector with old dates
      newSerotypes[i]=_serotypeSymptomatic[subject][counter];      //vector with old serotypes
      counter++;
    }
  }
  for (i=0;i<oldNoInfections+1;i++){
    _daySymptomaticSerKnown[subject][i]=newDates[i];
    _serotypeSymptomatic[subject][i]=newSerotypes[i];
  }
  
  counter=0;
  for (i=0;i<noAsympInfec+noSympInfec1+noSympInfec2+1;i++){
    for (ser=0;ser<_nSerAssay;ser++){
      if(i==0){
        _newIndEffectsR[i][ser]=_tmpRvec[ser];
        _newIndEffectsS[i][ser]=_tmpSvec[ser];
        _newIndEffectsT[i][ser]=_tmpTvec[ser];
      }else{
        aa=_newIndEffectsR[i][ser];
        bb=_IndividualImpactTempTiterByInfec[subject][counter][ser];
        cc=_newIndEffectsT[i][ser];
        _newIndEffectsR[i][ser]=_IndividualImpactTempTiterByInfec[subject][counter][ser];
        _newIndEffectsS[i][ser]=_IndividualImpactSlopeByInfec[subject][counter][ser];
        _newIndEffectsT[i][ser]=_IndividualImpactTiterByInfec[subject][counter][ser];
        if(ser==_nSerAssay-1){counter++;};
      }
    }
  }
  for (i=0;i<noAsympInfec+noSympInfec1+noSympInfec2+1;i++){
    for (ser=0;ser<_nSerAssay;ser++){
      _IndividualImpactTiterByInfec[subject][i][ser]=_newIndEffectsT[i][ser];
      _IndividualImpactTempTiterByInfec[subject][i][ser]=_newIndEffectsR[i][ser];
      _IndividualImpactSlopeByInfec[subject][i][ser]=_newIndEffectsS[i][ser];
    }
  }
  
  _freqSymptomaticSerKnown[subject]++;
  _numberNotMissingSymptomaticSerotype++;
}

void drawSer(int date, int *serOut, double *probOut)
{
  double samp=_nSer*runif();
  int ser=int(samp);
  (*serOut)=ser;
  (*probOut)=1;
  
}

int chooseSerotypeSubject(int subject, int date)
{
  int newSerotype=999;
  int noInfectionsSerKnown=_freqSymptomaticSerKnown[subject];
  int noInfectionsSerUnknown=_freqSymptomaticSerUnknown[subject];
  double probSer;
  int indic=0;
  int i;
  while(indic==0){
    int indicator=0;
    drawSer(date,&newSerotype,&probSer);
    if(noInfectionsSerKnown>0){
      for (i=0; i<noInfectionsSerKnown;i++)
        { if(_serotypeSymptomatic[subject][i]==newSerotype){indicator=1;}
        }
    }
    if(noInfectionsSerUnknown>0){
      for (i=0; i<noInfectionsSerUnknown;i++)
        {  if(_augmentedSerotypeSymptomatic[subject][i]==newSerotype){indicator=1;}
        }
    }
    if(indicator==0){indic=1;}
  }
  
  return newSerotype;
}

void addSymptomaticInfectionSerUnknown(int subject,int newDate, int Serotype)
{   int oldNoInfections,i,counter,ser;
  oldNoInfections=_freqSymptomaticSerUnknown[subject];
  int noAsympInfec=_numberAsymptomatic[subject];
  int noSympInfec1=_freqSymptomaticSerKnown[subject];
  int noSympInfec2=_freqSymptomaticSerUnknown[subject];
  int newDates[oldNoInfections+1];
  int newSerotypes[oldNoInfections+1];
  double *newIndEffectsT[noAsympInfec+noSympInfec1+noSympInfec2+1];
  double *newIndEffectsR[noAsympInfec+noSympInfec1+noSympInfec2+1];
  double *newIndEffectsS[noAsympInfec+noSympInfec1+noSympInfec2+1];
  for (i=0;i<noAsympInfec+noSympInfec1+noSympInfec2+1;i++){
    newIndEffectsT[i]=new double[_nSerAssay];
    newIndEffectsR[i]=new double[_nSerAssay];
    newIndEffectsS[i]=new double[_nSerAssay];
  }
  
  counter=0;
  for (i=0;i<oldNoInfections+1;i++){
    if(i==0){
      newDates[i]=newDate;
      newSerotypes[i]=Serotype;
    }else{
      newDates[i]=_daySymptomaticSerUnknown[subject][counter];               //vector with old dates
      newSerotypes[i]=_augmentedSerotypeSymptomatic[subject][counter];      //vector with old serotypes
      counter++;
    }
  }
  if(oldNoInfections>0){
    delete _augmentedSerotypeSymptomatic[subject];
    delete _daySymptomaticSerUnknown[subject];
  }
  
  _daySymptomaticSerUnknown[subject]=new int[oldNoInfections+1];
  _augmentedSerotypeSymptomatic[subject]=new int[oldNoInfections+1];
  
  for (i=0;i<oldNoInfections+1;i++){
    _daySymptomaticSerUnknown[subject][i]=newDates[i];
    _augmentedSerotypeSymptomatic[subject][i]=newSerotypes[i];
  }
  
  counter=0;
  for (i=0;i<noAsympInfec+noSympInfec1+noSympInfec2+1;i++){
    for (ser=0;ser<_nSerAssay;ser++){
      if(i==noSympInfec1){
        newIndEffectsR[i][ser]=_tmpRvec[ser];
        newIndEffectsS[i][ser]=_tmpSvec[ser];
        newIndEffectsT[i][ser]=_tmpTvec[ser];
      }else{
        newIndEffectsR[i][ser]=_IndividualImpactTempTiterByInfec[subject][counter][ser];
        newIndEffectsS[i][ser]=_IndividualImpactSlopeByInfec[subject][counter][ser];
        newIndEffectsT[i][ser]=_IndividualImpactTiterByInfec[subject][counter][ser];
        if(ser==_nSerAssay-1){counter++;};
      }
    }
  }
  delete[] _IndividualImpactTempTiterByInfec[subject];
  delete[] _IndividualImpactTiterByInfec[subject];
  delete[] _IndividualImpactSlopeByInfec[subject];
  _IndividualImpactTempTiterByInfec[subject]=new double*[noAsympInfec+noSympInfec1+noSympInfec2+1];
  _IndividualImpactSlopeByInfec[subject]=new double*[noAsympInfec+noSympInfec1+noSympInfec2+1];
  _IndividualImpactTiterByInfec[subject]=new double*[noAsympInfec+noSympInfec1+noSympInfec2+1];
  for (i=0;i<noAsympInfec+noSympInfec1+noSympInfec2+1;i++){
    _IndividualImpactTiterByInfec[subject][i]=new double[_nSerAssay];
    _IndividualImpactTempTiterByInfec[subject][i]=new double[_nSerAssay];
    _IndividualImpactSlopeByInfec[subject][i]=new double[_nSerAssay];
  }
  for (i=0;i<noAsympInfec+noSympInfec1+noSympInfec2+1;i++){
    for (ser=0;ser<_nSerAssay;ser++){
      _IndividualImpactTiterByInfec[subject][i][ser]=newIndEffectsT[i][ser];
      _IndividualImpactTempTiterByInfec[subject][i][ser]=newIndEffectsR[i][ser];
      _IndividualImpactSlopeByInfec[subject][i][ser]=newIndEffectsS[i][ser];
    }
  }
  
  int _newidsMissingSymptomaticSerotype[_numberMissingSymptomaticSerotype+1];
  int _newranksMissingSymptomaticSerotype[_numberMissingSymptomaticSerotype+1];
  
  for (i=0;i<_numberMissingSymptomaticSerotype;i++)
  {   _newidsMissingSymptomaticSerotype[i]=_idsMissingSymptomaticSerotype[i];
    _newranksMissingSymptomaticSerotype[i]=_ranksMissingSymptomaticSerotype[i];
  }
  _newidsMissingSymptomaticSerotype[_numberMissingSymptomaticSerotype]=subject;
  _newranksMissingSymptomaticSerotype[_numberMissingSymptomaticSerotype]=oldNoInfections;
  
  delete _idsMissingSymptomaticSerotype;
  delete _ranksMissingSymptomaticSerotype;
  
  _idsMissingSymptomaticSerotype=new int[_numberMissingSymptomaticSerotype+1];
  _ranksMissingSymptomaticSerotype=new int[_numberMissingSymptomaticSerotype+1];
  
  for (i=0;i<_numberMissingSymptomaticSerotype+1;i++)
  {   _idsMissingSymptomaticSerotype[i]=_newidsMissingSymptomaticSerotype[i];
    _ranksMissingSymptomaticSerotype[i]=_newranksMissingSymptomaticSerotype[i];
  }
  
  _freqSymptomaticSerUnknown[subject]++;
  _numberMissingSymptomaticSerotype++;
}



void drawDateAsymptomatic(int dayIn, int dayOut, int *dateOut, double *probOut)
{
  int day;
  double totProb=_epiCDF[dayOut-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]-_epiCDF[dayIn-_minDayInfected-int(_incubationPeriod)-int(_delayRise)];
  double aa,bb;
  aa=_epiCDF[dayOut-_minDayInfected];
  bb=_epiCDF[dayIn-_minDayInfected];
  int totDay=dayOut-dayIn;
  double u=runif();
  double cumprob=0;
  for (day=0;day<totDay;day++){
    cumprob+=_epiPDF[dayIn-_minDayInfected+day]/totProb;
    if(u<cumprob) break;
  }
  
  (*dateOut)=dayIn+day;
  double totProbEscape=-_parameter[12]*(_epiCDF[dayIn+day-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]-_epiCDF[dayIn-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]);
  
  (*probOut)=totProbEscape+log(1-exp(-_parameter[12]*_epiPDF[dayIn-_minDayInfected+day-int(_incubationPeriod)-int(_delayRise)]));
}


double probaSer(int date, int ser)
{
  return 1;
}


double probaSerOverall(int minDay, int maxDay, int ser)
{
  return (1.0/_nSer);
}

void drawSerOverall(int minDay, int maxDay,int *serOut, double *probOut)
{
  double samp=_nSer*runif();
  int ser=int(samp);
  (*serOut)=ser;
  (*probOut)=-9999;
}

void changeSymptomaticSerotype(int subject, int serotype, int infectionNumber)
{   _augmentedSerotypeSymptomatic[subject][infectionNumber]=serotype;
}



//========== titers ==========
double cumulativeInfectionProbability(int subject, int ser, int studyPeriodOnly)
{
  int infecting, pastSymp, infecNo, infecType,infecSer,minDay,maxDay,noInfectionsSerKnown,noInfectionsSerUnknown,noInfectionsAsymp;
  double probInf,timeRisk, delay, delay2;
  noInfectionsSerKnown=_freqSymptomaticSerKnown[subject];
  noInfectionsSerUnknown=_freqSymptomaticSerUnknown[subject];
  noInfectionsAsymp=_numberAsymptomatic[subject];
  int totInf=noInfectionsSerKnown+noInfectionsSerUnknown+noInfectionsAsymp;
  
  int aaa, bb;
  double dates[totInf];
  int type[totInf];
  int sers[totInf];
  int infecNoByType[totInf];
  int counter=0;
  for (infecNo=0;infecNo<noInfectionsSerKnown; infecNo++)
  {   dates[counter]=double(_daySymptomaticSerKnown[subject][infecNo])+_parameter[1]+_IndividualImpactDelay[subject];
    aaa=_IndividualImpactDelay[subject];
    bb=dates[counter];
    type[counter]=0;
    sers[counter]=_serotypeSymptomatic[subject][infecNo];
    infecNoByType[counter]=infecNo;
    counter++;
  }
  for (infecNo=0;infecNo<noInfectionsSerUnknown; infecNo++)
  {   dates[counter]=double(_daySymptomaticSerUnknown[subject][infecNo])+_parameter[1]+_IndividualImpactDelay[subject];
    type[counter]=1;
    sers[counter]=_augmentedSerotypeSymptomatic[subject][infecNo];
    infecNoByType[counter]=infecNo;
    counter++;
  }
  
  for (infecNo=0;infecNo<noInfectionsAsymp; infecNo++)
  {   dates[counter]=double(_dayAsymptomatic[subject][infecNo])+_parameter[1];
    type[counter]=2;
    sers[counter]=_serotypeAsymptomatic[subject][infecNo];
    infecNoByType[counter]=infecNo;
    counter++;
  }
  int idx[totInf];
  int aa;
  for (aa=0;aa<totInf;aa++){idx[aa]=aa;}
  
  sort_indices_by_vector(dates,totInf, idx);
  
  //Pre first infection
  probInf=0;
  maxDay=_maxDayStudy[subject]-1;
  minDay=_minDayStudy[subject]+1;
  if(totInf!=0){maxDay=dates[idx[0]];}
  
  probInf+=-_parameter[12]*(_epiCDF[maxDay-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]-_epiCDF[minDay-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]);
  
  double negative_infinity = - std::numeric_limits<double>::infinity();
  if(probInf==negative_infinity){
    
  }
  if (totInf==0) return probInf;
  
  //Impact of infections
  for (infecNo=0;infecNo<totInf;infecNo++){
    maxDay=_maxDayStudy[subject]-1;
    if(infecNo!=totInf-1){maxDay=dates[idx[infecNo+1]];}
    minDay=dates[idx[infecNo]];
    delay=double(maxDay-minDay);
    infecting=0;
    if(ser==sers[idx[infecNo]]){infecting=1;}
    
    if(infecting==1){
      if(type[idx[infecNo]]==2){
        probInf+=log(1-exp(-_parameter[12]*_epiPDF[minDay-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]));
      }else{probInf+=log(1-exp(-_parameter[12]*_epiPDF[minDay-_minDayInfected-int(_incubationPeriod)]));
      };
      if(probInf==negative_infinity){
        
      }
      break;
    }
    
    delay2=delay-maxTimeBetweenInfections;
    if(delay2<0){
      delay2=0;
      continue;}
    if(type[idx[infecNo]]==2){
      probInf+=-_parameter[12]*(_epiCDF[maxDay-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]-_epiCDF[minDay+maxTimeBetweenInfections-_minDayInfected-int(_incubationPeriod)]);}else{
        probInf+=-_parameter[12]*(_epiCDF[maxDay-_minDayInfected-int(_incubationPeriod)]-_epiCDF[minDay+maxTimeBetweenInfections-_minDayInfected-int(_incubationPeriod)]);
      }
  }
  return probInf;
}


void computeProbaInfectionHistoryIndividual(int subject, int studyPeriodOnly)
{
  int ser;
  for (ser=0;ser<_nSer;ser++){
    _probabInfecHistoryIndividual[subject][ser]=cumulativeInfectionProbability(subject, ser, studyPeriodOnly);
  }
}


void computeTiterIndividual(int subject)
{
  int testNo, ser, noTests;
  noTests=_numberOfTests[subject];
  for (testNo=0;testNo<noTests;testNo++){
    for (ser=0; ser<_nSerAssay; ser++){
      _meanTiterindividual[subject][testNo][ser]=computeTiterIndividualAtDay(subject, _dayOfTesting[subject][testNo],ser);
    }
  }
}


void maxDelaySympInfec()
{
  int subject, testNo,noTests, infecNo;
  double delay, timepre;
  int a;
  
  for (subject=0;subject<_numberOfSubject;subject++){
    noTests=_numberOfTests[subject];
    int noInfectionsSerKnown=_freqSymptomaticSerKnown[subject];
    int noInfectionsSerUnknown=_freqSymptomaticSerUnknown[subject];
    for (testNo=0;testNo<noTests;testNo++)
    {   _minDelaySinceInfection[subject][testNo]=9999;
      _minTimPreInfection[subject][testNo]=9999;
      for (infecNo=0;infecNo<noInfectionsSerKnown;infecNo++)
      {
        delay=_dayOfTesting[subject][testNo]-_daySymptomaticSerKnown[subject][infecNo];
        if(delay<_minDelaySinceInfection[subject][testNo]&delay>0){
          _minDelaySinceInfection[subject][testNo]=delay;
        }
        if(-delay<_minTimPreInfection[subject][testNo]&delay<0){
          _minTimPreInfection[subject][testNo]=-delay;
          a=_minTimPreInfection[subject][testNo];
        }
      }
      for (infecNo=0;infecNo<noInfectionsSerUnknown;infecNo++)
      {
        delay=_dayOfTesting[subject][testNo]-_daySymptomaticSerUnknown[subject][infecNo];
        if(delay<_minDelaySinceInfection[subject][testNo]&delay>0){
          _minDelaySinceInfection[subject][testNo]=delay;
        }
        if(-delay<_minTimPreInfection[subject][testNo]&delay<0){
          _minTimPreInfection[subject][testNo]=-delay;
        }
      }
      
    }
  }
}

void computeProbaIndividual(int subject, int sympOnly, int studyPeriodOnly, int vaccineBaselineOnly)
{
  computeTiterIndividual(subject);
  computeProbaInfectionHistoryIndividual(subject, studyPeriodOnly);
  
    double negative_infinity = - std::numeric_limits<double>::infinity();
    double infinity = std::numeric_limits<double>::infinity();
    
  double probaInd, finalTiter, obsLogTiter;
  int noTests, ser, testNo,i,b;
  int symp=_freqSymptomaticSerKnown[subject];
  int symp2=_freqSymptomaticSerUnknown[subject];
  int asymp=_numberAsymptomatic[subject];
  noTests=_numberOfTests[subject];
  probaInd=0;
  int indicatorVaccineAll=0;
  int indicatorAll=0;
  if(symp+symp2==0&sympOnly==1){indicatorAll=1;}

   
  int TotInf=symp+symp2+asymp;
  int inf, inf2;
    double firstInfec=9999;
  int indicator, indicatorVaccine;
  for (testNo=0;testNo<noTests;testNo++){
     
      if(vaccineBaselineOnly==1){
          
      }
      
      
    if(_dayOfTesting[subject][testNo]<_minDayStudy[subject])continue;
    b=_dayOfTesting[subject][testNo];
    indicator=0;
    indicatorVaccine=0;
    if(symp+symp2+asymp>0&vaccineBaselineOnly==1&_vaccineSubject[subject]==1){
      for (inf=0;inf<asymp; inf++){
          if(_dayAsymptomatic[subject][inf]<firstInfec){
              firstInfec=_dayAsymptomatic[subject][inf];
          }
          if(_dayOfTesting[subject][testNo]>_dayAsymptomatic[subject][inf]){
                  indicatorVaccine=1;
          }
      }
      for (inf=0;inf<symp; inf++){
          if(_daySymptomaticSerKnown[subject][inf]<firstInfec){
              firstInfec=_daySymptomaticSerKnown[subject][inf];
          }
          if(_dayOfTesting[subject][testNo]>_daySymptomaticSerKnown[subject][inf]){
              indicatorVaccine=1;
          }
      }
      for (inf=0;inf<symp2; inf++){
          if(_daySymptomaticSerUnknown[subject][inf]<firstInfec){
              firstInfec=_daySymptomaticSerUnknown[subject][inf];
          }
          if(_dayOfTesting[subject][testNo]>_daySymptomaticSerUnknown[subject][inf]){
              indicatorVaccine=1;
          }
      }
    }
      if(indicatorVaccine==1)break;
      
      int maxTest=noTests;
      int testNo2;
      
    if(symp+symp2>0&sympOnly==1){
      for (inf=0;inf<asymp; inf++){
        for (inf2=0;inf2<symp;inf2++){
          if(_dayOfTesting[subject][testNo]>_dayAsymptomatic[subject][inf]&_dayAsymptomatic[subject][inf]>_daySymptomaticSerKnown[subject][inf2]){
            indicator=1;
              break;
          };
        }
        for (inf2=0;inf2<symp2;inf2++){
          if(_dayOfTesting[subject][testNo]>_dayAsymptomatic[subject][inf]&_dayAsymptomatic[subject][inf]>_daySymptomaticSerUnknown[subject][inf2]){
            indicator=1;
              break;
          };
        }
      }
    }
      if(indicator==1)break;
      
      //There is a test that occurs XX days after an infection (without an intervening infection)
      if(symp+symp2>0&sympOnly==1){
          for (testNo2=0;testNo2<noTests;testNo2++){
              for (inf=0;inf<asymp; inf++){
                  for (inf2=0;inf2<symp;inf2++){
                      if(_dayOfTesting[subject][testNo2]>_dayAsymptomatic[subject][inf]&_dayAsymptomatic[subject][inf]>_daySymptomaticSerKnown[subject][inf2]&_dayOfTesting[subject][testNo2]-_daySymptomaticSerKnown[subject][inf2]<_maxDelaySympCalc){
                          indicator=1;
                          break;
                      };
                  }
                  for (inf2=0;inf2<symp2;inf2++){
                      if(_dayOfTesting[subject][testNo2]>_dayAsymptomatic[subject][inf]&_dayAsymptomatic[subject][inf]>_daySymptomaticSerUnknown[subject][inf2]&_dayOfTesting[subject][testNo2]-_daySymptomaticSerUnknown[subject][inf2]<_maxDelaySympCalc){
                          indicator=1;
                          break;
                      };
                  }
              }
          }
      }
      if(indicator==1)break;
      
    
//    If test within 2 days of 'titer rise day' don't include as we're not modelling the rate of rise
    indicator=0;
    double day;
    for (inf=0;inf<symp; inf++){
      day=_daySymptomaticSerKnown[subject][inf]+_parameter[1]+_IndividualImpactDelay[subject];
      if(_dayOfTesting[subject][testNo]>=day&_dayOfTesting[subject][testNo]<2+day){
        indicator=1;
      };
    }
    for (inf=0;inf<symp2;inf++){
      day=_daySymptomaticSerUnknown[subject][inf]+_parameter[1]+_IndividualImpactDelay[subject];
      if(_dayOfTesting[subject][testNo]>=day&_dayOfTesting[subject][testNo]<2+day){
        indicator=1;
      }
    }
    if(indicator==1)continue;

      
    double upperBound,lowerBound,diff;
    double serImpact;
    for (ser=0; ser<_nSerAssay; ser++)
    {   finalTiter=_meanTiterindividual[subject][testNo][ser];
      serImpact=0;
      if(ser==1){serImpact=_parameter[16];}
      if(ser==2){serImpact=_parameter[17];}
      if(ser==3){serImpact=_parameter[18];}
      if(finalTiter==0){serImpact=0;}
      finalTiter=max(0.0,finalTiter+(serImpact));
      obsLogTiter=_resultOfTesting[subject][testNo][ser];
      if(_discreteTits==1){
        //                if(obsLogTiter>10){
        //                    lowerBound=pnorm((10-finalTiter)/_parameter[8]);
        //                    diff=max(0.00000001,1-lowerBound);
        //                }else
        if(obsLogTiter<log(10.0)){
          upperBound=pnorm((log(10.0)-finalTiter)/_parameter[8]);
          diff=max(1e-10,upperBound);}else{
            //                        upperBound=pnorm((floor(obsLogTiter)+1-finalTiter)/_parameter[8]);
            //                        lowerBound=pnorm((floor(obsLogTiter)-finalTiter)/_parameter[8]);
            //
            //                        diff=max(0.00000001,upperBound-lowerBound);
            diff=max(1e-10,dnorm(obsLogTiter,finalTiter,_parameter[8]));
          }
          probaInd+=log(diff);
      }else{
        probaInd+=log(max(1e-10,dnorm(obsLogTiter,finalTiter,_parameter[8])));
      }

      if(probaInd==infinity|probaInd==negative_infinity){
        
      }
    }
  }
    if(indicatorVaccineAll==0&indicatorAll==0){
        for (ser=0; ser<_nSer; ser++){
          probaInd+=_probabInfecHistoryIndividual[subject][ser];
        }
    }

    if(probaInd==infinity|probaInd==negative_infinity){
    
    }

    double aa;
    //Add individual level impact - slow decay
    if(_individualLevelEffectsSlowDecay==1){
        for (ser=0; ser<_nSerAssay; ser++){
            double par1=_parameter[3]*_parameter[3]/_parameter[19];
            double par2=_parameter[3]/_parameter[19];
            aa=_IndividualImpactSlowDecay[subject][ser];
            probaInd+=log(dgamma(_IndividualImpactSlowDecay[subject][ser],par1,par2));
            if(probaInd==infinity|probaInd==negative_infinity){
                
            }
        }
    }
    
  //Add indivudal level impact - infections
  if(((sympOnly==0)|(symp2+symp>0))&_individualLevelEffectsAllByInfec==1&vaccineBaselineOnly==0&indicatorAll==0){
    
      if(_vaccineSubject[subject]==0){
        for (i=0;i<TotInf;i++){
            for (ser=0; ser<_nSerAssay; ser++){
                double par1=_parameter[11]*_parameter[11]/_parameter[13];
                double par2=_parameter[11]/_parameter[13];
                probaInd+=log(dgamma(_IndividualImpactTiterByInfec[subject][i][ser],par1,par2));

                par1=_parameter[0]*_parameter[0]/_parameter[14];
                par2=_parameter[0]/_parameter[14];
                probaInd+=log(dgamma(_IndividualImpactTempTiterByInfec[subject][i][ser],par1,par2));

                par1=_parameter[2]*_parameter[2]/_parameter[15];
                par2=_parameter[2]/_parameter[15];
                probaInd+=log(dgamma(_IndividualImpactSlopeByInfec[subject][i][ser],par1,par2));
            }
        }
      }
      if(_vaccineSubject[subject]==1&_baselineNaive[subject]==1){
          for (i=0;i<TotInf;i++){
              for (ser=0; ser<_nSerAssay; ser++){
                  double mean=_parameter[54];//*_parameter[11];
                  double par1=mean*mean/_parameter[55];
                  double par2=mean/_parameter[55];
                  probaInd+=log(dgamma(_IndividualImpactTiterByInfec[subject][i][ser],par1,par2));
                  
                  mean=_parameter[56];//*_parameter[0];
                  par1=mean*mean/_parameter[57];
                  par2=mean/_parameter[57];
                  probaInd+=log(dgamma(_IndividualImpactTempTiterByInfec[subject][i][ser],par1,par2));
                  
                  mean=_parameter[58];//*_parameter[2];
                  par1=mean*mean/_parameter[59];
                  par2=mean/_parameter[59];
                  probaInd+=log(dgamma(_IndividualImpactSlopeByInfec[subject][i][ser],par1,par2));
              }
          }
      }
      if(_vaccineSubject[subject]==1&_baselineNaive[subject]==0){
          for (i=0;i<TotInf;i++){
              for (ser=0; ser<_nSerAssay; ser++){
                  double mean=_parameter[64];//*_parameter[11];
                  double par1=mean*mean/_parameter[65];
                  double par2=mean/_parameter[65];
                  probaInd+=log(dgamma(_IndividualImpactTiterByInfec[subject][i][ser],par1,par2));
                  
                  mean=_parameter[66];//*_parameter[0];
                  par1=mean*mean/_parameter[67];
                  par2=mean/_parameter[67];
                  probaInd+=log(dgamma(_IndividualImpactTempTiterByInfec[subject][i][ser],par1,par2));
                  
                  mean=_parameter[68];//*_parameter[2];
                  par1=mean*mean/_parameter[69];
                  par2=mean/_parameter[69];
                  probaInd+=log(dgamma(_IndividualImpactSlopeByInfec[subject][i][ser],par1,par2));
              }
          }
      }
  }
    double bb;
  if(_vaccineSubject[subject]==1&_individualLevelEffectsAllByInfecVac==1&_baselineNaive[subject]==0&indicatorAll==0&indicatorVaccineAll==0){
    for (i=0;i<_numberVaccines[subject];i++){
      for (ser=0; ser<_nSerAssay; ser++){
        double m=_parameter[34];
        double par1=m*m/_parameter[35];
        double par2=m/_parameter[35];
        probaInd+=log(dgamma(_IndividualImpactTiterByInfecVac[subject][i][ser],par1,par2));
        if(probaInd==infinity|probaInd==negative_infinity){
          
        }
        m=_parameter[36];
        par1=m*m/_parameter[14];
        par2=m/_parameter[14];
        probaInd+=log(dgamma(_IndividualImpactTempTiterByInfecVac[subject][i][ser],par1,par2));
        if(probaInd==infinity|probaInd==negative_infinity){
          
        }
        m=_parameter[38];
        par1=m*m/_parameter[15];
        par2=m/_parameter[15];
        probaInd+=log(dgamma(_IndividualImpactSlopeByInfecVac[subject][i][ser],par1,par2));
        
        if(probaInd==infinity|probaInd==negative_infinity){
          
        }
      }
    }
  }
    if(_vaccineSubject[subject]==1&_individualLevelEffectsAllByInfecVac==1&_baselineNaive[subject]==1){
        for (i=0;i<_numberVaccines[subject];i++){
            for (ser=0; ser<_nSerAssay; ser++){
                double m=_parameter[44];
                double par1=m*m/_parameter[13];
                double par2=m/_parameter[13];
                probaInd+=log(dgamma(_IndividualImpactTiterByInfecVac[subject][i][ser],par1,par2));
                if(probaInd==infinity|probaInd==negative_infinity){
                    
                }
                
                m=_parameter[46];
                par1=m*m/_parameter[14];
                par2=m/_parameter[14];
                probaInd+=log(dgamma(_IndividualImpactTempTiterByInfecVac[subject][i][ser],par1,par2));
                if(probaInd==infinity|probaInd==negative_infinity){
                    
                }
                
                m=_parameter[48];
                par1=m*m/_parameter[15];
                par2=m/_parameter[15];
                probaInd+=log(dgamma(_IndividualImpactSlopeByInfecVac[subject][i][ser],par1,par2));
                if(probaInd==infinity|probaInd==negative_infinity){
                    
                }
            }
        }
    }
  _probaTiterIndividual[subject]=probaInd;

}


void computeProbaAll(double *globalLogLik,int SympOnly, int studyPeriodOnly, int vaccineBaselineOnly)
{   int subject;
  
  for(subject=0;subject<_numberOfSubject;subject++){
    computeProbaIndividual(subject, SympOnly, studyPeriodOnly,vaccineBaselineOnly);
//          cout<<subject<<" "<<_probaTiterIndividual[subject]<<"/n";
  }
  
  (*globalLogLik)=0;
  for(subject=0;subject<_numberOfSubject;subject++)
  {
    (*globalLogLik)+=_probaTiterIndividual[subject];
    double negative_infinity = - std::numeric_limits<double>::infinity();
    if((*globalLogLik)==negative_infinity){
      
    }
    double infinity = std::numeric_limits<double>::infinity();
    if((*globalLogLik)==infinity){
      
    }
  }
}













//================================
//========== build data ==========
void buildSubjectData(const char *dataFile)
{ifstream infile(dataFile);
  
  int subjectID,indid,nTests,nSymp,nVacc,nSympX,counter,i,ser,mDay;
  
  for(subjectID=0;subjectID<_numberOfSubject;subjectID++)
  {
    int minDayNaive,tmp;
    
    infile>>indid;
    infile>>_baselineNaive[indid];
    infile>>_vaccineSubject[indid];
    infile>>_numberOfTests[indid];
    infile>>_freqSymptomaticSerKnown[indid];
    infile>>_freqSymptomaticSerUnknown[indid];
    infile>>_minDayStudy[indid];
    infile>>_maxDayStudy[indid];
    infile>>_dayVacc[indid][0]; // HAVE FINAL DOSE AS THE DATE YOU USE!!
    
    if(_vaccineSubject[indid]==0){_minDayStudy[indid]=max(_minDayInfected,_minDayStudy[indid]-600);} //PD3 only for vaccinees

    if(_vaccineSubject[indid]==1){_numberVaccines[indid]=1;} //Just use first point for whole first year
    if(_vaccineSubject[indid]==0){_numberVaccines[indid]=0;} //Just use first point for whole first year
    if(_baselineNaive[indid]==2){_baselineNaive[indid]=1;} //Make indeterminate ones seropositive
    _baselineNaive[indid]=1-_baselineNaive[indid];  //Input is actually if seropositive

    nSymp=_freqSymptomaticSerKnown[indid];
    _daySymptomaticSerKnown[indid]=new int[_nSer];
    _serotypeSymptomatic[indid]=new int[_nSer];
    
    nSympX=_freqSymptomaticSerUnknown[subjectID];
    _daySymptomaticSerUnknown[subjectID]=new int[_nSer];
    _augmentedSerotypeSymptomatic[subjectID]=new int[_nSer];
    
      for(i=0;i<nSymp;i++){
          _daySymptomaticSerKnown[indid][i]=0;
          _serotypeSymptomatic[indid][i]=0;
      }
      for(i=0;i<nSympX;i++){
          _daySymptomaticSerUnknown[indid][i]=0;
          _augmentedSerotypeSymptomatic[indid][i]=0;
      }
    
  }
  
  for(subjectID=0;subjectID<_numberOfSubject;subjectID++)
  {   nTests=_numberOfTests[subjectID];
    _resultOfTesting[subjectID]=new double*[nTests];
    _meanTiterindividual[subjectID]=new double*[nTests];
    _dayOfTesting[subjectID]=new int[nTests];
    _minDelaySinceInfection[subjectID]=new double[nTests];
    _minTimPreInfection[subjectID]=new double[nTests];
    for(counter=0;counter<nTests;counter++)
    {   _resultOfTesting[subjectID][counter]=new double[4];
      _meanTiterindividual[subjectID][counter]=new double[_nSerAssay];
      _dayOfTesting[subjectID][counter]=0;
      _minDelaySinceInfection[subjectID][counter]=9999;
      _minTimPreInfection[subjectID][counter]=9999;
      for (i=0; i<4; i++)
      {   _resultOfTesting[subjectID][counter][i]=0;
      }
    for (i=0; i<_nSerAssay; i++)
    { _meanTiterindividual[subjectID][counter][i]=0;
    }
    }
  }
  
}

void buildTestingData(const char *dataFile)
{ifstream infile(dataFile);
  
  int indid,testno,testLine;
  for(testLine=0;testLine<_totalTests;testLine++)
  {
    infile>>indid;
    infile>>testno;
    infile>>_dayOfTesting[indid][testno-1];
    infile>>_resultOfTesting[indid][testno-1][0];
    infile>>_resultOfTesting[indid][testno-1][1];
    infile>>_resultOfTesting[indid][testno-1][2];
    infile>>_resultOfTesting[indid][testno-1][3];
      if(_calcMeanTiterImpact==1){
        _resultOfTesting[indid][testno-1][0]=(_resultOfTesting[indid][testno-1][0]+_resultOfTesting[indid][testno-1][1]+_resultOfTesting[indid][testno-1][2]+_resultOfTesting[indid][testno-1][3])/4;
          _resultOfTesting[indid][testno-1][1]=-999;
          _resultOfTesting[indid][testno-1][2]=-999;
          _resultOfTesting[indid][testno-1][3]=-999;
      }

  }
}

void buildSymptomaticSerKnownData(const char *dataFile)
{ifstream infile(dataFile);
  
  int indid,testno,testLine,ser,day;
  for(testLine=0;testLine<_numberNotMissingSymptomaticSerotype;testLine++)
  {
    infile>>indid;
    infile>>testno;
    infile>>day;
    infile>>ser;
    _daySymptomaticSerKnown[indid][testno-1]=day;
    _serotypeSymptomatic[indid][testno-1]=ser-1;
  }
}


void buildSymptomaticSerUnknownData(const char *dataFile)
{ifstream infile(dataFile);
  
  int indid,testno,testLine,i,j,k,l;
  for(testLine=0;testLine<_numberMissingSymptomaticSerotype;testLine++)
  {
    infile>>indid;
    infile>>testno;
    infile>>_daySymptomaticSerUnknown[indid][testno-1];
    infile>>i;
    _idsMissingSymptomaticSerotype[testLine]=indid;
    _ranksMissingSymptomaticSerotype[testLine]=testno-1;
  }
}

void buildPreviousAugmentedInfections(const char *dataFile)
{ifstream infile(dataFile);
  
  int indid,ser,testLine,date, minRise;
  double par1,par2;
  for(testLine=0;testLine<10000;testLine++){
    infile>>indid;
    if(indid==-999) break;
    infile>>date;
    infile>>ser;
      minRise=1;
    while(minRise==1){
        for (ser=0; ser<_nSerAssay;ser++){
            par1=_parameter[11]*_parameter[11]/_parameter[13];
            par2=_parameter[11]/_parameter[13];
            _tmpTvec[ser]=max(1e-10,rgamma(par1,par2));
        
            par1=_parameter[0]*_parameter[0]/_parameter[14];
            par2=_parameter[0]/_parameter[14];
            _tmpRvec[ser]=max(1e-10,rgamma(par1,par2));
            
            par1=_parameter[2]*_parameter[2]/_parameter[15];
            par2=_parameter[2]/_parameter[15];
            _tmpSvec[ser]=max(_minSlope,rgamma(par1,par2));
        }
        minRise=checkMinRise(indid);
    }
    addAsymptomaticInfection(indid,date,ser);
  }
}

void buildPreviousIndPars(const char *dataFile)
{ifstream infile(dataFile);
  
  int indid,testLine;
  double nan;
  for(testLine=0;testLine<10000;testLine++){
    infile>>indid;
    if(indid==-999) break;
    infile>>_IndividualImpactTiter[indid];
    infile>>_IndividualImpactDelay[indid];
    infile>>_IndividualImpactTempTiter[indid];
    infile>>_IndividualImpactSlope[indid];
    infile>>_baselineTiter[indid][0];
    infile>>_baselineTiter[indid][1];
    infile>>_baselineTiter[indid][2];
    infile>>_baselineTiter[indid][3];
    infile>>_baselineNaive[indid];
    infile>>nan;

  }
}





////// Make observed titers from simulated (for simulated)
void removeAllInfectionsForSim()
{
  _augmentedSerotypeSymptomatic=new int*[_numberOfSubject];
  _idsMissingSymptomaticSerotype=new int[_numberMissingSymptomaticSerotype];
  _ranksMissingSymptomaticSerotype=new int[_numberMissingSymptomaticSerotype];
  
  int subject, noInfec, i;
  for (subject=0;subject<_numberOfSubject;subject++)
  {
    _freqSymptomaticSerKnown[subject]=0;
    _freqSymptomaticSerUnknown[subject]=0;
    noInfec=_numberAsymptomatic[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++)
      {   removeAsymptomaticInfection(subject, i);}
    }
  }
  
}

void addAsymptomaticInfectionForSim(double probaInfOb, double probaSerOb)
{
  int subject,ser, day,i, daySinceInfec,a;
  double totTiter, u,aa;
  
  for(subject=0;subject<_numberOfSubject;subject++){
    
    int minDay=_minDayStudy[subject]+1;
    int maxDay=_maxDayStudy[subject]-1;
    daySinceInfec=9999;
    for (day=minDay;day<maxDay;day++){
      daySinceInfec++;
      if(day<(_minDayInfected+int(_incubationPeriod)+int(_delayRise))) continue;
      for (ser=0;ser<_nSer;ser++){
        int noInfec=_numberAsymptomatic[subject];
        int symp1=_freqSymptomaticSerKnown[subject];
        int symp2=_freqSymptomaticSerUnknown[subject];
        int indicator=0;
        for (i=0; i<noInfec;i++){if(_serotypeAsymptomatic[subject][i]==ser){indicator=1;}}
        for (i=0; i<symp1;i++){if(_serotypeSymptomatic[subject][i]==ser){indicator=1;}}
        for (i=0; i<symp2;i++){if(_augmentedSerotypeSymptomatic[subject][i]==ser){indicator=1;}}
        if (indicator==1) continue;
        if (daySinceInfec<maxTimeBetweenInfections) continue;
        if (_vaccineSubject[subject]==1&day<_dayVacc[subject][0]) continue;
        
        //                titerImpact=1;
        
        double probInfec=1-exp(-_parameter[12]*_epiPDF[day-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]);
        u=runif();
        if(u>probInfec){
          daySinceInfec++;
          continue;
        }
        
        daySinceInfec=0;
          double u=runif();
          double v=runif();
        
          double par1R,par2R,par1S,par2S,par1T,par2T;
        int sero, minRise;
          
          if(_vaccineSubject[subject]==0){
              par1T=_parameter[11]*_parameter[11]/_parameter[13];
              par2T=_parameter[11]/_parameter[13];
              par1R=_parameter[0]*_parameter[0]/_parameter[14];
            par2R=_parameter[0]/_parameter[14];
            par1S=_parameter[2]*_parameter[2]/_parameter[15];
            par2S=_parameter[2]/_parameter[15];
            
          }
          if(_vaccineSubject[subject]==1&_baselineNaive[subject]==1){
              double mean=_parameter[54];//*_parameter[11];
              par1T=mean*mean/_parameter[55];
              par2T=mean/_parameter[55];
              
              mean=_parameter[56];//*_parameter[0];
              par1R=mean*mean/_parameter[57];
              par2R=mean/_parameter[57];
              
              mean=_parameter[58];//*_parameter[2];
              par1S=mean*mean/_parameter[59];
              par2S=mean/_parameter[59];
          }
          if(_vaccineSubject[subject]==1&_baselineNaive[subject]==0){
              double mean=_parameter[64];//*_parameter[11];
              par1T=mean*mean/_parameter[65];
              par2T=mean/_parameter[65];
              
              mean=_parameter[66];//*_parameter[0];
              par1R=mean*mean/_parameter[67];
              par2R=mean/_parameter[67];
              
              mean=_parameter[68];//*_parameter[2];
              par1S=mean*mean/_parameter[69];
              par2S=mean/_parameter[69];
          }
        
        int seri;
        double bb;
        for (seri=0;seri<_nSerAssay;seri++){
            bb=max(1e-10,rgamma(par1S,par2S));
          _tmpSvec[seri]=max(1e-10,rgamma(par1S,par2S));
          _tmpTvec[seri]=max(1e-10,rgamma(par1T,par2T));
          _tmpRvec[seri]=max(1e-10,rgamma(par1R,par2R));
        }
        
        if(u<probaInfOb){
          if(v<probaSerOb){
            addSymptomaticInfectionSerKnown(subject, day, ser);
              minRise=checkMinRise(subject,0);
              while(minRise==1){
                  _IndividualImpactTiterByInfec[subject][0][0]=max(1e-10,rgamma(par1T,par2T));
                  _IndividualImpactTempTiterByInfec[subject][0][0]=max(1e-10,rgamma(par1R,par2R));
                  _IndividualImpactSlopeByInfec[subject][0][0]=max(1e-10,rgamma(par1S,par2S));
                  minRise=checkMinRise(subject,0);
              }
          }else{
            addSymptomaticInfectionSerUnknown(subject,day,ser);
          }
        }else{
            int noAsympInfec=_numberAsymptomatic[subject];
            int noSympInfec1=_freqSymptomaticSerKnown[subject];
            int noSympInfec2=_freqSymptomaticSerUnknown[subject];
            addAsymptomaticInfection(subject, day, ser);
            minRise=checkMinRise(subject,0);
            if(minRise==1){
                removeAsymptomaticInfection(subject, 0);
            }
        }
          minRise=checkMinRise(subject,0);
          if(minRise==1){
              
              
          }
      }
    }
    //Add testing then too
    int totinfs=_freqSymptomaticSerKnown[subject];
    int noTests=_numberOfTests[subject];
    int i;
    int counter=0;
            for(i=0;i<totinfs;i++){
                if (counter==noTests-1) continue;
                if (_daySymptomaticSerKnown[subject][i]<_minDayStudy[subject]) continue;
                if(_daySymptomaticSerKnown[subject][i]-14>minDay){_dayOfTesting[subject][counter]=_daySymptomaticSerKnown[subject][i]-14;
                    counter++;}
                if (counter==noTests-1) continue;
                if(_daySymptomaticSerKnown[subject][i]+3<maxDay){
                    _dayOfTesting[subject][counter]=_daySymptomaticSerKnown[subject][i]+3;
                    counter++;}
                if (counter==noTests-1) continue;
                if(_daySymptomaticSerKnown[subject][i]+7<maxDay){
                    _dayOfTesting[subject][counter]=_daySymptomaticSerKnown[subject][i]+7;
                    counter++;}
            }
  }
}

void addSerotypesAll()
{
  
  int id, subject, i,newSerotype=999;
  for (id=0; id<_numberMissingSymptomaticSerotype; id++)
  {
    subject=_idsMissingSymptomaticSerotype[id];
    int infecNo=_ranksMissingSymptomaticSerotype[id];
    int date=_daySymptomaticSerUnknown[subject][infecNo];
    
    double probSer;
    int indic=0;
    int noInfec;
    while(indic==0){
      int indicator=0;
      drawSer(date,&newSerotype,&probSer);
      noInfec=_freqSymptomaticSerKnown[subject];
      if(noInfec>0){
        for (i=0; i<noInfec;i++)
          {   if(abs(_daySymptomaticSerKnown[subject][i]-date)<maxTimeBetweenInfections){indicator=1;}
          if(_serotypeSymptomatic[subject][i]==newSerotype){indicator=1;}
          }
      }
      noInfec=_freqSymptomaticSerUnknown[subject];
      if(noInfec>0){
        for (i=0; i<noInfec;i++)
        {   if(i==infecNo)continue;
        if(_augmentedSerotypeSymptomatic[subject][i]==newSerotype){indicator=1;}
        }
      }
      if(indicator==0){indic=1;}
    }
    _augmentedSerotypeSymptomatic[subject][infecNo]=newSerotype;
  }
}

void MakeSimulatedTiters()
{
  int subject, noTests, testNo, ser,e,f,g,h;
  double titer,serImpact,c,d;
  for (subject=0;subject<_numberOfSubject;subject++){
    noTests=_numberOfTests[subject];
    for (testNo=0;testNo<noTests;testNo++){
      for(ser=0;ser<_nSerAssay;ser++)
      {   serImpact=0;
        if(ser==1){serImpact=_parameter[16];}
        if(ser==2){serImpact=_parameter[17];}
        if(ser==3){serImpact=_parameter[18];}
        titer=_meanTiterindividual[subject][testNo][ser]+(serImpact)+rnorm()*_parameter[8];
        if(titer<0){titer=0;}
        _resultOfTesting[subject][testNo][ser]=titer;
      }
    }
  }
}












//=========== update parameter ============
double updateParameter(int parameterNumber)
{
  double oldValue=_parameter[parameterNumber];
  double newValue=oldValue*exp(_rateForRandomWalk[parameterNumber]*rnorm());

  if(parameterNumber>=16&parameterNumber<19){
    newValue=oldValue+_rateForRandomWalk[parameterNumber]*rnorm();
  }
  
  if(parameterNumber==20 & newValue<0) return 0; //Can't have negative rise in titers
  if(parameterNumber==20 & newValue>1) return 0; //Can't have greater risk for high titers
  
  if(parameterNumber==32 & newValue>2*3.14159265358979323846) return 0; //2pi
  if(parameterNumber==32 & newValue<0) return 0; //between 0-2pi
  if(parameterNumber==31 & newValue>1) return 0; //Creates negative seasonality
  
  double _logLikGlobalSymp=0;
  if(_sympOnly==1&(parameterNumber==0|parameterNumber==2|parameterNumber==11|parameterNumber==13|parameterNumber==14|parameterNumber==15|parameterNumber>53)){
    computeProbaAll(&_logLikGlobalSymp,1,0,0);
  }
    if(_sympOnlyEpiCurve==1&(parameterNumber==31|parameterNumber==32)){
      computeProbaAll(&_logLikGlobalSymp,1,0,0);
    }
  
  if(_vaccineParOnly==1&(parameterNumber>33&parameterNumber<50)){
    computeProbaAll(&_logLikGlobalSymp,0,0,1);
  }
  
  double logProposal=log(newValue)-log(oldValue);
  if((parameterNumber>=16&parameterNumber<19)){
    logProposal=0;
  }
  _parameter[parameterNumber]=newValue;
  
  if(_seperateVaccineParameters==0){
    if(parameterNumber==13){
        _parameter[55]=newValue;
        _parameter[65]=newValue;
    }
    if(parameterNumber==14){
        _parameter[57]=newValue;
        _parameter[67]=newValue;
        
    }
    if(parameterNumber==15){
        _parameter[59]=newValue;
        _parameter[69]=newValue;
    }
    if(parameterNumber==2){
        _parameter[58]=newValue;
        _parameter[68]=newValue;
    }
  }
  
  double logPrior=0;
  double _varLNormPrior=1;
  double _meanLNormPrior=0;
  logPrior=logdlnorm(newValue,_meanLNormPrior,_varLNormPrior)-logdlnorm(oldValue,_meanLNormPrior,_varLNormPrior);
  
  if((parameterNumber>=16&parameterNumber<19)){
    logPrior=log(dnorm(newValue,0,5))-log(dnorm(oldValue,0,5));
  }
      
  if(parameterNumber>20&parameterNumber<33){
    updateEpiPDF();
  }
    

  double Q;
  if(_sympOnly==1&(parameterNumber==0|parameterNumber==2|parameterNumber==11|parameterNumber==13|parameterNumber==14|parameterNumber==15|parameterNumber>53)){
    computeProbaAll(&_newLogLikGlobal,1,0,0);
    Q=_newLogLikGlobal-_logLikGlobalSymp+logProposal+logPrior;
  }else if(_sympOnlyEpiCurve==1&(parameterNumber==31|parameterNumber==32)){
    computeProbaAll(&_newLogLikGlobal,1,0,0);
    Q=_newLogLikGlobal-_logLikGlobalSymp+logProposal+logPrior;
  }else if(_vaccineParOnly==1&(parameterNumber>33&parameterNumber<50)){
      computeProbaAll(&_newLogLikGlobal,0,0,1);
      Q=_newLogLikGlobal-_logLikGlobalSymp+logProposal+logPrior;
  }else{
    computeProbaAll(&_newLogLikGlobal,0,0,0);
    Q=_newLogLikGlobal-_logLikGlobal+logProposal+logPrior;
  }
    
    int subject;
    double minRise;
    minRise=0;
    for (subject=0;subject<_numberOfSubject;subject++){
        minRise=checkMinRise(subject);
        if(minRise==1)break;
    }
    
    if(log(runif())<Q&minRise==0){
        if(_sympOnly==1&(parameterNumber==0|parameterNumber==2|parameterNumber==11|parameterNumber==13|parameterNumber==14|parameterNumber==15|parameterNumber>53)){
      computeProbaAll(&_newLogLikGlobal,0,0,0);
      _logLikGlobal=_newLogLikGlobal;
        }else if(_sympOnlyEpiCurve==1&(parameterNumber==31|parameterNumber==32)){
      computeProbaAll(&_newLogLikGlobal,0,0,0);
      _logLikGlobal=_newLogLikGlobal;
        }else if(_vaccineParOnly==1&(parameterNumber>33&parameterNumber<50)){
        computeProbaAll(&_newLogLikGlobal,0,0,0);
        _logLikGlobal=_newLogLikGlobal;
    }else{
      _logLikGlobal=_newLogLikGlobal;
    }
        
        int minRise2=0;
        for (subject=0;subject<_numberOfSubject;subject++){
            minRise2=checkMinRise(subject);
            if(minRise2==1)break;
        }
        if(minRise2==1){
            
        }
        
    return 1;}
  else{_parameter[parameterNumber]=oldValue;
    if(_seperateVaccineParameters==0){
      if(parameterNumber==13){
          _parameter[55]=oldValue;
          _parameter[65]=oldValue;
      }
      if(parameterNumber==14){
          _parameter[57]=oldValue;
          _parameter[67]=oldValue;
          
      }
      if(parameterNumber==15){
          _parameter[59]=oldValue;
          _parameter[69]=oldValue;
      }
        if(parameterNumber==2){
            _parameter[58]=oldValue;
            _parameter[68]=oldValue;
        }
    }
    if(parameterNumber>20&parameterNumber<33){
      updateEpiPDF();
    }
    computeProbaAll(&_logLikGlobal,0,0,0);
    return 0;}
}





//=========== change serotype - Symptomatic ============
double updateSymptomaticSerotype(int nbIter)
{
  int iter,subject, infecNo,z,noInfec,indicator, i, nai;
  double accepted=0,successfullyProposed=0;
  for(iter=0;iter<nbIter;iter++)
  {   if(_numberMissingSymptomaticSerotype==0) continue;
  z=int(runif()*_numberMissingSymptomaticSerotype);
  
  subject=_idsMissingSymptomaticSerotype[z];
  
  infecNo=_ranksMissingSymptomaticSerotype[z];
  int date=_daySymptomaticSerUnknown[subject][infecNo];
  int oldSerotype=_augmentedSerotypeSymptomatic[subject][infecNo];
  double oldProba=probaSer(date,oldSerotype);
  double probSer;
  int newSerotype;
  drawSer(date,&newSerotype,&probSer);
  if(oldSerotype==newSerotype) continue;
  double logProposal=log(oldProba)-log(probSer);
  
  indicator=0;
  noInfec=_freqSymptomaticSerKnown[subject];
  if(noInfec>0){
    for (i=0; i<noInfec;i++)
      {   if(_serotypeSymptomatic[subject][i]==newSerotype){indicator=1;}
      if(_daySymptomaticSerKnown[subject][i]<date){indicator=1;}
      }
  }
  noInfec=_freqSymptomaticSerUnknown[subject];
  if(noInfec>0){
    for (i=0; i<noInfec;i++)
    {   if(i==infecNo) continue;
    if(_augmentedSerotypeSymptomatic[subject][i]==newSerotype){indicator=1;}
    if(_daySymptomaticSerUnknown[subject][i]<date){indicator=1;}
    }
  }
  noInfec=_numberAsymptomatic[subject];
  if(noInfec>0){
    for (i=0; i<noInfec;i++)
      {    if(_serotypeAsymptomatic[subject][i]==newSerotype){indicator=1;}
      if(_dayAsymptomatic[subject][i]<date){indicator=1;}
      }
  }
  if(indicator==1) continue;
  
  double _oldProbaSubject=_probaTiterIndividual[subject];
  _augmentedSerotypeSymptomatic[subject][infecNo]=newSerotype;
  computeProbaIndividual(subject,0,0,0);
  double _newProbaSubject=_probaTiterIndividual[subject];
  double _newLogLikGlobal=_logLikGlobal-_oldProbaSubject+_newProbaSubject;
  
  if(log(runif())<_newLogLikGlobal-_logLikGlobal+logProposal)
  {   _logLikGlobal=_newLogLikGlobal;
    accepted++;
  }
  else{
    _augmentedSerotypeSymptomatic[subject][infecNo]=oldSerotype;
    computeProbaIndividual(subject,0,0,0);
  }
  successfullyProposed++;
  
  }
  return accepted/successfullyProposed;
}


//=========== change serotype - Asymptomatic ============
double updateAsymptomaticSerotype(int nbIter)
{
  int iter,subject, infecNo,noInfec,indicator, i;
  double accepted=0,successfullyProposed=0;
  for(iter=0;iter<nbIter;iter++)
  {   subject=int(runif()*_numberOfSubject);
    
    int symp=_freqSymptomaticSerKnown[subject];
    int symp2=_freqSymptomaticSerUnknown[subject];
    
    int currentNoInf=_numberAsymptomatic[subject];
    if (currentNoInf==0) continue;
    
    int infecNo=int(runif()*(currentNoInf));
    
    int oldSerotype=_serotypeAsymptomatic[subject][infecNo];
    int date=_dayAsymptomatic[subject][infecNo];
    double oldProba=probaSer(date,oldSerotype);
    double probSer;
    int newSerotype;
    drawSer(date,&newSerotype,&probSer);
    if(oldSerotype==newSerotype) continue;
    double logProposal=log(oldProba)-log(probSer);
    
    indicator=0;
    noInfec=_freqSymptomaticSerKnown[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++)
        {   if(_serotypeSymptomatic[subject][i]==newSerotype){indicator=1;}
        if(_daySymptomaticSerKnown[subject][i]<date){indicator=1;}
        }
    }
    noInfec=_freqSymptomaticSerUnknown[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++)
        {    if(_augmentedSerotypeSymptomatic[subject][i]==newSerotype){indicator=1;}
        if(_daySymptomaticSerUnknown[subject][i]<date){indicator=1;}
        }
    }
    noInfec=_numberAsymptomatic[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++)
      {   if(i==infecNo) continue;
      if(_serotypeAsymptomatic[subject][i]==newSerotype){indicator=1;}
      if(_dayAsymptomatic[subject][i]<date){indicator=1;}
      }
    }
    if(indicator==1) continue;
    
    double _oldProbaSubject=_probaTiterIndividual[subject];
    _serotypeAsymptomatic[subject][infecNo]=newSerotype;
    computeProbaIndividual(subject,0,0,0);
    double _newProbaSubject=_probaTiterIndividual[subject];
    double _newLogLikGlobal=_logLikGlobal-_oldProbaSubject+_newProbaSubject;
    
    if(log(runif())<_newLogLikGlobal-_logLikGlobal+logProposal)
    {   _logLikGlobal=_newLogLikGlobal;
      accepted++;
    }
    else{
      _serotypeAsymptomatic[subject][infecNo]=oldSerotype;
      computeProbaIndividual(subject,0,0,0);
    }
    successfullyProposed++;
  }
  return accepted/successfullyProposed;
}



////=========== Individual impact - permanent rise ============
double updateBaselineTiter(int nbIter, int specificSubjectInd=0, int specificSubject=-999, int checkMinRiseInd=1)
{
  int iter,subject, infecNo, ser, oldPrimary, newPrimary;
  double accepted=0,successfullyProposed=0,newBaseline,oldPar16,oldPar17,oldPar18;
  for(iter=0;iter<nbIter;iter++)
  {
      if(specificSubjectInd==1){subject=specificSubject;}else{
          subject=int(runif()*_numberOfSubject);}
    ser=int(runif()*_nSerAssay);
    
    if(_estNaiveBaseline==1){oldPrimary=_baselineNaive[subject];}
//    if(_estNaiveBaseline==0&_baselineNaive[subject]==1)continue;
    double oldBaseline=_baselineTiter[subject][ser];
    if(_estNaiveBaseline==1){newPrimary=drawPrimary();}
    newBaseline=drawBaselineTiter(_baselineNaive[subject]);
    double logProposal=0;
    
    double _oldProbaSubject=_probaTiterIndividual[subject];
    _baselineTiter[subject][ser]=newBaseline;
      
    computeProbaIndividual(subject,0,0,0);
      
      int minRise=0;
      if(checkMinRiseInd==1){minRise=checkMinRise(subject);}
      
    double _newProbaSubject=_probaTiterIndividual[subject];
    double _newLogLikGlobal=_logLikGlobal-_oldProbaSubject+_newProbaSubject;
    if(log(runif())<_newLogLikGlobal-_logLikGlobal+logProposal&minRise==0)
    {   _logLikGlobal=_newLogLikGlobal;
      accepted++;
    }
    else{
      _baselineTiter[subject][ser]=oldBaseline;
      if(_estNaiveBaseline==1){_baselineNaive[subject]=oldPrimary;}
      computeProbaIndividual(subject,0,0,0);
    }
    successfullyProposed++;
    
  }
  return accepted/successfullyProposed;
}

/// ============ Update slow decay of titers =====================
double updateSlowDecayTiters(int nbIter)
{
    int iter,subject, ser;
    double accepted=0,successfullyProposed=0,newDecay;
    for(iter=0;iter<nbIter;iter++)
    {   subject=int(runif()*_numberOfSubject);
        ser=int(runif()*_nSerAssay);
        
//        if(_baselineNaive[subject]==1)continue;
        double oldDecay=_IndividualImpactSlowDecay[subject][ser];

        double par1=_parameter[3]*_parameter[3]/_parameter[19];
        double par2=_parameter[3]/_parameter[19];
        newDecay=max(rgamma(par1,par2),1e-20);
        double logProposal=0;
        
        double _oldProbaSubject=_probaTiterIndividual[subject];
        _IndividualImpactSlowDecay[subject][ser]=newDecay;
        
        int minRise=0;
        minRise=checkMinRise(subject);
        
        computeProbaIndividual(subject,0,0,0);
        double _newProbaSubject=_probaTiterIndividual[subject];
        double _newLogLikGlobal=_logLikGlobal-_oldProbaSubject+_newProbaSubject;
        if(log(runif())<_newLogLikGlobal-_logLikGlobal+logProposal&minRise==0)
        {   _logLikGlobal=_newLogLikGlobal;
            accepted++;
        }
        else{
            _IndividualImpactSlowDecay[subject][ser]=oldDecay;
            computeProbaIndividual(subject,0,0,0);
        }
        successfullyProposed++;
    }
    return accepted/successfullyProposed;
}


////=========== Individual impact - all infecs have diff ind impacts ============
double updateIndividualImpactAllInfAllPar(int nbIter, int specificSubjectInd=0, int specificSubject=-999)
{
  int iter,subject, infecNo, Ninfec,ser;
  double accepted=0,successfullyProposed=0;
  double tmpRise,aa,bb, tmpLL,tmpLL2,tmpLL3;
  for(iter=0;iter<nbIter;iter++)
  {    if(specificSubjectInd==1){subject=specificSubject;}else{
      subject=int(runif()*_numberOfSubject);}
    ser=int(runif()*_nSerAssay);
    
    Ninfec=_freqSymptomaticSerKnown[subject]+_freqSymptomaticSerUnknown[subject]+_numberAsymptomatic[subject];
    if(Ninfec==0) continue;
    infecNo=int(runif()*Ninfec);
    
      double oldBaseline=_baselineTiter[subject][ser];
      double newBaseline=drawBaselineTiter(_baselineNaive[subject]);
      
    double oldIndividualImpactR=_IndividualImpactTempTiterByInfec[subject][infecNo][ser];
    double oldIndividualImpactS=_IndividualImpactSlopeByInfec[subject][infecNo][ser];
    double oldIndividualImpactT=_IndividualImpactTiterByInfec[subject][infecNo][ser];
      
    double par1R,par2R,par1S,par2S,par1T,par2T,par1R_Vac,par2R_Vac,par1S_Vac,par2S_Vac,par1T_Vac,par2T_Vac, mean;
    if(_vaccineSubject[subject]==0){
      par1T=_parameter[11]*_parameter[11]/_parameter[13];
      par2T=_parameter[11]/_parameter[13];
      par1R=_parameter[0]*_parameter[0]/_parameter[14];
      par2R=_parameter[0]/_parameter[14];
      par1S=_parameter[2]*_parameter[2]/_parameter[15];
      par2S=_parameter[2]/_parameter[15];
    }
      if(_vaccineSubject[subject]==1&_baselineNaive[subject]==1){
          double mean=_parameter[54];//*_parameter[11];
          par1T=mean*mean/_parameter[55];
          par2T=mean/_parameter[55];
          
          mean=_parameter[56];//*_parameter[0];
          par1R=mean*mean/_parameter[57];
          par2R=mean/_parameter[57];
          
          mean=_parameter[58];//*_parameter[2];
          par1S=mean*mean/_parameter[59];
          par2S=mean/_parameter[59];
      }
      if(_vaccineSubject[subject]==1&_baselineNaive[subject]==0){
          double mean=_parameter[64];//*_parameter[11];
          par1T=mean*mean/_parameter[65];
          par2T=mean/_parameter[65];
          
          mean=_parameter[66];//*_parameter[0];
          par1R=mean*mean/_parameter[67];
          par2R=mean/_parameter[67];
          
          mean=_parameter[68];//*_parameter[2];
          par1S=mean*mean/_parameter[69];
          par2S=mean/_parameter[69];
      }
    
    double newIndividualImpactR=max(1e-10,rgamma(par1R,par2R));
    double newIndividualImpactS=max(_minSlope,rgamma(par1S,par2S));
    double newIndividualImpactT=max(1e-10,rgamma(par1T,par2T));
    
      double newIndividualImpactR_Vac, newIndividualImpactS_Vac, newIndividualImpactT_Vac;
      double oldIndividualImpactR_Vac=0;
      double oldIndividualImpactS_Vac=0;
      double oldIndividualImpactT_Vac=0;
      if(_vaccineSubject[subject]==1){
          oldIndividualImpactR_Vac=_IndividualImpactTempTiterByInfecVac[subject][0][ser];
          oldIndividualImpactS_Vac=_IndividualImpactSlopeByInfecVac[subject][0][ser];
          oldIndividualImpactT_Vac=_IndividualImpactTiterByInfecVac[subject][0][ser];
          
          if(_baselineNaive[subject]==0){
              mean=_parameter[34];
              par1T_Vac=mean*mean/_parameter[35];
              par2T_Vac=mean/_parameter[35];
              mean=_parameter[36];
              par1R_Vac=mean*mean/_parameter[14];
              par2R_Vac=mean/_parameter[14];
              mean=_parameter[38];
              par1S_Vac=mean*mean/_parameter[15];
              par2S_Vac=mean/_parameter[15];
          }
          if(_baselineNaive[subject]==1){
              mean=_parameter[44];
              par1T_Vac=mean*mean/_parameter[13];
              par2T_Vac=mean/_parameter[13];
              mean=_parameter[46];
              par1R_Vac=mean*mean/_parameter[14];
              par2R_Vac=mean/_parameter[14];
              mean=_parameter[48];
              par1S_Vac=mean*mean/_parameter[15];
              par2S_Vac=mean/_parameter[15];
          }
          
          newIndividualImpactR_Vac=max(1e-10,rgamma(par1R_Vac,par2R_Vac));
          newIndividualImpactS_Vac=max(1e-10,rgamma(par1S_Vac,par2S_Vac));
          newIndividualImpactT_Vac=max(1e-10,rgamma(par1T_Vac,par2T_Vac));
      }
      
    double _oldProbaSubject=_probaTiterIndividual[subject];
      
    _IndividualImpactTempTiterByInfec[subject][infecNo][ser]=newIndividualImpactR;
    _IndividualImpactSlopeByInfec[subject][infecNo][ser]=newIndividualImpactS;
    _IndividualImpactTiterByInfec[subject][infecNo][ser]=newIndividualImpactT;
      
      if(_vaccineSubject[subject]==1){
          _IndividualImpactTempTiterByInfecVac[subject][infecNo][ser]=newIndividualImpactR_Vac;
          _IndividualImpactSlopeByInfecVac[subject][infecNo][ser]=newIndividualImpactS_Vac;
          _IndividualImpactTiterByInfecVac[subject][infecNo][ser]=newIndividualImpactT_Vac;
      }
      
    _baselineTiter[subject][ser]=newBaseline;
      
    computeProbaIndividual(subject,0,0,0);

    int minRise=0;
    if(specificSubjectInd==0){minRise=checkMinRise(subject);}
      
    double logProposalR=log(dgamma(oldIndividualImpactR, par1R,par2R))-log(dgamma(newIndividualImpactR, par1R,par2R));
    double logProposalS=log(dgamma(oldIndividualImpactS, par1S,par2S))-log(dgamma(newIndividualImpactS, par1S,par2S));
    double logProposalT=log(dgamma(oldIndividualImpactT, par1T,par2T))-log(dgamma(newIndividualImpactT, par1T,par2T));
      double logProposalVac=0;
      if(_vaccineSubject[subject]==1){
          double logProposalR_Vac=log(dgamma(oldIndividualImpactR_Vac, par1R_Vac,par2R_Vac))-log(dgamma(newIndividualImpactR_Vac, par1R_Vac,par2R_Vac));
          double logProposalS_Vac=log(dgamma(oldIndividualImpactS_Vac, par1S_Vac,par2S_Vac))-log(dgamma(newIndividualImpactS_Vac, par1S_Vac,par2S_Vac));
          double logProposalT_Vac=log(dgamma(oldIndividualImpactT_Vac, par1T_Vac,par2T_Vac))-log(dgamma(newIndividualImpactT_Vac, par1T_Vac,par2T_Vac));
          logProposalVac=logProposalR_Vac+logProposalS_Vac+logProposalT_Vac;
      }
    double logProposal=logProposalR+logProposalS+logProposalT+logProposalVac;
    
    double _newProbaSubject=_probaTiterIndividual[subject];
    double _newLogLikGlobal=_logLikGlobal-_oldProbaSubject+_newProbaSubject;
    
    if(log(runif())<_newLogLikGlobal-_logLikGlobal+logProposal&minRise==0)
    {   _logLikGlobal=_newLogLikGlobal;
      accepted++;
    }
    else{
      _IndividualImpactTempTiterByInfec[subject][infecNo][ser]=oldIndividualImpactR;
      _IndividualImpactSlopeByInfec[subject][infecNo][ser]=oldIndividualImpactS;
      _IndividualImpactTiterByInfec[subject][infecNo][ser]=oldIndividualImpactT;
        if(_vaccineSubject[subject]==1){
            _IndividualImpactTempTiterByInfecVac[subject][infecNo][ser]=oldIndividualImpactR_Vac;
            _IndividualImpactSlopeByInfecVac[subject][infecNo][ser]=oldIndividualImpactS_Vac;
            _IndividualImpactTiterByInfecVac[subject][infecNo][ser]=oldIndividualImpactT_Vac;
        }
        _baselineTiter[subject][ser]=oldBaseline;
        computeProbaIndividual(subject,0,0,0);
    }
    successfullyProposed++;
    
  }
  return accepted/successfullyProposed;
}


////=========== Individual impact - all vaccinations have diff ind impacts ============
double updateIndividualImpactAllInfAllParVac(int nbIter, int specificSubjectInd=0, int specificSubject=-999)
{
  int iter,subject, vaccNo, NVacc,ser;
  double accepted=0,successfullyProposed=0;
  double meanR, par1,par2;
  for(iter=0;iter<nbIter;iter++)
  {
      if(specificSubjectInd==1){subject=specificSubject;}else{
          subject=int(runif()*_numberOfSubject);}
    ser=int(runif()*_nSerAssay);
    
    if(_vaccineSubject[subject]==0) continue;
    
    NVacc=_numberVaccines[subject];
    if(NVacc==0) continue;
    vaccNo=int(runif()*NVacc);
      
    double oldIndividualImpactR=_IndividualImpactTempTiterByInfecVac[subject][vaccNo][ser];
    double oldIndividualImpactS=_IndividualImpactSlopeByInfecVac[subject][vaccNo][ser];
    double oldIndividualImpactT=_IndividualImpactTiterByInfecVac[subject][vaccNo][ser];
    
    double oldBaseline=_baselineTiter[subject][ser];
    double newBaseline=drawBaselineTiter(_baselineNaive[subject]);
    
    double mean, par1T, par2T, par1R, par2R, par1S, par2S;
      if(_baselineNaive[subject]==0){
        mean=_parameter[34];
        par1T=mean*mean/_parameter[35];
        par2T=mean/_parameter[35];
        mean=_parameter[36];
        par1R=mean*mean/_parameter[14];
        par2R=mean/_parameter[14];
        mean=_parameter[38];
        par1S=mean*mean/_parameter[15];
        par2S=mean/_parameter[15];
      }
      if(_baselineNaive[subject]==1){
          mean=_parameter[44];
          par1T=mean*mean/_parameter[13];
          par2T=mean/_parameter[13];
          mean=_parameter[46];
          par1R=mean*mean/_parameter[14];
          par2R=mean/_parameter[14];
          mean=_parameter[48];
          par1S=mean*mean/_parameter[15];
          par2S=mean/_parameter[15];
      }
    
    double newIndividualImpactR=max(1e-10,rgamma(par1R,par2R));
    double newIndividualImpactS=max(1e-10,rgamma(par1S,par2S));
    double newIndividualImpactT=max(1e-10,rgamma(par1T,par2T));
    
    double titpre, titpost;
    titpre=computeTiterIndividualAtDay(subject,_dayVacc[subject][vaccNo]-1,ser);
    titpost=computeTiterIndividualAtDay(subject,_dayVacc[subject][vaccNo]+365,ser);
    if(titpost-titpre<_minTempTiterRise){continue;}
      
    double tmpRise;
    tmpRise=newIndividualImpactR*exp(-newIndividualImpactS*365);
    if(_linearDecay){tmpRise=newIndividualImpactR+(-newIndividualImpactS*365);}
    if(tmpRise<_minTempTiterRise){continue;}
    
    _baselineTiter[subject][ser]=newBaseline;
      
    double _oldProbaSubject=_probaTiterIndividual[subject];
    _IndividualImpactTempTiterByInfecVac[subject][vaccNo][ser]=newIndividualImpactR;
    _IndividualImpactSlopeByInfecVac[subject][vaccNo][ser]=newIndividualImpactS;
    _IndividualImpactTiterByInfecVac[subject][vaccNo][ser]=newIndividualImpactT;
      
    computeProbaIndividual(subject,0,0,0);

    int minRise=0;
    if(specificSubjectInd==0){minRise=checkMinRise(subject);}
      
    double logProposal=0;
    double logProposalR=log(dgamma(oldIndividualImpactR, par1R,par2R))-log(dgamma(newIndividualImpactR, par1R,par2R));
    double logProposalS=log(dgamma(oldIndividualImpactS, par1S,par2S))-log(dgamma(newIndividualImpactS, par1S,par2S));
    double logProposalT=log(dgamma(oldIndividualImpactT, par1T,par2T))-log(dgamma(newIndividualImpactT, par1T,par2T));
    logProposal=logProposalR+logProposalS+logProposalT;
    
    double _newProbaSubject=_probaTiterIndividual[subject];
    double _newLogLikGlobal=_logLikGlobal-_oldProbaSubject+_newProbaSubject;
    
    if(log(runif())<_newLogLikGlobal-_logLikGlobal+logProposal&minRise==0)
    {   _logLikGlobal=_newLogLikGlobal;
      accepted++;
    }
    else{
        _IndividualImpactTempTiterByInfecVac[subject][vaccNo][ser]=oldIndividualImpactR;
        _IndividualImpactSlopeByInfecVac[subject][vaccNo][ser]=oldIndividualImpactS;
        _IndividualImpactTiterByInfecVac[subject][vaccNo][ser]=oldIndividualImpactT;
        _baselineTiter[subject][ser]=oldBaseline;
        computeProbaIndividual(subject,0,0,0);
    }
    successfullyProposed++;
    
  }
  return accepted/successfullyProposed;
}




////=========== Individual impact - delay ============
double updateIndividualImpactDelay(int nbIter, int specificSubjectInd=0, int specificSubject=-999)
{
  int iter,subject;
  double accepted=0,successfullyProposed=0;
  for(iter=0;iter<nbIter;iter++)
  { if(specificSubjectInd==1){subject=specificSubject;}else{subject=int(runif()*_numberOfSubject);}
    
    int symp=_freqSymptomaticSerKnown[subject];
    int symp2=_freqSymptomaticSerUnknown[subject];
    if(symp==0&symp2==0) continue;
    
    double oldIndividualImpact=_IndividualImpactDelay[subject];
    double newIndividualImpact=runif()*(2*_maxDelayRise)-_maxDelayRise;
    double logProposal=0;
    
    int indicator=0;
    int i;
    int noInfec=_freqSymptomaticSerKnown[subject];
    for (i=0; i<noInfec;i++){   if(_daySymptomaticSerKnown[subject][i]+newIndividualImpact-_incubationPeriod<_minDayInfected){indicator=1;}
    }
    noInfec=_freqSymptomaticSerUnknown[subject];
    for (i=0; i<noInfec;i++){if(_daySymptomaticSerUnknown[subject][i]+newIndividualImpact-_incubationPeriod<_minDayInfected){indicator=1;}
    }
    noInfec=_numberAsymptomatic[subject];
    for (i=0; i<noInfec;i++){if(_dayAsymptomatic[subject][i]+newIndividualImpact-_incubationPeriod-_delayRise<_minDayInfected){indicator=1;}
    }
    if(indicator==1) continue;
    
    double _oldProbaSubject=_probaTiterIndividual[subject];
    computeProbaIndividual(subject,0,0,0);
    _IndividualImpactDelay[subject]=newIndividualImpact;
    computeProbaIndividual(subject,0,0,0);
    double _newProbaSubject=_probaTiterIndividual[subject];
    _newLogLikGlobal=_logLikGlobal-_oldProbaSubject+_newProbaSubject;
    
      int minRise=0;
      if(specificSubjectInd==0){minRise=checkMinRise(subject);}
      
    if(log(runif())<_newLogLikGlobal-_logLikGlobal+logProposal&minRise==0)
    {   _logLikGlobal=_newLogLikGlobal;
      accepted++;
    }
    else{
      _IndividualImpactDelay[subject]=oldIndividualImpact;
      computeProbaIndividual(subject,0,0,0);
    }
    
    successfullyProposed++;
  }
  return accepted/successfullyProposed;
}


//=========== add/remove asymptomatic infections ============
void RJMCMCAddAsymptomaticInfection(double *accepted,double *proposed, int specificSubjectInd=0, int specificSubject=-999)
{   int iter,subject,noInfec,indicator,i,ser;
  int nbIter=1;
  double par1,par2,tmpR;
          double aaa,tmpLL;
  for(iter=0;iter<nbIter;iter++)
  {   if(specificSubjectInd==1){subject=specificSubject;}else{
      subject=int(runif()*_numberOfSubject);
  }
    
    double oldtestLL;
    
    int symp=_freqSymptomaticSerKnown[subject];
    int symp2=_freqSymptomaticSerUnknown[subject];
    
    int currentNoInf=_numberAsymptomatic[subject];
    int currentTotInf=currentNoInf+symp+symp2;
    
    if(currentTotInf==_nSer)continue;
      
    int minDay=_minDayStudy[subject]+1;
    int maxDay=_maxDayStudy[subject]-1;
    if(minDay>maxDay) continue;
    
    double probSer;
    int newSerotype;
    drawSerOverall(minDay, maxDay,&newSerotype,&probSer);
    probSer=1;
    
    int newDate;
    double probDay;
    drawDateAsymptomatic(minDay, maxDay, &newDate, &probDay);
    
    if(newDate+_parameter[1]-_incubationPeriod-_delayRise<_minDayInfected) continue;

    indicator=0;
    noInfec=_freqSymptomaticSerKnown[subject];
    
    if(noInfec>0){
      for (i=0; i<noInfec;i++)
        {   if(abs(_daySymptomaticSerKnown[subject][i]-newDate)<maxTimeBetweenInfections){indicator=1;}
        if(_serotypeSymptomatic[subject][i]==newSerotype){indicator=1;}
        }
    }
    noInfec=_freqSymptomaticSerUnknown[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++)
        {   if(abs(_daySymptomaticSerUnknown[subject][i]-newDate)<maxTimeBetweenInfections){indicator=1;}
        if(_augmentedSerotypeSymptomatic[subject][i]==newSerotype){indicator=1;}
        }
    }
    noInfec=_numberAsymptomatic[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++)
        {   if(abs(_dayAsymptomatic[subject][i]-newDate)<maxTimeBetweenInfections){indicator=1;}
        if(_serotypeAsymptomatic[subject][i]==newSerotype){indicator=1;}
        }
    }
    if(_vaccineSubject[subject]==1){
      if(abs(_dayVacc[subject][0]-newDate)<_minTimeInfecPostVaccine){indicator=1;}
    }
    if(indicator==1) continue;
      
    int inf;
    double meanR,probDecayOld,probDecayNew;
    int ser;
    
      if(_individualLevelEffectsAllByInfec==1){
        for (inf=0;inf<currentTotInf;inf++){
          for (ser=0; ser<_nSerAssay;ser++){
            _oldTiters[inf][ser]=_IndividualImpactTiterByInfec[subject][inf][ser];
            _oldTempTiters[inf][ser]=_IndividualImpactTempTiterByInfec[subject][inf][ser];
            _oldSlopes[inf][ser]=_IndividualImpactSlopeByInfec[subject][inf][ser];
          }
        }
      }
      
      if(_individualLevelEffectsAllByInfecVac==1&_vaccineSubject[subject]==1){
         for (ser=0; ser<_nSerAssay;ser++){
          _oldTitersVac[0][ser]=_IndividualImpactTiterByInfecVac[subject][0][ser];
          _oldTempTitersVac[0][ser]=_IndividualImpactTempTiterByInfecVac[subject][0][ser];
          _oldSlopesVac[0][ser]=_IndividualImpactSlopeByInfecVac[subject][0][ser];
        }
    }
    
    //decay and baseline
      probDecayOld=0;
      probDecayNew=0;
      par1=_parameter[3]*_parameter[3]/_parameter[19];
      par2=_parameter[3]/_parameter[19];
//      if(_estNaiveBaseline==0&_baselineNaive[subject]==0){
        for (ser=0; ser<_nSerAssay;ser++){
          _oldBaselines[ser]=_baselineTiter[subject][ser];
          _baselineTiter[subject][ser]=drawBaselineTiter(_baselineNaive[subject]);
            if(_individualLevelEffectsSlowDecay==1){
              _oldSlowDecay[ser]=_IndividualImpactSlowDecay[subject][ser];
              probDecayOld+=log(dgamma(_IndividualImpactSlowDecay[subject][ser],par1,par2));
                _IndividualImpactSlowDecay[subject][ser]=max(rgamma(par1,par2),1e-10);
              probDecayNew+=log(dgamma(_IndividualImpactSlowDecay[subject][ser],par1,par2));
            }
        }
//    }
    double _oldProbaSubject=_probaTiterIndividual[subject];

    double tmpT, tmpS;
      double tmpT2, tmpS2, tmpR2;
    double probCurrentTiter=0;
      if(_individualLevelEffectsAllByInfec==1&_vaccineSubject[subject]==0){
        for (inf=0;inf<currentTotInf;inf++){
          for (ser=0;ser<_nSerAssay;ser++){
            par1=_parameter[11]*_parameter[11]/_parameter[13];
            par2=_parameter[11]/_parameter[13];
            tmpT=_IndividualImpactTiterByInfec[subject][inf][ser];
            probCurrentTiter+=log(dgamma(tmpT,par1,par2));

            meanR=_parameter[0];
            par1=meanR*meanR/_parameter[14];
            par2=meanR/_parameter[14];
            tmpR=_IndividualImpactTempTiterByInfec[subject][inf][ser];
            probCurrentTiter+=log(dgamma(tmpR,par1,par2));

            par1=_parameter[2]*_parameter[2]/_parameter[15];
            par2=_parameter[2]/_parameter[15];
            tmpS=_IndividualImpactSlopeByInfec[subject][inf][ser];
            probCurrentTiter+=log(dgamma(tmpS,par1,par2));
          }
        }
      }
      if(_individualLevelEffectsAllByInfec==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==1){
          for (inf=0;inf<currentTotInf;inf++){
              for (ser=0;ser<_nSerAssay;ser++){
                  double mean=_parameter[54];//*_parameter[11];
                  par1=mean*mean/_parameter[55];
                  par2=mean/_parameter[55];
                  tmpT=_IndividualImpactTiterByInfec[subject][inf][ser];
                  probCurrentTiter+=log(dgamma(tmpT,par1,par2));
                  
                  meanR=_parameter[56];//*_parameter[0];
                  par1=meanR*meanR/_parameter[57];
                  par2=meanR/_parameter[57];
                  tmpR=_IndividualImpactTempTiterByInfec[subject][inf][ser];
                  probCurrentTiter+=log(dgamma(tmpR,par1,par2));
                  
                  mean=_parameter[58];//*_parameter[2];
                  par1=mean*mean/_parameter[59];
                  par2=mean/_parameter[59];
                  tmpS=_IndividualImpactSlopeByInfec[subject][inf][ser];
                  probCurrentTiter+=log(dgamma(tmpS,par1,par2));
              }
          }
      }
      if(_individualLevelEffectsAllByInfec==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==0){
          for (inf=0;inf<currentTotInf;inf++){
              for (ser=0;ser<_nSerAssay;ser++){
                  double mean=_parameter[64];//*_parameter[11];
                  par1=mean*mean/_parameter[65];
                  par2=mean/_parameter[65];
                  tmpT=_IndividualImpactTiterByInfec[subject][inf][ser];
                  probCurrentTiter+=log(dgamma(tmpT,par1,par2));
                  
                  meanR=_parameter[66];//*_parameter[0];
                  par1=meanR*meanR/_parameter[67];
                  par2=meanR/_parameter[67];
                  tmpR=_IndividualImpactTempTiterByInfec[subject][inf][ser];
                  probCurrentTiter+=log(dgamma(tmpR,par1,par2));
                  
                  mean=_parameter[68];//*_parameter[2];
                  par1=mean*mean/_parameter[69];
                  par2=mean/_parameter[69];
                  tmpS=_IndividualImpactSlopeByInfec[subject][inf][ser];
                  probCurrentTiter+=log(dgamma(tmpS,par1,par2));
              }
          }
      }
      if(_individualLevelEffectsAllByInfecVac==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==0){
          for (ser=0;ser<_nSerAssay;ser++){
              par1=_parameter[34]*_parameter[34]/_parameter[35];
              par2=_parameter[34]/_parameter[35];
              tmpT=_IndividualImpactTiterByInfecVac[subject][0][ser];
              probCurrentTiter+=log(dgamma(tmpT,par1,par2));

              meanR=_parameter[36];
              par1=meanR*meanR/_parameter[14];
              par2=meanR/_parameter[14];
              tmpR=_IndividualImpactTempTiterByInfecVac[subject][0][ser];
              probCurrentTiter+=log(dgamma(tmpR,par1,par2));

              par1=_parameter[38]*_parameter[38]/_parameter[15];
              par2=_parameter[38]/_parameter[15];
              tmpS=_IndividualImpactSlopeByInfecVac[subject][0][ser];
              probCurrentTiter+=log(dgamma(tmpS,par1,par2));
          }
      }
      if(_individualLevelEffectsAllByInfecVac==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==1){
          for (ser=0;ser<_nSerAssay;ser++){
              par1=_parameter[44]*_parameter[44]/_parameter[13];
              par2=_parameter[44]/_parameter[13];
              tmpT=_IndividualImpactTiterByInfecVac[subject][0][ser];
              probCurrentTiter+=log(dgamma(tmpT,par1,par2));
              
              meanR=_parameter[46];
              par1=meanR*meanR/_parameter[14];
              par2=meanR/_parameter[14];
              tmpR=_IndividualImpactTempTiterByInfecVac[subject][0][ser];
              probCurrentTiter+=log(dgamma(tmpR,par1,par2));
              
              par1=_parameter[48]*_parameter[48]/_parameter[15];
              par2=_parameter[48]/_parameter[15];
              tmpS=_IndividualImpactSlopeByInfecVac[subject][0][ser];
              probCurrentTiter+=log(dgamma(tmpS,par1,par2));
          }
      }

    int indicator;
    double probFutureTiter=0;
    double a1,a2,a3;
    int bb;
    int minRise=0;
      
    addAsymptomaticInfection(subject, newDate, newSerotype);
      
     if(_individualLevelEffectsAllByInfec==1&_vaccineSubject[subject]==0){
        for (inf=0;inf<currentTotInf+1;inf++){
          for (ser=0; ser<_nSerAssay;ser++){
                par1=_parameter[11]*_parameter[11]/_parameter[13];
                par2=_parameter[11]/_parameter[13];
                tmpT2=max(1e-10,rgamma(par1,par2));
                _IndividualImpactTiterByInfec[subject][inf][ser]=tmpT2;
                probFutureTiter+=log(dgamma(tmpT2,par1,par2));

                meanR=_parameter[0];
                par1=meanR*meanR/_parameter[14];
                par2=meanR/_parameter[14];
                tmpR2=max(1e-10,rgamma(par1,par2));
                _IndividualImpactTempTiterByInfec[subject][inf][ser]=tmpR2;
                probFutureTiter+=log(dgamma(tmpR2,par1,par2));

                par1=_parameter[2]*_parameter[2]/_parameter[15];
                par2=_parameter[2]/_parameter[15];
                tmpS2=max(_minSlope,rgamma(par1,par2));
                _IndividualImpactSlopeByInfec[subject][inf][ser]=tmpS2;
                probFutureTiter+=log(dgamma(tmpS2,par1,par2));
          }
        }
     }
      if(_individualLevelEffectsAllByInfec==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==1){
          for (inf=0;inf<currentTotInf+1;inf++){
              for (ser=0; ser<_nSerAssay;ser++){
                  double mean=_parameter[54];//*_parameter[11];
                  par1=mean*mean/_parameter[55];
                  par2=mean/_parameter[55];
                  tmpT2=max(1e-10,rgamma(par1,par2));
                  _IndividualImpactTiterByInfec[subject][inf][ser]=tmpT2;
                  probFutureTiter+=log(dgamma(tmpT2,par1,par2));
                  
                  meanR=_parameter[56];//*_parameter[0];
                  par1=meanR*meanR/_parameter[57];
                  par2=meanR/_parameter[57];
                  tmpR2=max(1e-10,rgamma(par1,par2));
                  _IndividualImpactTempTiterByInfec[subject][inf][ser]=tmpR2;
                  probFutureTiter+=log(dgamma(tmpR2,par1,par2));
                  
                  mean=_parameter[58];//*_parameter[2];
                  par1=mean*mean/_parameter[59];
                  par2=mean/_parameter[59];
                  tmpS2=max(_minSlope,rgamma(par1,par2));
                  _IndividualImpactSlopeByInfec[subject][inf][ser]=tmpS2;
                  probFutureTiter+=log(dgamma(tmpS2,par1,par2));
              }
          }
      }
      if(_individualLevelEffectsAllByInfec==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==0){
          for (inf=0;inf<currentTotInf+1;inf++){
              for (ser=0; ser<_nSerAssay;ser++){
                  double mean=_parameter[64];//*_parameter[11];
                  par1=mean*mean/_parameter[65];
                  par2=mean/_parameter[65];
                  tmpT2=max(1e-10,rgamma(par1,par2));
                  _IndividualImpactTiterByInfec[subject][inf][ser]=tmpT2;
                  probFutureTiter+=log(dgamma(tmpT2,par1,par2));
                  
                  meanR=_parameter[66];//*_parameter[0];
                  par1=meanR*meanR/_parameter[67];
                  par2=meanR/_parameter[67];
                  tmpR2=max(1e-10,rgamma(par1,par2));
                  _IndividualImpactTempTiterByInfec[subject][inf][ser]=tmpR2;
                  probFutureTiter+=log(dgamma(tmpR2,par1,par2));
                  
                  mean=_parameter[68];//*_parameter[2];
                  par1=mean*mean/_parameter[69];
                  par2=mean/_parameter[69];
                  tmpS2=max(_minSlope,rgamma(par1,par2));
                  _IndividualImpactSlopeByInfec[subject][inf][ser]=tmpS2;
                  probFutureTiter+=log(dgamma(tmpS2,par1,par2));
              }
          }
      }
      double tmpRise;
      int minRise2=0;
    if(_individualLevelEffectsAllByInfecVac==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==0){
         for (ser=0; ser<_nSerAssay;ser++){
             par1=_parameter[34]*_parameter[34]/_parameter[35];
             par2=_parameter[34]/_parameter[35];
             tmpT=max(1e-10,rgamma(par1,par2));
             _IndividualImpactTiterByInfecVac[subject][0][ser]=tmpT;
             probFutureTiter+=log(dgamma(tmpT,par1,par2));

             par1=_parameter[36]*_parameter[36]/_parameter[14];
             par2=_parameter[36]/_parameter[14];
             tmpR=max(1e-10,rgamma(par1,par2));
             _IndividualImpactTempTiterByInfecVac[subject][0][ser]=tmpR;
             probFutureTiter+=log(dgamma(tmpR,par1,par2));

             par1=_parameter[38]*_parameter[38]/_parameter[15];
             par2=_parameter[38]/_parameter[15];
             tmpS=max(_minSlope,rgamma(par1,par2));
             _IndividualImpactSlopeByInfecVac[subject][0][ser]=tmpS;
             probFutureTiter+=log(dgamma(tmpS,par1,par2));
             
             tmpRise=tmpR*exp(-tmpS*365);
             if(tmpRise<_minTempTiterRise){minRise2=1;}
          }
      }
      if(_individualLevelEffectsAllByInfecVac==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==1){
          for (ser=0; ser<_nSerAssay;ser++){
              par1=_parameter[44]*_parameter[44]/_parameter[13];
              par2=_parameter[44]/_parameter[13];
              tmpT=max(1e-10,rgamma(par1,par2));
              _IndividualImpactTiterByInfecVac[subject][0][ser]=tmpT;
              probFutureTiter+=log(dgamma(tmpT,par1,par2));
              
              par1=_parameter[46]*_parameter[46]/_parameter[14];
              par2=_parameter[46]/_parameter[14];
              tmpR=max(1e-10,rgamma(par1,par2));
              _IndividualImpactTempTiterByInfecVac[subject][0][ser]=tmpR;
              probFutureTiter+=log(dgamma(tmpR,par1,par2));
              
              par1=_parameter[48]*_parameter[48]/_parameter[15];
              par2=_parameter[48]/_parameter[15];
              tmpS=max(_minSlope,rgamma(par1,par2));
              _IndividualImpactSlopeByInfecVac[subject][0][ser]=tmpS;
              probFutureTiter+=log(dgamma(tmpS,par1,par2));
              
              tmpRise=tmpR*exp(-tmpS*365);
              if(tmpRise<_minTempTiterRise){minRise2=1;}
          }
      }

    computeProbaIndividual(subject,0,0,0);
      
    if(specificSubjectInd==0){minRise2=checkMinRise(subject);}
      
    double _newProbaSubject=_probaTiterIndividual[subject];
    
    double logProbaAddCase=probDay+probFutureTiter+probDecayNew;
    double logProbaRemoveCase=log(1)+probCurrentTiter+probDecayOld;
    double logProposal=logProbaRemoveCase-logProbaAddCase;
    
    double _newLogLikGlobal;
    _newLogLikGlobal=_logLikGlobal-_oldProbaSubject+_newProbaSubject;
            
    if(log(runif())<_newLogLikGlobal-_logLikGlobal+logProposal&minRise==0&minRise2==0)
    {   _logLikGlobal=_newLogLikGlobal;
      (*accepted)++;
    }
    else{
      removeAsymptomaticInfection(subject, 0);
        if(_individualLevelEffectsAllByInfec==1){
          for (inf=0;inf<currentTotInf;inf++){
            for (ser=0;ser<_nSerAssay;ser++){
              _IndividualImpactTiterByInfec[subject][inf][ser]=_oldTiters[inf][ser];
              _IndividualImpactTempTiterByInfec[subject][inf][ser]=_oldTempTiters[inf][ser];
              _IndividualImpactSlopeByInfec[subject][inf][ser]=_oldSlopes[inf][ser];
            }
          }
        }
        if(_individualLevelEffectsAllByInfecVac==1&_vaccineSubject[subject]==1){
            for (ser=0;ser<_nSerAssay;ser++){
                _IndividualImpactTiterByInfecVac[subject][0][ser]=_oldTitersVac[0][ser];
                _IndividualImpactTempTiterByInfecVac[subject][0][ser]=_oldTempTitersVac[0][ser];
                _IndividualImpactSlopeByInfecVac[subject][0][ser]=_oldSlopesVac[0][ser];
            }
        }
      for (ser=0;ser<_nSerAssay;ser++){
        _baselineTiter[subject][ser]=_oldBaselines[ser];
          if(_individualLevelEffectsSlowDecay==1){_IndividualImpactSlowDecay[subject][ser]=_oldSlowDecay[ser];}
      }
      computeProbaIndividual(subject,0,0,0);
    }

    (*proposed)++;
  }
}




void RJMCMCRemoveAsymptomaticInfection(double *accepted,double *proposed, int specificSubjectInd=0,int specificSubject=-999)
{   int iter,subject, ser,i, inf;
  double par1,par2,probDecayOld,probDecayNew;
  int nbIter=1;
  for(iter=0;iter<nbIter;iter++){
      if(specificSubjectInd==1){subject=specificSubject;}else{
         subject=int(runif()*_numberOfSubject);
      }
      
    int symp=_freqSymptomaticSerKnown[subject];
    int symp2=_freqSymptomaticSerUnknown[subject];
    
    int currentNoInf=_numberAsymptomatic[subject];
    int currentTotInf=currentNoInf+symp+symp2;
    if (currentNoInf==0) continue;
    
    int infecNo=int(runif()*(currentNoInf));
    int Date=_dayAsymptomatic[subject][infecNo];
    int Serotype=_serotypeAsymptomatic[subject][infecNo];
      
    if(_individualLevelEffectsAllByInfec==1){
        for (ser=0;ser<_nSerAssay;ser++){
          _tmpTvec[ser]=_IndividualImpactTiterByInfec[subject][symp+symp2+infecNo][ser];
          _tmpRvec[ser]=_IndividualImpactTempTiterByInfec[subject][symp+symp2+infecNo][ser];
          _tmpSvec[ser]=_IndividualImpactSlopeByInfec[subject][symp+symp2+infecNo][ser];
        }
    }
      if(_individualLevelEffectsAllByInfecVac==1&_vaccineSubject[subject]==1){
          for (ser=0; ser<_nSerAssay;ser++){
              _oldTitersVac[0][ser]=_IndividualImpactTiterByInfecVac[subject][0][ser];
              _oldTempTitersVac[0][ser]=_IndividualImpactTempTiterByInfecVac[subject][0][ser];
              _oldSlopesVac[0][ser]=_IndividualImpactSlopeByInfecVac[subject][0][ser];
          }
      }
      
      //decay and baseline
      probDecayOld=0;
      probDecayNew=0;
      par1=_parameter[3]*_parameter[3]/_parameter[19];
      par2=_parameter[3]/_parameter[19];
//      if(_estNaiveBaseline==0&_baselineNaive[subject]==0){
          for (ser=0; ser<_nSerAssay;ser++){
              _oldBaselines[ser]=_baselineTiter[subject][ser];
              _baselineTiter[subject][ser]=drawBaselineTiter(_baselineNaive[subject]);
              if(_individualLevelEffectsSlowDecay==1){
                  _oldSlowDecay[ser]=_IndividualImpactSlowDecay[subject][ser];
                  probDecayOld+=log(dgamma(_IndividualImpactSlowDecay[subject][ser],par1,par2));
                  _IndividualImpactSlowDecay[subject][ser]=max(rgamma(par1,par2),1e-20);
                  probDecayNew+=log(dgamma(_IndividualImpactSlowDecay[subject][ser],par1,par2));
              }
          }
//      }
      
    double baseRise, meanR, tmpR, tmpS, tmpT;
    double probCurrentTiter=0;
    if(_individualLevelEffectsAllByInfec==1&_vaccineSubject[subject]==0){
      for (inf=0;inf<currentTotInf;inf++){
        for (ser=0;ser<_nSerAssay;ser++){
          par1=_parameter[11]*_parameter[11]/_parameter[13];
          par2=_parameter[11]/_parameter[13];
          probCurrentTiter+=log(dgamma(_IndividualImpactTiterByInfec[subject][inf][ser],par1,par2));
          
          meanR=_parameter[0];
          par1=meanR*meanR/_parameter[14];
          par2=meanR/_parameter[14];
          probCurrentTiter+=log(dgamma(_IndividualImpactTempTiterByInfec[subject][inf][ser],par1,par2));
          
          par1=_parameter[2]*_parameter[2]/_parameter[15];
          par2=_parameter[2]/_parameter[15];
          probCurrentTiter+=log(dgamma(_IndividualImpactSlopeByInfec[subject][inf][ser],par1,par2));
        }
      }
    }
      if(_individualLevelEffectsAllByInfec==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==1){
          for (inf=0;inf<currentTotInf;inf++){
              for (ser=0;ser<_nSerAssay;ser++){
                  double mean=_parameter[54];//*_parameter[11];
                  par1=mean*mean/_parameter[55];
                  par2=mean/_parameter[55];
                  probCurrentTiter+=log(dgamma(_IndividualImpactTiterByInfec[subject][inf][ser],par1,par2));
                  
                  meanR=_parameter[56];//*_parameter[0];
                  par1=meanR*meanR/_parameter[57];
                  par2=meanR/_parameter[57];
                  probCurrentTiter+=log(dgamma(_IndividualImpactTempTiterByInfec[subject][inf][ser],par1,par2));
                  
                  mean=_parameter[58];//*_parameter[2];
                  par1=mean*mean/_parameter[59];
                  par2=mean/_parameter[59];
                  probCurrentTiter+=log(dgamma(_IndividualImpactSlopeByInfec[subject][inf][ser],par1,par2));
              }
          }
      }
      if(_individualLevelEffectsAllByInfec==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==0){
          for (inf=0;inf<currentTotInf;inf++){
              for (ser=0;ser<_nSerAssay;ser++){
                  double mean=_parameter[64];//*_parameter[11];
                  par1=mean*mean/_parameter[65];
                  par2=mean/_parameter[65];
                  probCurrentTiter+=log(dgamma(_IndividualImpactTiterByInfec[subject][inf][ser],par1,par2));
                  
                  meanR=_parameter[66];//*_parameter[0];
                  par1=meanR*meanR/_parameter[67];
                  par2=meanR/_parameter[67];
                  probCurrentTiter+=log(dgamma(_IndividualImpactTempTiterByInfec[subject][inf][ser],par1,par2));
                  
                  mean=_parameter[68];//*_parameter[2];
                  par1=mean*mean/_parameter[69];
                  par2=mean/_parameter[69];
                  probCurrentTiter+=log(dgamma(_IndividualImpactSlopeByInfec[subject][inf][ser],par1,par2));
              }
          }
      }
      if(_individualLevelEffectsAllByInfecVac==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==0){
          for (ser=0;ser<_nSerAssay;ser++){
              par1=_parameter[34]*_parameter[34]/_parameter[35];
              par2=_parameter[34]/_parameter[35];
              tmpT=_IndividualImpactTiterByInfecVac[subject][0][ser];
              probCurrentTiter+=log(dgamma(tmpT,par1,par2));

              meanR=_parameter[36];
              par1=meanR*meanR/_parameter[14];
              par2=meanR/_parameter[14];
              tmpR=_IndividualImpactTempTiterByInfecVac[subject][0][ser];
              probCurrentTiter+=log(dgamma(tmpR,par1,par2));

              par1=_parameter[38]*_parameter[38]/_parameter[15];
              par2=_parameter[38]/_parameter[15];
              tmpS=_IndividualImpactSlopeByInfecVac[subject][0][ser];
              probCurrentTiter+=log(dgamma(tmpS,par1,par2));
          }
      }
      if(_individualLevelEffectsAllByInfecVac==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==1){
          for (ser=0;ser<_nSerAssay;ser++){
              par1=_parameter[44]*_parameter[44]/_parameter[13];
              par2=_parameter[44]/_parameter[13];
              tmpT=_IndividualImpactTiterByInfecVac[subject][0][ser];
              probCurrentTiter+=log(dgamma(tmpT,par1,par2));
              
              meanR=_parameter[46];
              par1=meanR*meanR/_parameter[14];
              par2=meanR/_parameter[14];
              tmpR=_IndividualImpactTempTiterByInfecVac[subject][0][ser];
              probCurrentTiter+=log(dgamma(tmpR,par1,par2));
              
              par1=_parameter[48]*_parameter[48]/_parameter[15];
              par2=_parameter[48]/_parameter[15];
              tmpS=_IndividualImpactSlopeByInfecVac[subject][0][ser];
              probCurrentTiter+=log(dgamma(tmpS,par1,par2));
          }
      }
      
    double _oldProbaSubject=_probaTiterIndividual[subject];
    removeAsymptomaticInfection(subject, infecNo);
    
    if(_individualLevelEffectsAllByInfec==1){
        for (inf=0;inf<currentTotInf-1;inf++){
          for (ser=0;ser<_nSerAssay;ser++){
            _oldTiters[inf][ser]=_IndividualImpactTiterByInfec[subject][inf][ser];
            _oldTempTiters[inf][ser]=_IndividualImpactTempTiterByInfec[subject][inf][ser];
            _oldSlopes[inf][ser]=_IndividualImpactSlopeByInfec[subject][inf][ser];
          }
        }
    }
    
    double probFutureTiter=0;
    double indR, indS, indT;
    int minRise=0;
      if(_individualLevelEffectsAllByInfec==1&_vaccineSubject[subject]==0){
        for (inf=0;inf<currentTotInf-1;inf++){
          for (ser=0;ser<_nSerAssay;ser++){
            par1=_parameter[11]*_parameter[11]/_parameter[13];
            par2=_parameter[11]/_parameter[13];
            indT=max(1e-10,rgamma(par1,par2));
            _IndividualImpactTiterByInfec[subject][inf][ser]=indT;
            probFutureTiter+=log(dgamma(indT,par1,par2));
            
            meanR=_parameter[0];
            par1=meanR*meanR/_parameter[14];
            par2=meanR/_parameter[14];
            indR=max(1e-10,rgamma(par1,par2));
            _IndividualImpactTempTiterByInfec[subject][inf][ser]=indR;
            probFutureTiter+=log(dgamma(indR,par1,par2));
            
            par1=_parameter[2]*_parameter[2]/_parameter[15];
            par2=_parameter[2]/_parameter[15];
            indS=max(_minSlope,rgamma(par1,par2));
            _IndividualImpactSlopeByInfec[subject][inf][ser]=indS;
            probFutureTiter+=log(dgamma(indS,par1,par2));
          }
        }
    }
      if(_individualLevelEffectsAllByInfec==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==1){
          for (inf=0;inf<currentTotInf-1;inf++){
              for (ser=0;ser<_nSerAssay;ser++){
                  double mean=_parameter[54];//*_parameter[11];
                  par1=mean*mean/_parameter[55];
                  par2=mean/_parameter[55];
                  indT=max(1e-10,rgamma(par1,par2));
                  _IndividualImpactTiterByInfec[subject][inf][ser]=indT;
                  probFutureTiter+=log(dgamma(indT,par1,par2));
                  
                  meanR=_parameter[56];//*_parameter[0];
                  par1=meanR*meanR/_parameter[57];
                  par2=meanR/_parameter[57];
                  indR=max(1e-10,rgamma(par1,par2));
                  _IndividualImpactTempTiterByInfec[subject][inf][ser]=indR;
                  probFutureTiter+=log(dgamma(indR,par1,par2));
                  
                  mean=_parameter[58];//*_parameter[2];
                  par1=mean*mean/_parameter[59];
                  par2=mean/_parameter[59];
                  indS=max(_minSlope,rgamma(par1,par2));
                  _IndividualImpactSlopeByInfec[subject][inf][ser]=indS;
                  probFutureTiter+=log(dgamma(indS,par1,par2));
              }
          }
      }
      if(_individualLevelEffectsAllByInfec==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==0){
          for (inf=0;inf<currentTotInf-1;inf++){
              for (ser=0;ser<_nSerAssay;ser++){
                  double mean=_parameter[64];//*_parameter[11];
                  par1=mean*mean/_parameter[65];
                  par2=mean/_parameter[65];
                  indT=max(1e-10,rgamma(par1,par2));
                  _IndividualImpactTiterByInfec[subject][inf][ser]=indT;
                  probFutureTiter+=log(dgamma(indT,par1,par2));
                  
                  meanR=_parameter[66];//*_parameter[0];
                  par1=meanR*meanR/_parameter[67];
                  par2=meanR/_parameter[67];
                  indR=max(1e-10,rgamma(par1,par2));
                  _IndividualImpactTempTiterByInfec[subject][inf][ser]=indR;
                  probFutureTiter+=log(dgamma(indR,par1,par2));
                  
                  mean=_parameter[68];//*_parameter[2];
                  par1=mean*mean/_parameter[69];
                  par2=mean/_parameter[69];
                  indS=max(_minSlope,rgamma(par1,par2));
                  _IndividualImpactSlopeByInfec[subject][inf][ser]=indS;
                  probFutureTiter+=log(dgamma(indS,par1,par2));
              }
          }
      }
      double tmpRise;
      int minRise2=0;
      if(_individualLevelEffectsAllByInfecVac==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==0){
          for (ser=0; ser<_nSerAssay;ser++){
              par1=_parameter[34]*_parameter[34]/_parameter[35];
              par2=_parameter[34]/_parameter[35];
              tmpT=max(1e-10,rgamma(par1,par2));
              _IndividualImpactTiterByInfecVac[subject][0][ser]=tmpT;
              probFutureTiter+=log(dgamma(tmpT,par1,par2));

              par1=_parameter[36]*_parameter[36]/_parameter[14];
              par2=_parameter[36]/_parameter[14];
              tmpR=max(1e-10,rgamma(par1,par2));
              _IndividualImpactTempTiterByInfecVac[subject][0][ser]=tmpR;
              probFutureTiter+=log(dgamma(tmpR,par1,par2));

              par1=_parameter[38]*_parameter[38]/_parameter[15];
              par2=_parameter[38]/_parameter[15];
              tmpS=max(_minSlope,rgamma(par1,par2));
              _IndividualImpactSlopeByInfecVac[subject][0][ser]=tmpS;
              probFutureTiter+=log(dgamma(tmpS,par1,par2));
              
              tmpRise=tmpR*exp(-tmpS*365);
              if(tmpRise<_minTempTiterRise){minRise2=1;}
          }
      }
      if(_individualLevelEffectsAllByInfecVac==1&_vaccineSubject[subject]==1&_baselineNaive[subject]==1){
          for (ser=0; ser<_nSerAssay;ser++){
              par1=_parameter[44]*_parameter[44]/_parameter[13];
              par2=_parameter[44]/_parameter[13];
              tmpT=max(1e-10,rgamma(par1,par2));
              _IndividualImpactTiterByInfecVac[subject][0][ser]=tmpT;
              probFutureTiter+=log(dgamma(tmpT,par1,par2));
              
              par1=_parameter[46]*_parameter[46]/_parameter[14];
              par2=_parameter[46]/_parameter[14];
              tmpR=max(1e-10,rgamma(par1,par2));
              _IndividualImpactTempTiterByInfecVac[subject][0][ser]=tmpR;
              probFutureTiter+=log(dgamma(tmpR,par1,par2));
              
              par1=_parameter[48]*_parameter[48]/_parameter[15];
              par2=_parameter[48]/_parameter[15];
              tmpS=max(_minSlope,rgamma(par1,par2));
              _IndividualImpactSlopeByInfecVac[subject][0][ser]=tmpS;
              probFutureTiter+=log(dgamma(tmpS,par1,par2));
              
              tmpRise=tmpR*exp(-tmpS*365);
              if(tmpRise<_minTempTiterRise){minRise2=1;}
          }
      }
      
    computeProbaIndividual(subject,0,0,0);
    
    if(specificSubjectInd==0){minRise2=checkMinRise(subject);}
      
    double _newProbaSubject=_probaTiterIndividual[subject];
    
    int minDay=_minDayStudy[subject]+1;
    int maxDay=_maxDayStudy[subject]-1;
    double totProbEscape=-_parameter[12]*(_epiCDF[Date-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]-_epiCDF[minDay-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]);
    double probDay=totProbEscape+log(1-exp(-_parameter[12]*_epiPDF[Date-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]));
      
    double logProbaAddCase=probDay+probCurrentTiter+probDecayOld;
    double logProbaRemoveCase=log(1)+probFutureTiter+probDecayNew;
    double logProposal=logProbaAddCase-logProbaRemoveCase;
    double _newLogLikGlobal=_logLikGlobal-_oldProbaSubject+_newProbaSubject;
    
    if(log(runif())<_newLogLikGlobal-_logLikGlobal+logProposal&minRise==0&minRise2==0)
    {   _logLikGlobal=_newLogLikGlobal;
      (*accepted)++;
    }
    else{
        if(_individualLevelEffectsAllByInfec==1){
          for (inf=0;inf<currentTotInf-1;inf++){
            for(ser=0;ser<_nSerAssay;ser++){
              _IndividualImpactTiterByInfec[subject][inf][ser]=_oldTiters[inf][ser];
              _IndividualImpactTempTiterByInfec[subject][inf][ser]=_oldTempTiters[inf][ser];
              _IndividualImpactSlopeByInfec[subject][inf][ser]=_oldSlopes[inf][ser];
            }
          }
        }
        if(_individualLevelEffectsAllByInfecVac==1&_vaccineSubject[subject]==1){
            for (ser=0;ser<_nSerAssay;ser++){
                _IndividualImpactTiterByInfecVac[subject][0][ser]=_oldTitersVac[0][ser];
                _IndividualImpactTempTiterByInfecVac[subject][0][ser]=_oldTempTitersVac[0][ser];
                _IndividualImpactSlopeByInfecVac[subject][0][ser]=_oldSlopesVac[0][ser];
            }
        }
      for (ser=0;ser<_nSerAssay;ser++){
        _baselineTiter[subject][ser]=_oldBaselines[ser];
          if(_individualLevelEffectsSlowDecay==1){_IndividualImpactSlowDecay[subject][ser]=_oldSlowDecay[ser];}
      }
      addAsymptomaticInfection(subject, Date, Serotype);
      computeProbaIndividual(subject,0,0,0);
    }
    (*proposed)++;
  }
}

////=========== update asymptomatic date ============
void updateAsymptomaticDate(int nbIter,double *accepted,double *proposed,int specificSubjectInd=0,int specificSubject=-999)
{   int iter,subject, indicator, noInfec,i;
  int a,b;
  for(iter=0;iter<nbIter;iter++)
  {   if(specificSubjectInd==1){subject=specificSubject;}else{
      subject=int(runif()*_numberOfSubject);
        }
    
    int symp=_freqSymptomaticSerKnown[subject];
    int symp2=_freqSymptomaticSerUnknown[subject];
    
    int currentNoInf=_numberAsymptomatic[subject];
    if (currentNoInf==0)continue;
    
    int infecNo=int(runif()*(currentNoInf));
    int ser=_serotypeAsymptomatic[subject][infecNo];
    int oldDate=_dayAsymptomatic[subject][infecNo];
    double _oldProbaSubject=_probaTiterIndividual[subject];
    
    int minDay=_minDayStudy[subject]+1;
    int maxDay=_maxDayStudy[subject]-1;
    double totProbEscape=-_parameter[12]*(_epiCDF[oldDate-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]-_epiCDF[minDay-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]);
    double oldProbDay=totProbEscape+log(1-exp(-_parameter[12]*_epiPDF[oldDate-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]));
    
    int newDate;
    double newProbDay;
    drawDateAsymptomatic(minDay, maxDay, &newDate, &newProbDay);
    
    if(newDate+_parameter[1]-_incubationPeriod-_delayRise<_minDayInfected)continue;
    
    indicator=0;
    noInfec=_freqSymptomaticSerKnown[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++)
        {   if(abs(_daySymptomaticSerKnown[subject][i]-newDate)<maxTimeBetweenInfections){indicator=1;}
        }
    }
    noInfec=_freqSymptomaticSerUnknown[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++)
        {   if(abs(_daySymptomaticSerUnknown[subject][i]-newDate)<maxTimeBetweenInfections){indicator=1;}
        }
    }
    noInfec=_numberAsymptomatic[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++)
        {   if((i!=infecNo)&(abs(_dayAsymptomatic[subject][i]-newDate)<maxTimeBetweenInfections)){indicator=1;}
        }
    }
    if(_vaccineSubject[subject]==1){
      if(abs(_dayVacc[subject][0]-newDate)<_minTimeInfecPostVaccine){indicator=1;}
    }
    if(indicator==1) continue;
      
    int sero;
    for (sero=0;sero<_nSerAssay;sero++){
      _tmpTvec[sero]=_IndividualImpactTiterByInfec[subject][symp+symp2+infecNo][sero];
      _tmpRvec[sero]=_IndividualImpactTempTiterByInfec[subject][symp+symp2+infecNo][sero];
      _tmpSvec[sero]=_IndividualImpactSlopeByInfec[subject][symp+symp2+infecNo][sero];
    }
    
    removeAsymptomaticInfection(subject, infecNo);
    addAsymptomaticInfection(subject,newDate,ser);
      
      computeProbaIndividual(subject,0,0,0);
      
      int minRise=0;
      minRise=checkMinRise(subject,1);
      
    
    double _newProbaSubject=_probaTiterIndividual[subject];
    double _newLogLikGlobal=_logLikGlobal-_oldProbaSubject+_newProbaSubject;
    
    double logProbaAddCase=newProbDay;
    double logProbaRemoveCase=oldProbDay;
    double logProposal=logProbaRemoveCase-logProbaAddCase;
    
    if(log(runif())<_newLogLikGlobal-_logLikGlobal+logProposal&minRise==0)
    {   _logLikGlobal=_newLogLikGlobal;
      (*accepted)++;
    }
    else{
      removeAsymptomaticInfection(subject, 0);
      addAsymptomaticInfection(subject,oldDate,ser);
      computeProbaIndividual(subject,0,0,0);
    }
    (*proposed)++;
  }
}


////=========== update asymptomatic date ============
void updateAsymptomaticDateSmallChange(int nbIter,double *accepted,double *proposed)
{   int iter,subject, indicator, noInfec,i;
  int a,b;
  for(iter=0;iter<nbIter;iter++)
  {   subject=int(runif()*_numberOfSubject);
    
    int symp=_freqSymptomaticSerKnown[subject];
    int symp2=_freqSymptomaticSerUnknown[subject];
    
    int currentNoInf=_numberAsymptomatic[subject];
    if (currentNoInf==0) continue;
    
    int infecNo=int(runif()*(currentNoInf));
    int ser=_serotypeAsymptomatic[subject][infecNo];
    int oldDate=_dayAsymptomatic[subject][infecNo];
    double _oldProbaSubject=_probaTiterIndividual[subject];
    
    int minDay=oldDate-30;
    int maxDay=oldDate+30;
    
    double totProbEscape=-_parameter[12]*(_epiCDF[oldDate-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]-_epiCDF[minDay-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]);
    double oldProbDay=totProbEscape+log(1-exp(-_parameter[12]*_epiPDF[oldDate-_minDayInfected-int(_incubationPeriod)-int(_delayRise)]));
    int newDate;
    double newProbDay;
    drawDateAsymptomatic(minDay, maxDay, &newDate, &newProbDay);
    
    if(newDate<_minDayStudy[subject]+1)continue;
    if(newDate>_maxDayStudy[subject]-1)continue;
    if(newDate+_parameter[1]-_incubationPeriod-_delayRise<_minDayInfected)continue;
    //        if(_baselineNaive[subject]==1&newDate<minDay+365) continue;
    
    indicator=0;
    noInfec=_freqSymptomaticSerKnown[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++)
        {   if(abs(_daySymptomaticSerKnown[subject][i]-newDate)<maxTimeBetweenInfections){indicator=1;}
        }
    }
    noInfec=_freqSymptomaticSerUnknown[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++)
        {   if(abs(_daySymptomaticSerUnknown[subject][i]-newDate)<maxTimeBetweenInfections){indicator=1;}
        }
    }
    noInfec=_numberAsymptomatic[subject];
    if(noInfec>0){
      for (i=0; i<noInfec;i++)
        {   if((i!=infecNo)&(abs(_dayAsymptomatic[subject][i]-newDate)<maxTimeBetweenInfections)){indicator=1;}
        }
    }
    if(_vaccineSubject[subject]==1){
      if(abs(_dayVacc[subject][0]-newDate)<_minTimeInfecPostVaccine){indicator=1;}
    }
    if(indicator==1) continue;
      
    int sero;
    for (sero=0;sero<_nSerAssay;sero++){
      _tmpTvec[sero]=_IndividualImpactTiterByInfec[subject][symp+symp2+infecNo][sero];
      _tmpRvec[sero]=_IndividualImpactTempTiterByInfec[subject][symp+symp2+infecNo][sero];
      _tmpSvec[sero]=_IndividualImpactSlopeByInfec[subject][symp+symp2+infecNo][sero];
    }
    
    double oldPropR, newPropR, meanR, par1,par2, tmpR;
    removeAsymptomaticInfection(subject, infecNo);
    addAsymptomaticInfection(subject,newDate,ser);
    
      computeProbaIndividual(subject,0,0,0);
      
      int minRise=0;
      minRise=checkMinRise(subject,1);
      
    
    double _newProbaSubject=_probaTiterIndividual[subject];
    double _newLogLikGlobal=_logLikGlobal-_oldProbaSubject+_newProbaSubject;
    
    double logProbaAddCase=newProbDay;
    double logProbaRemoveCase=oldProbDay;
    double logProposal=logProbaRemoveCase-logProbaAddCase;
    
    if(log(runif())<_newLogLikGlobal-_logLikGlobal+logProposal&minRise==0)
    {   _logLikGlobal=_newLogLikGlobal;
      (*accepted)++;
    }
    else{
      removeAsymptomaticInfection(subject, 0);
      addAsymptomaticInfection(subject,oldDate,ser);
      computeProbaIndividual(subject,0,0,0);
    }
    (*proposed)++;
  }
}


















//========== RUN MCMC==============
void runMCMC(const char *outputParameters,const char *outputIndividualTiterSim,const char *outputIndividualImpactParameters,const char *outputIndividualInfectionHistories, const char *outputIndividualTiterTrajectories, int pas,int pasPrinting,int pasParameterPrinting, int numberOfIterationTot)
{
  ofstream outputParametersFile(outputParameters);
  ofstream outputIndividualTiterSimFile,outputIndividualImpactParametersFile,outputIndividualInfectionHistoriesFile, outputIndividualTiterTrajectoriesFile;
  
  outputIndividualTiterSimFile.open(outputIndividualTiterSim);
  outputIndividualImpactParametersFile.open(outputIndividualImpactParameters);
  outputIndividualInfectionHistoriesFile.open(outputIndividualInfectionHistories);
  outputIndividualTiterTrajectoriesFile.open(outputIndividualTiterTrajectories);
  
  int parameterNumber,iteration,i;
  
  double *numberOfMoveAccepted=new double[_numberOfParameters];
  double *numberOfMoveProposed=new double[_numberOfParameters];
  
  for(i=0;i<_numberOfParameters;i++)
  {numberOfMoveAccepted[i]=0;
    numberOfMoveProposed[i]=0;}
  
  cout<<"=============MCMC=================\n";
  
  double accepted=0;
  double acceptedIndividualLevelDelay=0;
  double acceptedSerotypeSymptomaticAugmented=0;
  double acceptedSerotypeAsymptomaticAugmented=0;
  double acceptedIndividualLevelAllInfAllPar=0;
  double acceptedIndividualLevelAllInfAllParVac=0;
  computeProbaAll(&_logLikGlobal,0,0,0);
  cout<<"check1\n";
  int numberOfIteration=int(numberOfIterationTot/pas);
  for(iteration=0;iteration<numberOfIteration;iteration++)
  {
    int addAsymps, subject, noInfec, i;
    if(_simThenEstimate==1&int(iteration/_iterationPerSim)*_iterationPerSim==iteration){
      
      for (i=0;i<_numberOfParameters;i++){
        _parameter[i]=_initialParameter[i];
      }
      if(_individualLevelEffectsD==1){
        for(subject=0;subject<_numberOfSubject;subject++){
          _IndividualImpactDelay[subject]=0;
        }}
      double probaInfOb=0.3; //0.3
      double probaSerOb=1.0; //0.8
      removeAllInfectionsForSim();
      _numberMissingSymptomaticSerotype=0;
      _numberNotMissingSymptomaticSerotype=0;
      _numberOfSubject=_numberOfSubjectSim;
      cout<<"check1\n";
      double probNaiveBaseline=0.2;
      int ser;
      for(i=0;i<_numberOfSubject;i++){
        for (ser=0;ser<_nSerAssay;ser++){
          _baselineTiter[i][ser]=2+runif()*1;
          if(_baselineNaive[i]==1){_baselineTiter[i][ser]=-999;}
        }
      }
      updateEpiPDF();
      addAsymptomaticInfectionForSim(probaInfOb,probaSerOb);
//    //Add vaccine response
        int minRise;
        for (subject=0;subject<_numberOfSubject;subject++){
    minRise=1;
        while(minRise==1){
        if(_vaccineSubject[subject]==1&_individualLevelEffectsAllByInfecVac==1&_baselineNaive[subject]==0){
          for (i=0;i<_numberVaccines[subject];i++){
            for (ser=0;ser<_nSerAssay;ser++){
              double par1=_parameter[34]*_parameter[34]/_parameter[35];
              double par2=_parameter[34]/_parameter[35];
              _IndividualImpactTiterByInfecVac[subject][i][ser]=max(1e-10,rgamma(par1,par2));
              
              par1=_parameter[36]*_parameter[36]/_parameter[14];
              par2=_parameter[36]/_parameter[14];
              _IndividualImpactTempTiterByInfecVac[subject][i][ser]=max(1e-10,rgamma(par1,par2));
              
              par1=_parameter[38]*_parameter[38]/_parameter[15];
              par2=_parameter[38]/_parameter[15];
              _IndividualImpactSlopeByInfecVac[subject][i][ser]=max(1e-10,rgamma(par1,par2));
            }
          }
        }
        if(_vaccineSubject[subject]==1&_individualLevelEffectsAllByInfecVac==1&_baselineNaive[subject]==1){
            for (i=0;i<_numberVaccines[subject];i++){
                for (ser=0;ser<_nSerAssay;ser++){
                    double par1=_parameter[44]*_parameter[44]/_parameter[13];
                    double par2=_parameter[44]/_parameter[13];
                    _IndividualImpactTiterByInfecVac[subject][i][ser]=max(1e-10,rgamma(par1,par2));
                    
                    par1=_parameter[46]*_parameter[46]/_parameter[14];
                    par2=_parameter[46]/_parameter[14];
                    _IndividualImpactTempTiterByInfecVac[subject][i][ser]=max(1e-10,rgamma(par1,par2));
                    
                    par1=_parameter[48]*_parameter[48]/_parameter[15];
                    par2=_parameter[48]/_parameter[15];
                    _IndividualImpactSlopeByInfecVac[subject][i][ser]=max(1e-10,rgamma(par1,par2));
                }
            }
        }
        minRise=checkMinRise(subject,0);
        }
    }
      
        
      // Record parameters
      cout<<iteration<<":"<<_logLikGlobal/1000<<"; ";
      outputParametersFile<<iteration<<" "<<_logLikGlobal<<" ";
      for(i=0;i<8;i++)
      {   outputParametersFile<<_parameter[i]<<" ";
        cout<<_parameter[i]<<" ; ";}
      cout<<"\n";
      for(i=8;i<16;i++)
      {   outputParametersFile<<_parameter[i]<<" ";
        cout<<_parameter[i]<<" ; ";}
      cout<<"\n";
      for(i=16;i<_numberOfParameters;i++)
      {   outputParametersFile<<_parameter[i]<<" ";
        cout<<_parameter[i]<<" ; ";}
      outputParametersFile<<_nAsymptomatic<<" ";
      outputParametersFile<<"\n"<<flush;
      
      computeProbaAll(&_logLikGlobal,0,0,0);
      MakeSimulatedTiters();

        
      int day;
      for (subject=0; subject<_numberOfSubject; subject++)
      {   outputIndividualTiterTrajectoriesFile<<-999<<" "<<iteration<<" ";
        outputIndividualTiterTrajectoriesFile<<-888<<" "<<subject<<" ";
        for(ser=0;ser<_nSerAssay;ser++){
          for(day=_minDayInfected;day<_maxDay;day++){
            outputIndividualTiterTrajectoriesFile<< computeTiterIndividualAtDay(subject,day,ser)<<" ";
          }
        }
        outputIndividualTiterTrajectoriesFile<<"\n"<<flush;
      }
      
      cout<<"NumberSymptomaticSerKnown: "<<_numberNotMissingSymptomaticSerotype<<"\n";
      cout<<"NumberSymptomaticSerUnknown: "<<_numberMissingSymptomaticSerotype<<"\n";
      cout<<"TrueNumberAsymptomatic: "<<_nAsymptomatic<<"\n";
      
      int aa,bb;
      aa=_numberAsymptomatic[9];
      for (subject=0;subject<_numberOfSubject;subject++){
        noInfec=_numberAsymptomatic[subject];
        _IndividualImpactDelay[subject]=0;
        if(noInfec>0){
          for (i=0; i<noInfec;i++){
            aa=_numberAsymptomatic[9];
            removeAsymptomaticInfection(subject, 0);
            bb=_numberAsymptomatic[9];
          }
        }
      }
      maxDelaySympInfec();
    }
    computeProbaAll(&_logLikGlobal,0,0,0);
      
      cout<<"LL2: "<<_logLikGlobal<<"\n";
    
      
///Give model a helping hand to find starting infections
      double acceptedRemoveInfection=0;
      double acceptedAddInfection=0;
      double proposeRemoveInfection=0;
      double proposeAddInfection=0;
      double acceptedChangeDate=0;
      double proposedChangeDate=0;
      double acceptedChangeDateSmallChange=0;
      double proposedChangeDateSmallChange=0;
      double acceptedBaseline=0;
      double acceptedSlowDecay=0;
      int kRJMCMC;
      int noRJMCMC;
      int minRiseInd;
      if(iteration==0){
          for (i=0;i<_numberOfSubject;i++){
              minRiseInd=checkMinRise(i,1);
              while(minRiseInd==1){
                  if(_individualLevelEffectsD==1){acceptedIndividualLevelDelay=updateIndividualImpactDelay(10,1,i);}
                  if(_individualLevelEffectsAllByInfec==1){acceptedIndividualLevelAllInfAllPar=updateIndividualImpactAllInfAllPar(10,1,i);}
                  acceptedIndividualLevelAllInfAllParVac=updateIndividualImpactAllInfAllParVac(10,1,i);
                  for(kRJMCMC=0;kRJMCMC<10;kRJMCMC++){
                      if(runif()<0.5){
                          RJMCMCAddAsymptomaticInfection(&acceptedAddInfection,&proposeAddInfection,1,i);
                      }else{
                          RJMCMCRemoveAsymptomaticInfection(&acceptedRemoveInfection,&proposeRemoveInfection,1,i);
                      }
                  }
                  updateAsymptomaticDate(10,&acceptedChangeDate,&proposedChangeDate,1,i);
                  acceptedBaseline=updateBaselineTiter(10,1,i,0);
                  minRiseInd=checkMinRise(i,1);
              }
          }
      }
      
      int check;
      
      double tmpLL,aa;
    for(i=0;i<pas;i++){
      int selectedParameter;
      for(selectedParameter=0;selectedParameter<_numberOfSelectedParameter;selectedParameter++)
      {   parameterNumber=_idOfSelectedParameter[selectedParameter];
         
          if(_stopEstimatingPars==1&iteration>_nbIterRJMCMCpostDelay&(parameterNumber==0|parameterNumber==2|parameterNumber==3|parameterNumber==11|parameterNumber==13|parameterNumber==14|parameterNumber==15|parameterNumber==19|parameterNumber==31|parameterNumber==32|parameterNumber==34|parameterNumber==35|parameterNumber==36|parameterNumber==37|parameterNumber==38|parameterNumber==39)) continue;
          
        accepted=updateParameter(parameterNumber);
          
        numberOfMoveAccepted[parameterNumber]+=accepted;
        numberOfMoveProposed[parameterNumber]++;
      }
    }
      

      
//    // Update serotype of symptomatic infections where not known
    if(_numberMissingSymptomaticSerotype>0&_nSerAssay>1){
      acceptedSerotypeSymptomaticAugmented=updateSymptomaticSerotype(_numberMissingSymptomaticSerotype);
    }

    // Update serotype of asymptomatic infections
    if(_nSerAssay>1){
      acceptedSerotypeAsymptomaticAugmented=updateAsymptomaticSerotype(200);
    }

      
    acceptedRemoveInfection=0;
    acceptedAddInfection=0;
    proposeRemoveInfection=0;
    proposeAddInfection=0;

    if(iteration<_nbIterRJMCMCDelay){noRJMCMC=_nbIterRJMCMC;}else{noRJMCMC=_nbIterRJMCMCpostDelay;}
    for(kRJMCMC=0;kRJMCMC<noRJMCMC;kRJMCMC++){
      if(runif()<0.5){
        RJMCMCAddAsymptomaticInfection(&acceptedAddInfection,&proposeAddInfection);
      }else{
        RJMCMCRemoveAsymptomaticInfection(&acceptedRemoveInfection,&proposeRemoveInfection);
      }
    }
      

    // Change dates of asymptomatic infections
      if(iteration>_nbIterRJMCMCDelay){
    acceptedChangeDate=0;
    proposedChangeDate=0;
    updateAsymptomaticDate(200,&acceptedChangeDate,&proposedChangeDate);
      acceptedChangeDateSmallChange=0;
    proposedChangeDateSmallChange=0;
    updateAsymptomaticDateSmallChange(200,&acceptedChangeDateSmallChange,&proposedChangeDateSmallChange);
      }

    // Impact individual level effects
    if(_individualLevelEffectsAllByInfec==1){acceptedIndividualLevelAllInfAllPar=updateIndividualImpactAllInfAllPar(300);}
    
    if(_individualLevelEffectsAllByInfecVac==1){acceptedIndividualLevelAllInfAllParVac=updateIndividualImpactAllInfAllParVac(300);}

    if(_individualLevelEffectsD==1){acceptedIndividualLevelDelay=updateIndividualImpactDelay(100);}
      
    acceptedBaseline=updateBaselineTiter(300,0,-999,1);
      
      if(_individualLevelEffectsSlowDecay==1){acceptedSlowDecay=updateSlowDecayTiters(300);}
      // Record parameters and simulations
    if (int(iteration/pasParameterPrinting)*pasParameterPrinting==iteration){
      cout<<iteration<<":"<<_logLikGlobal/1000<<"; ";
      outputParametersFile<<iteration<<" "<<_logLikGlobal<<" ";
      for(i=0;i<5;i++)
      {   outputParametersFile<<_parameter[i]<<" ";
        cout<<_parameter[i]<<" ; ";}
      cout<<"\n";
      for(i=5;i<10;i++)
      {   outputParametersFile<<_parameter[i]<<" ";
        cout<<_parameter[i]<<" ; ";}
      cout<<"\n";
      for(i=10;i<15;i++)
      {   outputParametersFile<<_parameter[i]<<" ";
        cout<<_parameter[i]<<" ; ";}
      cout<<"\n";
      for(i=15;i<20;i++)
      {   outputParametersFile<<_parameter[i]<<" ";
        cout<<_parameter[i]<<" ; ";}
      cout<<"\n";
      for(i=20;i<25;i++)
      {   outputParametersFile<<_parameter[i]<<" ";
        cout<<_parameter[i]<<" ; ";}
        cout<<"\n";
      for(i=25;i<30;i++)
        {   outputParametersFile<<_parameter[i]<<" ";
            cout<<_parameter[i]<<" ; ";}
        cout<<"\n";
        for(i=30;i<35;i++)
        {   outputParametersFile<<_parameter[i]<<" ";
            cout<<_parameter[i]<<" ; ";}
        cout<<"\n";
        for(i=35;i<40;i++)
        {   outputParametersFile<<_parameter[i]<<" ";
            cout<<_parameter[i]<<" ; ";}
        cout<<"\n";
        for(i=40;i<45;i++)
        {   outputParametersFile<<_parameter[i]<<" ";
            cout<<_parameter[i]<<" ; ";}
        cout<<"\n";
        for(i=45;i<50;i++)
        {   outputParametersFile<<_parameter[i]<<" ";
            cout<<_parameter[i]<<" ; ";}
        cout<<"\n";
        for(i=50;i<55;i++)
        {   outputParametersFile<<_parameter[i]<<" ";
            cout<<_parameter[i]<<" ; ";}
        cout<<"\n";
        for(i=55;i<_numberOfParameters;i++)
        {   outputParametersFile<<_parameter[i]<<" ";
            cout<<_parameter[i]<<" ; ";}
        cout<<"\n";
      outputParametersFile<<_nAsymptomatic<<" ";
      outputParametersFile<<"\n"<<flush;
      cout<<"\n";
      cout<<"NumberAsymptomatic: "<<_nAsymptomatic<<"\n";
      cout<<"ChangedAsympSerotype: "<<acceptedSerotypeAsymptomaticAugmented<<"\n";
      cout<<"ChangedSympSerotype: "<<acceptedSerotypeSymptomaticAugmented<<"\n";
      if (_individualLevelEffectsAllByInfec==1){cout<<"ChangedIndividualLevelByInfAllPar: "<<acceptedIndividualLevelAllInfAllPar<<"\n";}
      if (_individualLevelEffectsAllByInfecVac==1){cout<<"ChangedIndividualLevelByInfAllParVac: "<<acceptedIndividualLevelAllInfAllParVac<<"\n";}
      if (_individualLevelEffectsD==1){cout<<"ChangedIndividualLevelDelay: "<<acceptedIndividualLevelDelay<<"\n";}
      cout<<"ChangedBaseline: "<<acceptedBaseline<<"\n";
        if(_individualLevelEffectsSlowDecay==1){cout<<"ChangedSlowDecay: "<<acceptedSlowDecay<<"\n";}
      cout<<"ChangedAsymptomaticDate: "<<acceptedChangeDate<<"/"<<proposedChangeDate<<"\n";
      cout<<"ChangedAsymptomaticDateSmallChange: "<<acceptedChangeDateSmallChange<<"/"<<proposedChangeDateSmallChange<<"\n";
      cout<<"RemoveInfection: "<<acceptedRemoveInfection<<"/"<<proposeRemoveInfection<<"__AddInfection: "<<acceptedAddInfection<<"/"<<proposeAddInfection;
      cout<<"\n";
    }
    
    
    if (int(iteration/pasPrinting)*pasPrinting==iteration)
    {
      
      int maxCol=14*_numberOfSubject;
      int i,j, subject, noInfec, noTests, ser;
      int counter;
           
    
    //Update Individual Impacts
    for (subject=0; subject<_numberOfSubject; subject++){
            outputIndividualImpactParametersFile<<-999<<" "<<iteration<<" ";
            outputIndividualImpactParametersFile<<-888<<" "<<subject<<" ";
            outputIndividualImpactParametersFile<<-666<<" "<<_IndividualImpactDelay[subject]<<" ";
            outputIndividualImpactParametersFile<<-333<<" "<<_baselineTiter[subject][0]<<" ";
            outputIndividualImpactParametersFile<<_baselineTiter[subject][1]<<" ";
            outputIndividualImpactParametersFile<<_baselineTiter[subject][2]<<" ";
            outputIndividualImpactParametersFile<<_baselineTiter[subject][3]<<" ";
            outputIndividualImpactParametersFile<<-222<<" "<<_baselineNaive[subject]<<" ";
        for(ser=0;ser<_nSerAssay;ser++){
            outputIndividualImpactParametersFile<<_IndividualImpactSlowDecay[subject][ser]<<" ";
        }
        if(_vaccineSubject[subject]==1){
             for(ser=0;ser<_nSerAssay;ser++){
               outputIndividualImpactParametersFile<<_IndividualImpactTiterByInfecVac[subject][0][ser]<<" ";
               outputIndividualImpactParametersFile<<_IndividualImpactTempTiterByInfecVac[subject][0][ser]<<" ";
               outputIndividualImpactParametersFile<<_IndividualImpactSlopeByInfecVac[subject][0][ser]<<" ";
             }
         }
         if(_vaccineSubject[subject]==0){
           for(ser=0;ser<_nSerAssay;ser++){
             outputIndividualImpactParametersFile<<-1234<<" ";
             outputIndividualImpactParametersFile<<-1234<<" ";
             outputIndividualImpactParametersFile<<-1234<<" ";
           }
         }
         counter=0;
        noInfec=_freqSymptomaticSerKnown[subject];
           if(noInfec>0){
             for (j=0; j<noInfec;j++)
             {   outputIndividualImpactParametersFile<<_daySymptomaticSerKnown[subject][j]<<" ";
                 outputIndividualImpactParametersFile<<_serotypeSymptomatic[subject][j]<<" ";
                 outputIndividualImpactParametersFile<<-1111<<" ";
                 for(ser=0;ser<_nSerAssay;ser++){
                   outputIndividualImpactParametersFile<<_IndividualImpactTiterByInfec[subject][counter][ser]<<" ";
                   outputIndividualImpactParametersFile<<_IndividualImpactTempTiterByInfec[subject][counter][ser]<<" ";
                   outputIndividualImpactParametersFile<<_IndividualImpactSlopeByInfec[subject][counter][ser]<<" ";
                 }
                 counter++;
             }
           }
           noInfec=_freqSymptomaticSerUnknown[subject];
           if(noInfec>0){
             for (j=0; j<noInfec;j++)
             {   outputIndividualImpactParametersFile<<_daySymptomaticSerUnknown[subject][j]<<" ";
                 outputIndividualImpactParametersFile<<_augmentedSerotypeSymptomatic[subject][j]<<" ";
                 outputIndividualImpactParametersFile<<-2222<<" ";
                 for(ser=0;ser<_nSerAssay;ser++){
                   outputIndividualImpactParametersFile<<_IndividualImpactTiterByInfec[subject][counter][ser]<<" ";
                   outputIndividualImpactParametersFile<<_IndividualImpactTempTiterByInfec[subject][counter][ser]<<" ";
                   outputIndividualImpactParametersFile<<_IndividualImpactSlopeByInfec[subject][counter][ser]<<" ";
                 }
                 counter++;
             }
           }
           noInfec=_numberAsymptomatic[subject];
           if(noInfec>0){
             for (j=0; j<noInfec;j++)
             {   outputIndividualImpactParametersFile<<_dayAsymptomatic[subject][j]<<" ";
                 outputIndividualImpactParametersFile<<_serotypeAsymptomatic[subject][j]<<" ";
                 outputIndividualImpactParametersFile<<-3333<<" ";
                 for(ser=0;ser<_nSerAssay;ser++){
                   outputIndividualImpactParametersFile<<_IndividualImpactTiterByInfec[subject][counter][ser]<<" ";
                   outputIndividualImpactParametersFile<<_IndividualImpactTempTiterByInfec[subject][counter][ser]<<" ";
                   outputIndividualImpactParametersFile<<_IndividualImpactSlopeByInfec[subject][counter][ser]<<" ";
                 }
                 counter++;
             }
           }
            for (i=counter;i<_nSer;i++){
                outputIndividualImpactParametersFile<<-101<<" ";
                outputIndividualImpactParametersFile<<-101<<" ";
                outputIndividualImpactParametersFile<<-101<<" ";
                for(ser=0;ser<_nSerAssay;ser++){
                    outputIndividualImpactParametersFile<<-101<<" "<<-101<<" "<<-101<<" ";
                    }
            }
    }
    outputIndividualImpactParametersFile<<"\n"<<flush;

    //Update Individual Trajectories
      int day;
      for (subject=0; subject<_numberOfSubject; subject++)
      {   outputIndividualTiterTrajectoriesFile<<-999<<" "<<iteration<<" ";
        outputIndividualTiterTrajectoriesFile<<-888<<" "<<subject<<" ";
        for(ser=0;ser<_nSerAssay;ser++){
          for(day=_minDayInfected;day<_maxDay;day++){
            outputIndividualTiterTrajectoriesFile<< computeTiterIndividualAtDay(subject,day,ser)<<" ";
          }
        }
        outputIndividualTiterTrajectoriesFile<<"\n"<<flush;
      }
    }
  }
}




//=========== main ============
int main(int argc, const char * argv[])
{
  _numberOfSubject=611;
  _numberOfSubjectSim=611;
  _totalTests=4127;
  _numberMissingSymptomaticSerotype=10;
  _numberNotMissingSymptomaticSerotype=104;
  _nSer=4;
  _nSerAssay=4;
  _minDayInfected=-359;
  _maxDay=2925;
  
  _calcMeanTiterImpact=1; //Whether to use mean baseline titer for pars16-18
  _discreteTits=0;  /// Lower bounds and upper bounds so discrete is 1
  _estNaiveBaseline=0;
  _addPreviousAugmentedInfs=0;
  _temporaryTitersAdd=0;
  _vaccineParOnly=1;  //When updating parameters, only use subset for vaccine parameters
  _sympOnly=1;
  _sympOnlyEpiCurve=0;
  _simThenEstimate=0;
  _individualLevelEffectsAllByInfec=1;
  _individualLevelEffectsAllByInfecVac=1;
  _individualLevelEffectsSlowDecay=0;
  _individualLevelEffectsD=1;
  _seperateVaccineParameters=0;
  _linearDecay=0;
  _stopEstimatingPars=0;
    
  _minTimeInfecPostVaccine=10; //Minimum time since vaccineation for first infection
  _incubationPeriod=10; //time from infection to rise
  _delayRise=4; //Time from symptom onset to rise
  maxTimeBetweenInfections=300;
  _nbIterRJMCMC=0;
  _nbIterRJMCMCDelay=1;
  _nbIterRJMCMCpostDelay=1000;
  _maxDelaySympCalc=10000;
  _maxDelayPreSympCalc=30;
  _maxDelayRise=5;
  _minTempTiterRise=0.05;  ///Minimum rise in titers after 365 days
  _minSlope=1e-10;
  _maxErrorTiter=1.5; //max difference in titers between observed and fitted when fitting model
  
  int pas=1;
  int pasParameterPrinting=1;
  int pasPrinting=1000;
  int numberOfIteration=2000000;
  _iterationPerSim=1000000;
  
  _numberOfParameters=70;
  _parameter=new double[_numberOfParameters];
  _parameter[0]=4;//Temporary Rise (R)
  _parameter[1]=0;   //Delay (D) of rise
  _parameter[2]=0.02;//0.02;//Slope of decay (S)
  _parameter[3]=0;   //Long-term decay of titers
  _parameter[4]=0;//T impact if infecting
  _parameter[5]=0;  //T impact infecting and non-primary
  _parameter[6]=0; //Value under limit of detection
   _parameter[7]=1;  //S for non-primary infections
  _parameter[8]=0.3;//Measurement error
  _parameter[9]=1;    //R for non-primary infections
  _parameter[10]=0;  //NA
  _parameter[11]=1;//Permanent Rise (T)
    _parameter[12]=0.5;//FOI
    _parameter[13]=0.1;//Individual level effect (T)
  _parameter[14]=2;//Individual level effect (R)
  _parameter[15]=2e-03;//Individual level effect (S)
  _parameter[16]=0;//DENV2
  _parameter[17]=0;//DENV3
  _parameter[18]=0;//DENV4
  _parameter[19]=0.01;   // Individual level variance slow decay - Variance
  _parameter[20]=1;   //NA
  _parameter[21]=1.21519e-01;//Y0
  _parameter[22]=1.0;   //Y1
  _parameter[23]=1.00375e+00;//Y2
  _parameter[24]=1.08185e+00;//Y3
  _parameter[25]=1.38062e+00;//Y4
  _parameter[26]=1.12476e+00;//Y5
  _parameter[27]=8.25141e-01;   //Y6
  _parameter[28]=8.52732e-01;   //Y7
  _parameter[29]=1.69038e+00;   //Y8
  _parameter[30]=1.80913e+00;   //Y9
  _parameter[31]=0;//0.1; //1.005925e-01;//Amplitude
  _parameter[32]=0;//6.43387e-01;//Delay
  _parameter[33]=1;   //NA
  _parameter[34]=1e-10;   //Vaccine T - 1
  _parameter[35]=1e-10;   //Vaccine T - 2
  _parameter[36]=2;   //Vaccine R - 1
  _parameter[37]=1;   //Vaccine R - 2
  _parameter[38]=0.005;   //Vaccine S - 1
  _parameter[39]=1;   //Vaccine S - 2
  _parameter[40]=1;   //NA
  _parameter[41]=1;   //NA
  _parameter[42]=1;   //NA
  _parameter[43]=1;   //NA
  _parameter[44]=1;   //Vaccine T - 1 Seronaive==1
  _parameter[45]=1;   //Vaccine T - 2
  _parameter[46]=4;   //Vaccine R - 1
  _parameter[47]=1;   //Vaccine R - 2
  _parameter[48]=0.01;   //Vaccine S - 1
  _parameter[49]=1;   //Vaccine S - 2
    _parameter[50]=1;   //NA
    _parameter[51]=1;   //NA
    _parameter[52]=1;   //NA
    _parameter[53]=1;   //NA
    _parameter[54]=1;   //Natural Infection Vaccineees T - 1
    _parameter[55]=_parameter[13];   //Natural Infection Vaccineees T - 2
    _parameter[56]=4;   //Natural Infection Vaccineees R - 1
    _parameter[57]=_parameter[14];   //Natural Infection Vaccineees R - 2
    _parameter[58]=0.01;   //Natural Infection Vaccineees S - 1
    _parameter[59]=_parameter[15];   //Natural Infection Vaccineees S - 2
    _parameter[60]=1;   //NA
    _parameter[61]=1;   //NA
    _parameter[62]=1;   //NA
    _parameter[63]=1;   //NA
    _parameter[64]=1;   //Natural Infection Seropositive Vaccineees T - 1
    _parameter[65]=_parameter[13];   //Natural Seropositive Infection Vaccineees T - 2
    _parameter[66]=3;   //Natural Infection Seropositive Vaccineees R - 1
    _parameter[67]=_parameter[14];   //Natural Infection Seropositive Vaccineees R - 2
    _parameter[68]=0.01;   //Natural Infection Seropositive Vaccineees S - 1
    _parameter[69]=_parameter[15];   //Natural Infection Seropositive Vaccineees S - 2

    
  _rateForRandomWalk=new double[_numberOfParameters];
  _rateForRandomWalk[0]=0.8;
  _rateForRandomWalk[1]=0.01;
  _rateForRandomWalk[2]=0.2;
  _rateForRandomWalk[3]=0.05;
  _rateForRandomWalk[4]=0.2;
  _rateForRandomWalk[5]=0.5;
  _rateForRandomWalk[6]=0.01;
  _rateForRandomWalk[7]=0.01;
  _rateForRandomWalk[8]=0.1;
  _rateForRandomWalk[9]=0.005;
  _rateForRandomWalk[10]=0.05;
  _rateForRandomWalk[11]=0.2;
  _rateForRandomWalk[12]=0.2;
  _rateForRandomWalk[13]=0.4;
  _rateForRandomWalk[14]=0.4;
  _rateForRandomWalk[15]=0.15;
  _rateForRandomWalk[16]=0.07;
  _rateForRandomWalk[17]=0.07;
  _rateForRandomWalk[18]=0.07;
  _rateForRandomWalk[19]=0.15;
  _rateForRandomWalk[20]=0.2;
  _rateForRandomWalk[21]=1.5;
  _rateForRandomWalk[22]=0.5;
  _rateForRandomWalk[23]=0.5;
  _rateForRandomWalk[24]=0.5;
  _rateForRandomWalk[25]=0.5;
  _rateForRandomWalk[26]=0.5;
  _rateForRandomWalk[27]=0.5;
  _rateForRandomWalk[28]=0.7;
  _rateForRandomWalk[29]=0.5;
  _rateForRandomWalk[30]=5;
  _rateForRandomWalk[31]=1.5;
  _rateForRandomWalk[32]=0.4;
  _rateForRandomWalk[34]=0.1;
  _rateForRandomWalk[35]=0.3;
  _rateForRandomWalk[36]=0.01;
  _rateForRandomWalk[37]=0.5;
  _rateForRandomWalk[38]=0.1;
  _rateForRandomWalk[39]=0.1;
    _rateForRandomWalk[40]=0.1;
    _rateForRandomWalk[41]=0.1;
    _rateForRandomWalk[42]=0.1;
    _rateForRandomWalk[43]=0.9;
    _rateForRandomWalk[44]=0.5;
    _rateForRandomWalk[45]=0.5;
    _rateForRandomWalk[46]=0.1;
    _rateForRandomWalk[47]=0.4;
    _rateForRandomWalk[48]=0.2;
    _rateForRandomWalk[49]=0.5;
    _rateForRandomWalk[50]=0.1;
    _rateForRandomWalk[51]=0.1;
    _rateForRandomWalk[52]=0.1;
    _rateForRandomWalk[53]=2;
    _rateForRandomWalk[54]=0.5;
    _rateForRandomWalk[55]=1;
    _rateForRandomWalk[56]=0.05;
    _rateForRandomWalk[57]=1;
    _rateForRandomWalk[58]=0.5;
    _rateForRandomWalk[59]=0.5;
    _rateForRandomWalk[60]=0.1;
    _rateForRandomWalk[61]=0.1;
    _rateForRandomWalk[62]=0.1;
    _rateForRandomWalk[63]=0.9;
    _rateForRandomWalk[64]=0.3;
    _rateForRandomWalk[65]=1;
    _rateForRandomWalk[66]=0.05;
    _rateForRandomWalk[67]=0.8;
    _rateForRandomWalk[68]=0.4;
    _rateForRandomWalk[69]=0.5;
    
  
  _selectedParameter=new int[_numberOfParameters];
  int parameterNumber;
  for(parameterNumber=0;parameterNumber<_numberOfParameters;parameterNumber++)
  {_selectedParameter[parameterNumber]=1;}
  
  if(_individualLevelEffectsD==1){
    _parameter[1]=0;
    _selectedParameter[1]=0;
  }
  //    // Sometimes off
//      _selectedParameter[0]=0;
//    _selectedParameter[2]=0;
//      _selectedParameter[8]=0;
//      _selectedParameter[11]=0;
//      _selectedParameter[12]=0;
//      _selectedParameter[13]=0;
//      _selectedParameter[14]=0;
//      _selectedParameter[15]=0;
//
//      _selectedParameter[21]=0;
//      _selectedParameter[23]=0;
//      _selectedParameter[24]=0;
//      _selectedParameter[25]=0;
//      _selectedParameter[26]=0;
//      _selectedParameter[27]=0;
//      _selectedParameter[28]=0;
//      _selectedParameter[29]=0;
//      _selectedParameter[30]=0;
//
      _selectedParameter[31]=0;
      _selectedParameter[32]=0;
//
    _selectedParameter[34]=0;
    _selectedParameter[35]=0;
//    _selectedParameter[36]=0;
        _selectedParameter[37]=0;
//    _selectedParameter[38]=0;
        _selectedParameter[39]=0;
//    _selectedParameter[44]=0;
    _selectedParameter[45]=0;
//    _selectedParameter[46]=0;
        _selectedParameter[47]=0;
//    _selectedParameter[48]=0;
        _selectedParameter[49]=0;
//
//    _selectedParameter[54]=0;
//    _selectedParameter[56]=0;
    _selectedParameter[58]=0;
//    _selectedParameter[64]=0;
//    _selectedParameter[66]=0;
    _selectedParameter[68]=0;
    
    
  if(_seperateVaccineParameters==0){
    _selectedParameter[55]=0;
    _selectedParameter[57]=0;
    _selectedParameter[59]=0;
    _selectedParameter[65]=0;
    _selectedParameter[67]=0;
    _selectedParameter[69]=0;
  }
  //
  if(_individualLevelEffectsT==0&_individualLevelEffectsTByInfec==0&_individualLevelEffectsAllByInfec==0){_selectedParameter[13]=0;_selectedParameter[55]=0;_selectedParameter[65]=0;}
  if(_individualLevelEffectsR==0&_individualLevelEffectsRByInfec==0&_individualLevelEffectsAllByInfec==0){_selectedParameter[14]=0;_selectedParameter[57]=0;_selectedParameter[67]=0;}
    if(_individualLevelEffectsS==0&_individualLevelEffectsSByInfec==0&_individualLevelEffectsAllByInfec==0){_selectedParameter[15]=0;_selectedParameter[59]=0;_selectedParameter[69]=0;}
    if(_individualLevelEffectsAllByInfecVac==0){_selectedParameter[35]=0;_selectedParameter[37]=0;_selectedParameter[39]=0;_selectedParameter[45]=0;_selectedParameter[47]=0;_selectedParameter[49]=0;}
  if(_calcMeanTiterImpact==1){
      _parameter[16]=0;
      _parameter[17]=0;
      _parameter[18]=0;
      _parameter[4]=0;
    _selectedParameter[4]=0;
      _selectedParameter[16]=0;
    _selectedParameter[17]=0;
    _selectedParameter[18]=0;
//      _discreteTits=0;
      _nSerAssay=1;
  }
    if(_individualLevelEffectsSlowDecay==0){
//        _selectedParameter[3]=0;
        _selectedParameter[19]=0;
    }
  
  //Always Off
  _selectedParameter[1]=0;
  _selectedParameter[3]=0;
  _selectedParameter[4]=0;
    _selectedParameter[5]=0;
  _selectedParameter[6]=0;
  _selectedParameter[7]=0;
  _selectedParameter[9]=0;
  _selectedParameter[10]=0;
    _selectedParameter[16]=0;
    _selectedParameter[17]=0;
    _selectedParameter[18]=0;
     _selectedParameter[19]=0;
  _selectedParameter[20]=0;
  _selectedParameter[22]=0;
  _selectedParameter[33]=0;
  _selectedParameter[40]=0;
  _selectedParameter[41]=0;
  _selectedParameter[42]=0;
  _selectedParameter[43]=0;
    _selectedParameter[50]=0;
    _selectedParameter[51]=0;
    _selectedParameter[52]=0;
    _selectedParameter[53]=0;
    _selectedParameter[60]=0;
    _selectedParameter[61]=0;
    _selectedParameter[62]=0;
    _selectedParameter[63]=0;
    _selectedParameter[35]=0;
    _selectedParameter[37]=0;
    _selectedParameter[39]=0;
    _selectedParameter[45]=0;
    _selectedParameter[47]=0;
    _selectedParameter[49]=0;
    _selectedParameter[55]=0;
    _selectedParameter[57]=0;
    _selectedParameter[59]=0;
    _selectedParameter[65]=0;
    _selectedParameter[67]=0;
    _selectedParameter[69]=0;
  
  
  _numberOfSelectedParameter=0;
  for(parameterNumber=0;parameterNumber<_numberOfParameters;parameterNumber++)
  {if(_selectedParameter[parameterNumber]==1) _numberOfSelectedParameter++;}
  
  _idOfSelectedParameter=new int[_numberOfSelectedParameter];
  int currentSelectedId=0;
  for(parameterNumber=0;parameterNumber<_numberOfParameters;parameterNumber++)
  {   if(_selectedParameter[parameterNumber]==1)
  {_idOfSelectedParameter[currentSelectedId]=parameterNumber;
    currentSelectedId++;
  }
  }
  
  _initialParameter=new double[_numberOfParameters];
  for(parameterNumber=0;parameterNumber<_numberOfParameters;parameterNumber++)
  {_initialParameter[parameterNumber]=_parameter[parameterNumber];}
  
  //=========== Paths to data================
  char *fileID;
  if(_simThenEstimate==1){
    fileID="ModelOutput_SIM1";}else{ ///CHANGE THIS FOR NEW NAMES
      fileID="ModelOutput_DataTest";}
    
// Work computer
    char* pathData="/Users/henriksalje/Documents/R/DengueTiters/CPCCohort/Data/cppInput/";
     char* pathOutput="/Users/henriksalje/Documents/RThisComputer/NMCTiterProject/";

    std::string chemin(pathData);
    std::string cheminOutput(pathOutput);
    std::string dataFile("");
    std::string dataTestingFile("");
    std::string dataCaseSerKnownFile("");
    std::string dataCaseSerUnknownFile("");
    std::string dataEpiPDFFile("");
    std::string dataPrevAugInfsFile("");
    std::string dataPrevIndParsFile("");
    
    dataFile+=chemin+"SubjectData_v2.txt";
    dataTestingFile+=chemin+"TiterData_v2.txt";
    dataCaseSerKnownFile+=chemin+"CaseDataSerKnown_v2.txt";
    dataCaseSerUnknownFile+=chemin+"CaseDataSerUnknown_v2.txt";
    //    dataEpiPDFFile+=chemin+"finEpiPDF.txt";
    dataPrevAugInfsFile+=chemin+"BringInPreviousInfecs.txt";
    //    dataPrevIndParsFile+=chemin+"Data76_IndividualParameters.txt";
    
    
    //===============data================
    _vaccineSubject=new int[_numberOfSubject];
    _numberVaccines=new int[_numberOfSubject];
    _maxDayStudy=new int[_numberOfSubject];
    _minDayStudy=new int[_numberOfSubject];
    _numberOfTests=new int[_numberOfSubject];
    _numberOfHITests=new int[_numberOfSubject];
    _serostatusSubject=new int[_numberOfSubject];
    _dayOfTesting=new int*[_numberOfSubject];
    _dayVacc=new int*[_numberOfSubject];
    _dayOfHITesting=new int*[_numberOfSubject];
    _minDelaySinceInfection=new double*[_numberOfSubject];
    _minTimPreInfection=new double*[_numberOfSubject];
    _resultOfTesting=new double**[_numberOfSubject];
    _resultOfHITesting=new double**[_numberOfSubject];
    _daySymptomaticSerKnown=new int*[_numberOfSubject];
    _daySymptomaticSerUnknown=new int*[_numberOfSubject];
    _serotypeSymptomatic=new int*[_numberOfSubject];
    _freqSymptomaticSerKnown=new int[_numberOfSubject];
    _freqSymptomaticSerUnknown=new int[_numberOfSubject];
    int i,j;
    for (i=0;i<_numberOfSubject;i++){
      _dayVacc[i]=new int[1];
    }
    
    
    //========== augmented data ==========
    _numberAsymptomatic=new int[_numberOfSubject];
    _serotypeAsymptomatic=new int*[_numberOfSubject];
    _dayAsymptomatic=new int*[_numberOfSubject];
    _augmentedSerotypeSymptomatic=new int*[_numberOfSubject];
    _idsMissingSymptomaticSerotype=new int[_numberMissingSymptomaticSerotype];
    _ranksMissingSymptomaticSerotype=new int[_numberMissingSymptomaticSerotype];
    _baselineTiter=new double*[_numberOfSubject];
    _IndividualImpactSlowDecay=new double*[_numberOfSubject];
    _baselineHITiter=new double*[_numberOfSubject];
    _baselineNaive=new int[_numberOfSubject];
    _IndividualImpactTiter=new double[_numberOfSubject];
    _IndividualImpactTiterByInfec=new double**[_numberOfSubject];
    _IndividualImpactTempTiterByInfec=new double**[_numberOfSubject];
    _IndividualImpactSlopeByInfec=new double**[_numberOfSubject];
    _IndividualImpactTiterByInfecVac=new double**[_numberOfSubject];
    _IndividualImpactTempTiterByInfecVac=new double**[_numberOfSubject];
    _IndividualImpactSlopeByInfecVac=new double**[_numberOfSubject];
    _IndividualImpactTempTiter=new double[_numberOfSubject];
    _IndividualImpactSlope=new double[_numberOfSubject];
    _IndividualImpactDelay= new double[_numberOfSubject];
    _normFactor= new double[_numberOfSubject];
    _nAsymptomatic=0;
    for (i=0;i<_numberOfSubject;i++){
      _baselineTiter[i]=new double[_nSerAssay];
      _IndividualImpactSlowDecay[i]=new double[_nSerAssay];
      _baselineHITiter[i]=new double[_nSerAssay];
    }
    for (i=0;i<_numberOfSubject;i++){
      _IndividualImpactTiterByInfec[i]=new double*[_nSer];
      _IndividualImpactTempTiterByInfec[i]=new double*[_nSer];
      _IndividualImpactSlopeByInfec[i]=new double*[_nSer];
      _IndividualImpactTiterByInfecVac[i]=new double*[_nSer];
      _IndividualImpactTempTiterByInfecVac[i]=new double*[_nSer];
      _IndividualImpactSlopeByInfecVac[i]=new double*[_nSer];
      _dayAsymptomatic[i]=new int[_nSer];
      _serotypeAsymptomatic[i]= new int[_nSer];
      for(j=0;j<_nSer;j++){
        _IndividualImpactTiterByInfec[i][j]=new double[_nSerAssay];
        _IndividualImpactTempTiterByInfec[i][j]=new double[_nSerAssay];
        _IndividualImpactSlopeByInfec[i][j]=new double[_nSerAssay];
        _IndividualImpactTiterByInfecVac[i][j]=new double[_nSerAssay];
        _IndividualImpactTempTiterByInfecVac[i][j]=new double[_nSerAssay];
        _IndividualImpactSlopeByInfecVac[i][j]=new double[_nSerAssay];
      }
    }
    
    //========== parameter-linked ==========
    _meanTiterindividual=new double**[_numberOfSubject];
    _probaTiterIndividual=new double[_numberOfSubject];
    _probabInfecHistoryIndividual=new double*[_numberOfSubject];
    for (i=0;i<_numberOfSubject;i++)
    {   _probabInfecHistoryIndividual[i]=new double[_nSer];
      _IndividualImpactTiter[i]=1;
      _IndividualImpactTempTiter[i]=1;
      _IndividualImpactSlope[i]=1;
      _IndividualImpactDelay[i]=0;
      _normFactor[i]=1;
    }
    _tmpTvec=new double[_nSerAssay];
    _tmpRvec=new double[_nSerAssay];
    _tmpSvec=new double[_nSerAssay];
    _tmpTvecVac=new double[_nSerAssay];
    _tmpRvecVac=new double[_nSerAssay];
    _tmpSvecVac=new double[_nSerAssay];
    _oldBaselines=new double[_nSerAssay];
    _oldSlowDecay=new double[_nSerAssay];
    _oldTiters=new double*[_nSer];
    _oldTempTiters=new double*[_nSer];
    _oldSlopes=new double*[_nSer];
    _oldTitersVac=new double*[_nSer];
    _oldTempTitersVac=new double*[_nSer];
    _oldSlopesVac=new double*[_nSer];
    _newIndEffectsT=new double*[_nSer];
    _newIndEffectsR=new double*[_nSer];
    _newIndEffectsS=new double*[_nSer];
//    _newIndividualImpact=new double[_nSerAssay];
//    _newIndividualImpactR=new double[_nSerAssay];
//    _newIndividualImpactS=new double[_nSerAssay];
    for (i=0;i<_nSer;i++){
      _oldTiters[i]=new double[_nSerAssay];
      _oldTempTiters[i]=new double[_nSerAssay];
      _oldSlopes[i]=new double[_nSerAssay];
        _oldTitersVac[i]=new double[_nSerAssay];
        _oldTempTitersVac[i]=new double[_nSerAssay];
        _oldSlopesVac[i]=new double[_nSerAssay];
      _newIndEffectsT[i]=new double[_nSerAssay];
      _newIndEffectsR[i]=new double[_nSerAssay];
      _newIndEffectsS[i]=new double[_nSerAssay];
    }
    //========== outputs ==========
    std::string outputParametersFile("");
    outputParametersFile+=cheminOutput+"outputParameters_"+fileID+".txt";
    
    std::string outputIndividualTiterSimFile("");
    outputIndividualTiterSimFile+=cheminOutput+"outputIndividualTiterSim_"+fileID+".txt";
    
    std::string outputIndividualImpactParametersFile("");
    outputIndividualImpactParametersFile+=cheminOutput+"outputIndividualImpactParameters_"+fileID+".txt";
    
    std::string outputIndividualInfectionHistoriesFile("");
    outputIndividualInfectionHistoriesFile+=cheminOutput+"outputIndividualInfectionHistories_"+fileID+".txt";
    
    std::string outputIndividualTiterTrajectoriesFile("");
    outputIndividualTiterTrajectoriesFile+=cheminOutput+"outputIndividualTiterTrajectories_"+fileID+".txt";
    
    
    
    //=========== build it and run it ==========
    buildSubjectData(dataFile.c_str());
    buildTestingData(dataTestingFile.c_str());
    buildSymptomaticSerKnownData(dataCaseSerKnownFile.c_str());
    buildSymptomaticSerUnknownData(dataCaseSerUnknownFile.c_str());
    updateEpiPDF();
    
    double meanR, par1,par2, indR,dd;
    int ser,  totInfec, infecNo, minRise;
    int check, subject;
    
    for (subject=0;subject<_numberOfSubject;subject++){
        minRise=1;
        while(minRise==1){
            for (ser=0;ser<_nSerAssay;ser++){
                _baselineTiter[subject][ser]=drawBaselineTiter(_baselineNaive[subject]);
                if(_individualLevelEffectsSlowDecay==1){
                    par1=_parameter[3]*_parameter[3]/_parameter[19];
                    par2=_parameter[3]/_parameter[19];
                    _IndividualImpactSlowDecay[subject][ser]=max(rgamma(par1,par2),1e-20);
                }
            }
            _numberAsymptomatic[subject]=0;
            totInfec=_freqSymptomaticSerKnown[subject]+_freqSymptomaticSerUnknown[subject];
            if(_vaccineSubject[subject]==0){
                for (infecNo=0;infecNo<totInfec;infecNo++){
                    for (ser=0;ser<_nSerAssay;ser++){
                        double mean=_parameter[11];//_parameter[11]*_parameter[54];
                        par1=mean*mean/_parameter[13];
                        par2=mean/_parameter[13];
                        _IndividualImpactTiterByInfec[subject][infecNo][ser]=max(1e-10,rgamma(par1,par2));
                    
                        mean=_parameter[0];//_parameter[0]*_parameter[56];
                        par1=mean*mean/_parameter[14];
                        par2=mean/_parameter[14];
                        _IndividualImpactTempTiterByInfec[subject][infecNo][ser]=max(1e-10,rgamma(par1,par2));
                        
                        mean=_parameter[2];//_parameter[2]*_parameter[58];
                        par1=mean*mean/_parameter[15];
                        par2=mean/_parameter[15];
                        _IndividualImpactSlopeByInfec[subject][infecNo][ser]=max(1e-10,rgamma(par1,par2));
                    }
                }
            }
            if(_vaccineSubject[subject]==1&_baselineNaive[subject]==1){
                for (infecNo=0;infecNo<totInfec;infecNo++){
                    for (ser=0;ser<_nSerAssay;ser++){
                        double mean=_parameter[54];//*_parameter[11];
                        par1=mean*mean/_parameter[55];
                        par2=mean/_parameter[55];
                        _IndividualImpactTiterByInfec[subject][infecNo][ser]=max(1e-10,rgamma(par1,par2));
                        
                        mean=_parameter[56];//*_parameter[0];
                        par1=mean*mean/_parameter[57];
                        par2=mean/_parameter[57];
                        _IndividualImpactTempTiterByInfec[subject][infecNo][ser]=max(1e-10,rgamma(par1,par2));
                        
                        mean=_parameter[58];//*_parameter[2];
                        par1=mean*mean/_parameter[59];
                        par2=mean/_parameter[59];
                        _IndividualImpactSlopeByInfec[subject][infecNo][ser]=max(1e-10,rgamma(par1,par2));
                    }
                }
            }
            if(_vaccineSubject[subject]==1&_baselineNaive[subject]==0){
                for (infecNo=0;infecNo<totInfec;infecNo++){
                    for (ser=0;ser<_nSerAssay;ser++){
                        double mean=_parameter[64];//*_parameter[11];
                        par1=mean*mean/_parameter[65];
                        par2=mean/_parameter[65];
                        _IndividualImpactTiterByInfec[subject][infecNo][ser]=max(1e-10,rgamma(par1,par2));
                        
                        mean=_parameter[66];//*_parameter[0];
                        par1=mean*mean/_parameter[67];
                        par2=mean/_parameter[67];
                        _IndividualImpactTempTiterByInfec[subject][infecNo][ser]=max(1e-10,rgamma(par1,par2));
                        
                        mean=_parameter[68];//*_parameter[2];
                        par1=mean*mean/_parameter[69];
                        par2=mean/_parameter[69];
                        _IndividualImpactSlopeByInfec[subject][infecNo][ser]=max(1e-10,rgamma(par1,par2));
                    }
                }
            }
            if(_vaccineSubject[subject]==1&_individualLevelEffectsAllByInfecVac==1&_baselineNaive[subject]==0){
              for (i=0;i<_numberVaccines[subject];i++){
                for (ser=0;ser<_nSerAssay;ser++){
                  double par1=_parameter[34]*_parameter[34]/_parameter[35];
                  double par2=_parameter[34]/_parameter[35];
                  _IndividualImpactTiterByInfecVac[subject][i][ser]=max(1e-10,rgamma(par1,par2));
                  
                  par1=_parameter[36]*_parameter[36]/_parameter[14];
                  par2=_parameter[36]/_parameter[14];
                  _IndividualImpactTempTiterByInfecVac[subject][i][ser]=max(1e-10,rgamma(par1,par2));
                  
                  par1=_parameter[38]*_parameter[38]/_parameter[15];
                  par2=_parameter[38]/_parameter[15];
                  _IndividualImpactSlopeByInfecVac[subject][i][ser]=max(1e-10,rgamma(par1,par2));
                }
              }
            }
            if(_vaccineSubject[subject]==1&_individualLevelEffectsAllByInfecVac==1&_baselineNaive[subject]==1){
                for (i=0;i<_numberVaccines[subject];i++){
                    for (ser=0;ser<_nSerAssay;ser++){
                        double par1=_parameter[44]*_parameter[44]/_parameter[13];
                        double par2=_parameter[44]/_parameter[13];
                        _IndividualImpactTiterByInfecVac[subject][i][ser]=max(1e-10,rgamma(par1,par2));
                        
                        par1=_parameter[46]*_parameter[46]/_parameter[14];
                        par2=_parameter[46]/_parameter[14];
                        _IndividualImpactTempTiterByInfecVac[subject][i][ser]=max(1e-10,rgamma(par1,par2));
                        
                        par1=_parameter[48]*_parameter[48]/_parameter[15];
                        par2=_parameter[48]/_parameter[15];
                        _IndividualImpactSlopeByInfecVac[subject][i][ser]=max(1e-10,rgamma(par1,par2));
                        
                    }
                }
            }
            minRise=checkMinRise(subject,0);
        }
    }
    
    
    
    if(_addPreviousAugmentedInfs==1){
      buildPreviousAugmentedInfections(dataPrevAugInfsFile.c_str());
    }
    
    maxDelaySympInfec();
    addSerotypesAll();
    computeProbaAll(&_logLikGlobal,0,0,0);
    runMCMC(outputParametersFile.c_str(),outputIndividualTiterSimFile.c_str(),outputIndividualImpactParametersFile.c_str(),outputIndividualInfectionHistoriesFile.c_str(),outputIndividualTiterTrajectoriesFile.c_str(),pas,pasPrinting,pasParameterPrinting,numberOfIteration);
    
    cout<<"DONE!!";
    
    return(0);
    
}
