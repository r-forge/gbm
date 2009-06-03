//------------------------------------------------------------------------------
//  GBM alteration by Daniel Edwards
//
//  File:       bisquare.h
//
//
//  Contains:   Distribution object to implement bisquare
//
//  History:    04/04/2008   Created
//                    
//------------------------------------------------------------------------------

#ifndef KCLASSCGBM_H
#define KCLASSCGBM_H

#include <algorithm>
#include "distribution.h"
#include "locationm.h"

class CKClass : public CDistribution
{

public:

    CKClass(int cNumClasses, int cRows);

    virtual ~CKClass();
    GBMRESULT UpdateParams(double *adF,
	                       double *adOffset,
						   double *adWeight,
	                       unsigned long cLength);

    GBMRESULT ComputeWorkingResponse(double *adY,
                                     double *adMisc,
                                     double *adOffset,
                                     double *adF, 
                                     double *adZ, 
                                     double *adWeight,
                                     bool *afInBag,
                                     unsigned long nTrain,
	                                 int cIdxOff);

    GBMRESULT InitF(double *adY, 
                    double *adMisc,
                    double *adOffset,
                    double *adWeight,
                    double &dInitF, 
                    unsigned long cLength);

    GBMRESULT FitBestConstant(double *adY,
                              double *adMisc,
                              double *adOffset,
                              double *adW,
                              double *adF,
                              double *adZ,
                              unsigned long *aiNodeAssign,
                              unsigned long nTrain,
                              VEC_P_NODETERMINAL vecpTermNodes,
                              unsigned long cTermNodes,
                              unsigned long cMinObsInNode,
                              bool *afInBag,
                              double *adFadj,
	                          int cIdxOff);

    double Deviance(double *adY,
                    double *adMisc,
                    double *adOffset,
                    double *adWeight,
                    double *adF,
                    unsigned long cLength,
	                int cIdxOff);

    double BagImprovement(double *adY,
                          double *adMisc,
                          double *adOffset,
                          double *adWeight,
                          double *adF,
                          double *adFadj,
                          bool *afInBag,
                          double dStepSize,
                          unsigned long nTrain);

private:
	int mcNumClasses;
	int mcRows;
	double *madProb; 
};

#endif // KCLASSCGBM_H



