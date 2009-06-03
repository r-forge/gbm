// GBM by Greg Ridgeway  Copyright (C) 2003

#include "huberized.h"

CHuberized::CHuberized()
{
}

CHuberized::~CHuberized()
{
}


GBMRESULT CHuberized::ComputeWorkingResponse
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adF, 
    double *adZ, 
    double *adWeight,
    bool *afInBag,
    unsigned long nTrain,
	int cIdxOff
)
{
    unsigned long i = 0;
    double dProb = 0.0;
    double dF = 0.0;

    for(i=0; i<nTrain; i++)
    {
        dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);
//        dProb = 1.0/(1.0+exp(-dF));
          if( (2*adY[i]-1)*dF < -1){
              dProb = -4 * (2*adY[i]-1) * dF;
              adZ[i] = -4 * (2*adY[i]-1);
          }
          else{
              dProb = ( 1 - (2*adY[i]-1)*dF )*(( 1 - (2*adY[i]-1)*dF ));
              adZ[i] = -2 * (2*adY[i]-1) + 2 * (2*adY[i]-1)*(2*adY[i]-1)*dF;
          }
//        adZ[i] = adY[i] - dProb;
    }

    return GBM_OK;
}



GBMRESULT CHuberized::InitF
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double &dInitF, 
    unsigned long cLength
)
{
    GBMRESULT hr = GBM_OK;

    unsigned long i=0;
    double dTemp=0.0;

    if(adOffset==NULL)
    {
        double dSum=0.0;
        for(i=0; i<cLength; i++)
        {
            dSum += adWeight[i]*adY[i];
            dTemp += adWeight[i];
        }
        dInitF = log(dSum/(dTemp-dSum));
    }
    else
    { 
        // Newton method for solving for F
        // should take about 3-6 iterations.
        double dNum=0.0;         // numerator
        double dDen=0.0;         // denominator
        double dNewtonStep=1.0;  // change
        dInitF = 0.0;
        while(dNewtonStep > 0.0001)
        {
            dNum=0.0;
            dDen=0.0;
            if( (2*adY[i]-1)*dInitF < -1){
            	for(i=0; i<cLength; i++)
            	{
//                dTemp = 
                  dTemp = -4 * (2*adY[i]-1) * (adOffset[i]+dInitF);
				  dNum += adWeight[i]*(adY[i]-dTemp);
			      dDen += adWeight[i] * (-4 * (2*adY[i]-1) );
				}
            }
            else
			{
               for(i=0; i < cLength; i++)
               {
               dTemp = ( 1 - (2*adY[i]-1)*(adOffset[i]+dInitF) )*(( 1 - (2*adY[i]-1)*(adOffset[i]+dInitF) ));
               dNum += adWeight[i]*(adY[i]-dTemp);
               dDen += adWeight[i]*2*dTemp*(1-(2*adY[i]-1)*dTemp);
               }
            }
            dNewtonStep = dNum/dDen;
            dInitF += dNewtonStep;
        }
    }

    return hr;
}



double CHuberized::Deviance
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double *adF,
    unsigned long cLength,
    int cIdxOff
)
{
   unsigned long i=0;
   double dL = 0.0;
   double dF = 0.0;
   double dW = 0.0;

   if(adOffset==NULL)
   {
	  if ( (2*adY[i]-1)*dF < -1 ){
          for(i=cIdxOff; i<cLength+cIdxOff; i++)
          {
             dL += adWeight[i]*((2*adY[i]-1)*(2*adY[i]-1) + 4*(2*adY[i]-1)*dF);
             dW += adWeight[i];
          }
      } // close if(adY[i]...
	  else{
         for(i=cIdxOff; i<cLength+cIdxOff; i++){
             dL += adWeight[i]*((2*adY[i]-1)*(2*adY[i]-1) -(1 - (2*adY[i]-1) * dF) * (1 - (2*adY[i]-1) * dF) );
         }
       }
   } // close if (adOffset...
   else
   {
      if ( (2*adY[i]-1)*dF < -1 ){
          for(i=cIdxOff; i<cLength+cIdxOff; i++)
          {
             dF = (2*adY[i]-1) + adOffset[i];
             dL += adWeight[i]*((2*adY[i]-1)*dF + 4*(2*adY[i]-1)*dF );
             dW += adWeight[i];
          }
       }
       else{
          for(i=cIdxOff; i<cLength+cIdxOff; i++)
             dF = (2*adY[i]-1) + adOffset[i];
             dL += adWeight[i] * ((2*adY[i]-1)*dF - (1 - (2*adY[i]-1) * dF) * (1 - (2*adY[i]-1) * dF) );
             dW += adWeight[i];
       }
   }

   return -2*dL/dW;
}


GBMRESULT CHuberized::FitBestConstant
(
    double *adY,
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
	int cIdxOff
)
{
    GBMRESULT hr = GBM_OK;

    unsigned long iObs = 0;
    unsigned long iNode = 0;
    vecdNum.resize(cTermNodes);
    vecdNum.assign(vecdNum.size(),0.0);
    vecdDen.resize(cTermNodes);
    vecdDen.assign(vecdDen.size(),0.0);

    for(iObs=0; iObs<nTrain; iObs++)
    {
        if(afInBag[iObs])
        {
            vecdNum[aiNodeAssign[iObs]] += adW[iObs]*adZ[iObs];
            vecdDen[aiNodeAssign[iObs]] += 
                adW[iObs]*(adY[iObs]-adZ[iObs])*(1-adY[iObs]+adZ[iObs]);
        }
    }

    for(iNode=0; iNode<cTermNodes; iNode++)
    {
        if(vecpTermNodes[iNode]!=NULL)
        {
            if(vecdDen[iNode] == 0)
            {
                vecpTermNodes[iNode]->dPrediction = 0.0;
            }
            else
            {
                vecpTermNodes[iNode]->dPrediction = 
                    vecdNum[iNode]/vecdDen[iNode];
            }
        }
    }

    return hr;
}


double CHuberized::BagImprovement
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double *adF,
    double *adFadj,
    bool *afInBag,
    double dStepSize,
    unsigned long nTrain
)
{
    double dReturnValue = 0.0;
    double dF = 0.0;
    double dW = 0.0;
    unsigned long i = 0;

    for(i=0; i<nTrain; i++)
    {
        if(!afInBag[i])
        {
            dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);

            if(adY[i]==1.0)
            {
                dReturnValue += adWeight[i]*dStepSize*adFadj[i];
            }
            if ( (2*adY[i]-1)*dF < -1 )
            {
                dReturnValue += adWeight[i]*
                                (-4*(2*adY[i]-1)*dF - 
                                 -4*(2*adY[i]-1)*(dF+dStepSize*adFadj[i]));
                dW += adWeight[i];
            }
            else{
                dReturnValue += adWeight[i]*
                                ( (1 - (2*adY[i]-1)*dF )*(1 - (2*adY[i]-1)*dF )*dF -
                                  (1 - (2*adY[i]-1)*dF+dStepSize*adFadj[i] )*(1 - (2*adY[i]-1)*dF+dStepSize*adFadj[i] )*dF+dStepSize*adFadj[i]
                                );
            }
        }
    }

    return dReturnValue/dW;
}


