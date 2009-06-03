//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "bisquare.h"

CBisquare::CBisquare(double adR)
{
    mdR = adR;

	double *adParams = new double[1];
	adParams[0] = adR;

	mpLocM = new CLocationM("bisquare", 1, adParams);

	delete[] adParams;
}

CBisquare::~CBisquare()
{
	delete mpLocM;
}


GBMRESULT CBisquare::ComputeWorkingResponse
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
	double dU = 0.0;

    if(adOffset == NULL)
    {
        for(i=0; i<nTrain; i++)
        {
			adZ[i] = 0;
		    dU = adY[i] - adF[i];	
			if (fabs(dU) < mdR){
				adZ[i] = dU * (1 - ((dU * dU) / (mdR * mdR))) * (1 - ((dU * dU) / (mdR * mdR)));
			}
            
        }
    }
    else
    {
        for(i=0; i<nTrain; i++)
        {
 			adZ[i] = 0;
		    dU = adY[i] - adOffset[i] - adF[i];	
			if (dU < mdR){
				adZ[i] = dU * (1 - ((dU * dU) / (mdR * mdR))) * (1 - ((dU * dU) / (mdR * mdR)));
			}
        }
    }

    return GBM_OK;
}


GBMRESULT CBisquare::InitF
(
    double *adY,
    double *adMisc,
    double *adOffset, 
    double *adWeight,
    double &dInitF, 
    unsigned long cLength
)
{

	// Local variables
	int ii;

	// Get objects to pass into the LocM function
	int iN = int(cLength);
	double *adArr = new double[iN];

	for (ii = 0; ii < iN; ii++)
	{
		double dOffset = (adOffset==NULL) ? 0.0 : adOffset[ii];
		adArr[ii] = adY[ii] - dOffset;
	}

	dInitF = mpLocM->LocationM(iN, adArr, adWeight);

    delete[] adArr;

    return GBM_OK;
}

double CBisquare::Deviance
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
    double dW = 0.0;
	double dU = 0.0;

    if(adOffset == NULL)
    {
        for(i=cIdxOff; i<cLength+cIdxOff; i++)
        {
			dU = adY[i] - adF[i];
			if (fabs(dU < mdR))
			{
				dU = (1 - ((dU * dU) / (mdR * mdR)));
				dL += adWeight[i] * (mdR * mdR / 6) * (1 -  dU * dU * dU);
			}

            dW += adWeight[i];
        }
    }
    else
    {
        for(i=cIdxOff; i<cLength+cIdxOff; i++)
        {
			dU = adY[i] - adOffset[i] - adF[i];
			if (fabs(dU < mdR))
			{
				dU = (1 - ((dU * dU) / (mdR * mdR)));
				dL += adWeight[i] * (mdR * mdR / 6) * (1 -  dU * dU * dU);
			}

            dW += adWeight[i];
        }
    }

    return dL/dW;
}


GBMRESULT CBisquare::FitBestConstant
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
   	// Local variables
    GBMRESULT hr = GBM_OK;
	unsigned long iNode = 0;
    unsigned long iObs = 0;
    double dOffset;

	// Call LocM for the array of values on each node
    for(iNode=0; iNode<cTermNodes; iNode++)
    {
        if(vecpTermNodes[iNode]->cN >= cMinObsInNode)
        {
			// Get the number of nodes here
			int iNumNodes = 0;
			for (iObs = 0; iObs < nTrain; iObs++)
			{
				if(afInBag[iObs] && (aiNodeAssign[iObs] == iNode))
                {
                    iNumNodes++;
                }
			}

			// Create the arrays to centre
			double *adArr = new double[iNumNodes];
			double *adWeight = new double[iNumNodes];

			int iIdx = 0;
			for(iObs=0; iObs<nTrain; iObs++)
            {
                if(afInBag[iObs] && (aiNodeAssign[iObs] == iNode))
                {
                    dOffset = (adOffset==NULL) ? 0.0 : adOffset[iObs];
                    adArr[iIdx] = adY[iObs] - dOffset - adF[iObs];
					adWeight[iIdx] = adW[iObs];

					iIdx++;
                }
            }

           	vecpTermNodes[iNode]->dPrediction = mpLocM->LocationM(iNumNodes, adArr, 
				                                                 adWeight); 

			delete[] adArr;
			delete[] adWeight;

        }
    }

    return hr;
}

double CBisquare::BagImprovement
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
	double dU = 0.0;
    double dV = 0.0;

    for(i=0; i<nTrain; i++)
    {
        if(!afInBag[i])
        {
            dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);
            
			dU = (adY[i] - dF) / mdR;
			dU = 1 - (dU * dU);

			dV = (adY[i] - dF - dStepSize * adFadj[i]) / mdR;
			dV = 1 - (dV * dV);

            dReturnValue += 
                adWeight[i] * (mdR * mdR / 6) * ((dV * dV * dV) - (dU * dU * dU));
            dW += adWeight[i];
        }
    }

    return dReturnValue/dW;
}




