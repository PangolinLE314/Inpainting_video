
//this function defines the patch match measure with which we compare patches

#include "patch_match_measure.h"

#define LAMDA_MEASURE 0

template <class T>
float ssd_patch_measure(nTupleVolume<T> *imgVolA, nTupleVolume<T> *imgVolB, nTupleVolume<T> *dispField, nTupleVolume<T> *occVol, int xA, int yA, int tA,
                        int xB, int yB, int tB, float minVal, const parameterStruct *params)
{
    //declarations
    int i,j,k, p,xAtemp, yAtemp, tAtemp, xBtemp, yBtemp, tBtemp;
    int xMinA,yMinA,tMinA;
    int xMinB,yMinB,tMinB;
    int sumOcc, occA;
    float tempFloat;
    float beta = 50.0;
    float *ptrA, *ptrB;

    float ssd = 0;

    if ( ((imgVolA->patchSizeX) != (imgVolB->patchSizeX)) || ((imgVolA->patchSizeY) != (imgVolB->patchSizeY)) ||
        ((imgVolA->patchSizeT) != (imgVolB->patchSizeT)) )
    {
        MY_PRINTF("Error in ssd_minimum_value, the patch sizes are not equal.\n");
        return -1;
    }

    if (params->patchIndexing == 0)
    {
        xMinA = xA-imgVolA->hPatchSizeX;
        yMinA = yA-imgVolA->hPatchSizeY;
        tMinA = tA-imgVolA->hPatchSizeT;
        xMinB = xB-imgVolB->hPatchSizeX;
        yMinB = yB-imgVolB->hPatchSizeY;
        tMinB = tB-imgVolB->hPatchSizeT;
    }
    else if(params->patchIndexing == 1)
    {
        xMinA = xA;
        yMinA = yA;
        tMinA = tA;
        xMinB = xB;
        yMinB = yB;
        tMinB = tB;
    }

    sumOcc = 0;

    if (params->partialComparison)
    {
        for (k=0; k<imgVolA->patchSizeT; k++)
            for (j=0; j<imgVolA->patchSizeY; j++)
                for (i=0; i<imgVolA->patchSizeX; i++)
                {
                    xAtemp = xMinA + i;
                    yAtemp = yMinA + j;
                    tAtemp = tMinA + k;
                    xBtemp = xMinB + i;
                    yBtemp = yMinB + j;
                    tBtemp = tMinB + k;
                    //do not compare the edges in any case
                    if ((!check_in_boundaries( imgVolA, xAtemp, yAtemp, tAtemp,params)))
                        continue;   //if we are not in the inner boundaries, do not compare

                    sumOcc = sumOcc + (int)(!(*(occVol->get_value_ptr(xAtemp,yAtemp,tAtemp,0))) == 1);

                }
    }
    else    //calculate patch size
    {
        int patchSizeXtemp = min_int(xA + imgVolA->hPatchSizeX,imgVolA->xSize-1) - max_int(xA - imgVolA->hPatchSizeX,0) + 1;
        int patchSizeYtemp = min_int(yA + imgVolA->hPatchSizeY,imgVolA->ySize-1) - max_int(yA - imgVolA->hPatchSizeY,0) + 1;
        int patchSizeTtemp = min_int(tA + imgVolA->hPatchSizeT,imgVolA->tSize-1) - max_int(tA - imgVolA->hPatchSizeT,0) + 1;

        sumOcc = patchSizeXtemp * patchSizeYtemp * patchSizeTtemp;
    }
    sumOcc = max_int(sumOcc,1);

    for (k=0; k<imgVolA->patchSizeT; k++)
        for (j=0; j<imgVolA->patchSizeY; j++)
            for (i=0; i<imgVolA->patchSizeX; i++)
            {
                xAtemp = xMinA + i;
                yAtemp = yMinA + j;
                tAtemp = tMinA + k;

                xBtemp = xMinB + i;
                yBtemp = yMinB + j;
                tBtemp = tMinB + k;


                //do not compare if we are not in the boundaries
                if ((!check_in_boundaries( imgVolA, xAtemp, yAtemp, tAtemp,params)))
                    occA = 1;
                else
                    occA = 0;
                if (occA == 1)
                    continue;   //we do not wish to compare this pixel
                /*if we want partial patch comparison*/
                if (params->partialComparison && occVol->xSize >0)
                    occA = (int)(*(occVol->get_value_ptr(xAtemp, yAtemp, tAtemp,0)) == 1);
                if (occA == 1)
                    continue;   //we do not wish to compare this pixel

                ptrA = imgVolA->get_begin_ptr(xAtemp, yAtemp, tAtemp);
                ptrB = imgVolB->get_begin_ptr(xBtemp, yBtemp, tBtemp);

                /* similarity */ //new
              for (p=0; p<imgVolA->nTupleSize; p++)
         //       for (p=0; p<1; p++)
                {
                    tempFloat = (*(ptrA+p)) - (*(ptrB+p));
                    ssd = ssd + (((tempFloat)*(tempFloat))/sumOcc);
                    //ssd = ssd + (abs(tempFloat))/sumOcc;
                }

                if( params->gradX != NULL)
                {
                    float normGradXtemp = *((params->normGradX) + (tAtemp)*(imgVolA->xSize)*(imgVolA->ySize)
                             + (xAtemp)*(imgVolA->ySize) + yAtemp) -
                            *((params->normGradX) + (tBtemp)*(imgVolB->xSize)*(imgVolB->ySize)
                             + (xBtemp)*(imgVolB->ySize) + yBtemp);

                    float normGradYtemp = *((params->normGradY) + (tAtemp)*(imgVolA->xSize)*(imgVolA->ySize)
                             + (xAtemp)*(imgVolA->ySize) + yAtemp) -
                            *((params->normGradY) + (tBtemp)*(imgVolB->xSize)*(imgVolB->ySize)
                             + (xBtemp)*(imgVolB->ySize) + yBtemp);

                    ssd = ssd + beta*normGradXtemp*normGradXtemp/sumOcc;
                    ssd = ssd + beta*normGradYtemp*normGradYtemp/sumOcc;
                }

        if ((minVal != -1) && (ssd > minVal))
                {
            return(-1);
                }
            }

    return(ssd);
}




template <class T>
float ssd_patch_measure_propagate(nTupleVolume<T> *imgVolA, nTupleVolume<T> *imgVolB, nTupleVolume<T> *dispField, nTupleVolume<T> *occVol, int xA, int yA, int tA,
                        int xB, int yB, int tB, float minVal, const parameterStruct *params, int dirPropagate)
{
 if(params->partialComparison)
     return ssd_patch_measure(imgVolA,imgVolB,dispField,occVol,xA,yA,tA,xB,yB,tB,minVal,params);
    float ssd = 0;
    int tMax=imgVolA->patchSizeT;
    int yMax=imgVolA->patchSizeY;
    int xMax=imgVolA->patchSizeX;

    int xMin=0;
    int yMin=0;
    int tMin=0;

    switch(dirPropagate){
    case -1: //Left
            ssd = (float)dispField->get_value((int)max_int(xA-1,0),yA,tA,3);
            if (ssd==FLT_MAX)
                return ssd_patch_measure(imgVolA,imgVolB,dispField,occVol,xA,yA,tA,xB,yB,tB,minVal,params);
            else{
                xMax=0;
                xMin=-1;
            }
        break;
    case -2: // Top
        ssd = (float)dispField->get_value(xA,(int)max_int(yA-1,0),tA,3);
        if (ssd==FLT_MAX)
            return ssd_patch_measure(imgVolA,imgVolB,dispField,occVol,xA,yA,tA,xB,yB,tB,minVal,params);
        else{
            yMax=0;
            yMin=-1;
        }
        break;
    case -3:
        ssd = (float)dispField->get_value(xA,yA,(int)max_int(tA-1,0),3);
        if (ssd==FLT_MAX)
            return ssd_patch_measure(imgVolA,imgVolB,dispField,occVol,xA,yA,tA,xB,yB,tB,minVal,params);
        else{
            tMax=0;
            tMin=-1;
        }
        break;
    case 1: // RIGHT
        ssd = (float)dispField->get_value((int)min_int(xA+1,(imgVolA->xSize)-1),yA,tA,3);
        if (ssd==FLT_MAX)
            return ssd_patch_measure(imgVolA,imgVolB,dispField,occVol,xA,yA,tA,xB,yB,tB,minVal,params);
        else{
            xMax=imgVolA->patchSizeX+1;
            xMin=xMax-1;
        }
        break;
    case 2:
        ssd = (float)dispField->get_value(xA,(int)min_int(yA+1,(imgVolA->ySize)-1),tA,3);
        if (ssd==FLT_MAX)
            return ssd_patch_measure(imgVolA,imgVolB,dispField,occVol,xA,yA,tA,xB,yB,tB,minVal,params);
        else{
            yMax=imgVolA->patchSizeY+1;
            yMin=yMax-1;
        }
        break;
    case 3:
        ssd = (float)dispField->get_value(xA,yA,(int)min_int(tA+1,(imgVolA->tSize)-1),3);
        if (ssd==FLT_MAX)
            return ssd_patch_measure(imgVolA,imgVolB,dispField,occVol,xA,yA,tA,xB,yB,tB,minVal,params);
        else{
            tMax=imgVolA->patchSizeT+1;
            tMin=tMax-1;
        }
        break;
    }
    //declarations
    int i,j,k, p,xAtemp, yAtemp, tAtemp, xBtemp, yBtemp, tBtemp;
    int xMinA,yMinA,tMinA;
    int xMinB,yMinB,tMinB;
    int sumOcc, occA;
    float tempFloat;
    float beta = 50.0;
    float *ptrA, *ptrB;

    if ( ((imgVolA->patchSizeX) != (imgVolB->patchSizeX)) || ((imgVolA->patchSizeY) != (imgVolB->patchSizeY)) ||
        ((imgVolA->patchSizeT) != (imgVolB->patchSizeT)) )
    {
        MY_PRINTF("Error in ssd_minimum_value, the patch sizes are not equal.\n");
        return -1;
    }

    if (params->patchIndexing == 0)
    {
        xMinA = xA-imgVolA->hPatchSizeX;
        yMinA = yA-imgVolA->hPatchSizeY;
        tMinA = tA-imgVolA->hPatchSizeT;
        xMinB = xB-imgVolB->hPatchSizeX;
        yMinB = yB-imgVolB->hPatchSizeY;
        tMinB = tB-imgVolB->hPatchSizeT;
    }
    else if(params->patchIndexing == 1)
    {
        xMinA = xA;
        yMinA = yA;
        tMinA = tA;
        xMinB = xB;
        yMinB = yB;
        tMinB = tB;
    }

    sumOcc = 0;

    if (params->partialComparison)
    {
        for (k=0; k<imgVolA->patchSizeT; k++)
            for (j=0; j<imgVolA->patchSizeY; j++)
                for (i=0; i<imgVolA->patchSizeX; i++)
                {
                    xAtemp = xMinA + i;
                    yAtemp = yMinA + j;
                    tAtemp = tMinA + k;
                    xBtemp = xMinB + i;
                    yBtemp = yMinB + j;
                    tBtemp = tMinB + k;
                    //do not compare the edges in any case
                    if ((!check_in_boundaries( imgVolA, xAtemp, yAtemp, tAtemp,params)))
                        continue;   //if we are not in the inner boundaries, do not compare

                    sumOcc = sumOcc + (int)(!(*(occVol->get_value_ptr(xAtemp,yAtemp,tAtemp,0))) == 1);

                }
    }
    else    //calculate patch size
    {
        int patchSizeXtemp = min_int(xA + imgVolA->hPatchSizeX,imgVolA->xSize-1) - max_int(xA - imgVolA->hPatchSizeX,0) + 1;
        int patchSizeYtemp = min_int(yA + imgVolA->hPatchSizeY,imgVolA->ySize-1) - max_int(yA - imgVolA->hPatchSizeY,0) + 1;
        int patchSizeTtemp = min_int(tA + imgVolA->hPatchSizeT,imgVolA->tSize-1) - max_int(tA - imgVolA->hPatchSizeT,0) + 1;

        sumOcc = patchSizeXtemp * patchSizeYtemp * patchSizeTtemp;
    }
    sumOcc = max_int(sumOcc,1);

    for (k=tMin; k<tMax; k++)
        for (j=yMin; j<yMax; j++)
            for (i=xMin; i<xMax; i++)
            {
                xAtemp = xMinA + i;
                yAtemp = yMinA + j;
                tAtemp = tMinA + k;

                xBtemp = xMinB + i;
                yBtemp = yMinB + j;
                tBtemp = tMinB + k;


                //do not compare if we are not in the boundaries
                if ((!check_in_boundaries( imgVolA, xAtemp, yAtemp, tAtemp,params)))
                    occA = 1;
                else
                    occA = 0;
                if (occA == 1)
                    continue;   //we do not wish to compare this pixel
                //if we want partial patch comparison
                if (params->partialComparison && occVol->xSize >0)
                    occA = (int)(*(occVol->get_value_ptr(xAtemp, yAtemp, tAtemp,0)) == 1);
                if (occA == 1)
                    continue;   //we do not wish to compare this pixel

                ptrA = imgVolA->get_begin_ptr(xAtemp, yAtemp, tAtemp);
                ptrB = imgVolB->get_begin_ptr(xBtemp, yBtemp, tBtemp);

                //similarity //new
              for (p=0; p<imgVolA->nTupleSize; p++)
           //     for (p=0; p<1; p++)
                {
                    tempFloat = (*(ptrA+p)) - (*(ptrB+p));
                    ssd = ssd - (((tempFloat)*(tempFloat))/sumOcc);
                }
                if( params->gradX != NULL)
                {
                    float normGradXtemp = *((params->normGradX) + (tAtemp)*(imgVolA->xSize)*(imgVolA->ySize)
                             + (xAtemp)*(imgVolA->ySize) + yAtemp) -
                            *((params->normGradX) + (tBtemp)*(imgVolB->xSize)*(imgVolB->ySize)
                             + (xBtemp)*(imgVolB->ySize) + yBtemp);

                    float normGradYtemp = *((params->normGradY) + (tAtemp)*(imgVolA->xSize)*(imgVolA->ySize)
                             + (xAtemp)*(imgVolA->ySize) + yAtemp) -
                            *((params->normGradY) + (tBtemp)*(imgVolB->xSize)*(imgVolB->ySize)
                             + (xBtemp)*(imgVolB->ySize) + yBtemp);

                    ssd = ssd - beta*normGradXtemp*normGradXtemp/sumOcc;
                    ssd = ssd - beta*normGradYtemp*normGradYtemp/sumOcc;
                }
            }

    switch(dirPropagate){
    case -1: //Left
            xMax=imgVolA->patchSizeX;
            xMin=xMax-1;
        break;
    case -2:
        yMax=imgVolA->patchSizeY;
        yMin=yMax-1;
        break;
    case -3:
        tMax=imgVolA->patchSizeT;
        tMin=tMax-1;
        break;
    case 1:
        xMax=1;
        xMin=0;
        break;
    case 2:
        yMax=1;
        yMin=0;
        break;
    case 3:
        tMax=1;
        tMin=0;
        break;
    }
    for (k=tMin; k<tMax; k++)
        for (j=yMin; j<yMax; j++)
            for (i=xMin; i<xMax; i++)
            {
                xAtemp = xMinA + i;
                yAtemp = yMinA + j;
                tAtemp = tMinA + k;

                xBtemp = xMinB + i;
                yBtemp = yMinB + j;
                tBtemp = tMinB + k;


                //do not compare if we are not in the boundaries
                if ((!check_in_boundaries( imgVolA, xAtemp, yAtemp, tAtemp,params)))
                    occA = 1;
                else
                    occA = 0;
                if (occA == 1)
                    continue;   //we do not wish to compare this pixel
                //if we want partial patch comparison
                if (params->partialComparison && occVol->xSize >0)
                    occA = (int)(*(occVol->get_value_ptr(xAtemp, yAtemp, tAtemp,0)) == 1);
                if (occA == 1)
                    continue;   //we do not wish to compare this pixel

                ptrA = imgVolA->get_begin_ptr(xAtemp, yAtemp, tAtemp);
                ptrB = imgVolB->get_begin_ptr(xBtemp, yBtemp, tBtemp);

                // similarity //new
              for (p=0; p<imgVolA->nTupleSize; p++)
             //   for (p=0; p<1; p++)
                {
                    tempFloat = (*(ptrA+p)) - (*(ptrB+p));
                    ssd = ssd + (((tempFloat)*(tempFloat))/sumOcc);
                }
                if( params->gradX != NULL)
                {
                    float normGradXtemp = *((params->normGradX) + (tAtemp)*(imgVolA->xSize)*(imgVolA->ySize)
                             + (xAtemp)*(imgVolA->ySize) + yAtemp) -
                            *((params->normGradX) + (tBtemp)*(imgVolB->xSize)*(imgVolB->ySize)
                             + (xBtemp)*(imgVolB->ySize) + yBtemp);

                    float normGradYtemp = *((params->normGradY) + (tAtemp)*(imgVolA->xSize)*(imgVolA->ySize)
                             + (xAtemp)*(imgVolA->ySize) + yAtemp) -
                            *((params->normGradY) + (tBtemp)*(imgVolB->xSize)*(imgVolB->ySize)
                             + (xBtemp)*(imgVolB->ySize) + yBtemp);

                    ssd = ssd + beta*normGradXtemp*normGradXtemp/sumOcc;
                    ssd = ssd + beta*normGradYtemp*normGradYtemp/sumOcc;
                }

        if ((minVal != -1) && (ssd > minVal))
                {
            return(-1);
                }
            }
    return(ssd);
}
