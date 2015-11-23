/*this function estimates the colours of pixels, given a displacement field
showing th*/

#include "reconstruct_video.h"

/*this fuction do a reconstruction of video after calculating ANN field.*/
template <class T>
void reconstruct_video(nTupleVolume<T>* imgVol, nTupleVolume<T>* occVol,
        nTupleVolume<T>* dispField, float sigmaColour, int useAllPatches, int reconstructionType)
{
	/*decalarations*/
    int xSize, ySize, tSize, nTupleSize; //img Size
    int i,j,k;                           //3-D index
    int iMin,iMax,jMin,jMax,kMin,kMax;   //Boundary of index Volume
    int ii,jj,kk, weightInd;
    int xDisp, yDisp, tDisp,xDispShift,yDispShift,tDispShift;
    int shiftI,shiftJ,shiftK;
    int hPatchSizeX,hPatchSizeY,hPatchSizeT; //half patchsize
    int hI,hJ,hK,doubleII,doubleJJ,doubleKK;
    int nbNeighbours;
    int correctInfo;
    float alpha, adaptiveSigma;
    float *weights,sumWeights, avgColourR, avgColourG, avgColourB, *colours;

	/*get image volume sizes*/
	xSize = imgVol->xSize;
	ySize = imgVol->ySize;
	tSize = imgVol->tSize;
	nTupleSize = imgVol->nTupleSize;

    hPatchSizeX = imgVol->hPatchSizeX;          //half patchsize
    hPatchSizeY = imgVol->hPatchSizeY;
    hPatchSizeT = imgVol->hPatchSizeT;
    
    /*allocate the (maximum) memory for the weights*/
    nbNeighbours = (imgVol->patchSizeX)*(imgVol->patchSizeY)*(imgVol->patchSizeT); //number of neighbor
    weights = (float*)malloc((size_t)(nbNeighbours*sizeof(float)));
    colours = (float*)malloc((size_t)(NCHANNELS*nbNeighbours*sizeof(float)));
    
	/*check certain parameters*/
	if( (imgVol->patchSizeX != (imgVol->patchSizeX)) || (imgVol->patchSizeY != (imgVol->patchSizeY)) ||
        (imgVol->patchSizeT != (imgVol->patchSizeT))  )	/*check that the patch sizes are equal*/ //what is this for?
	{
		MY_PRINTF("Error in estimate_colour, the size of the patches are not equal in the two image volumes.");
		return;
	}
	if ( ( imgVol->patchSizeX > imgVol->xSize) || ( imgVol->patchSizeY > imgVol->ySize) || ( imgVol->patchSizeT > imgVol->tSize) )	/*check that the patch size is less or equal to each dimension in the images*/
	{
		MY_PRINTF("Error in estimate_colour, the patch size is to large for one or more of the dimensions of the image volume.");
		return;
	}

    for (k=0; k<(occVol->tSize); k++)       //run all over occVolume, t dimension
        for (j=0; j<(occVol->ySize); j++)   // y dimension
            for (i=0; i<(occVol->xSize); i++)//x dimension
            {    
                if ( ((occVol->get_value(i,j,k,0)) == 0) || ((occVol->get_value(i,j,k,0) == 2) )  )     //do not inpaint
                    continue;
                else    /*an occluded pixel (therefore to be modified)*/
                {
                    if (reconstructionType == NEAREST_NEIGHBOUR_RECONSTRUCTION )    //nearest neighbor reconstruction
                    {
                        xDisp = i + (int)dispField->get_value(i,j,k,0);             // x Displacement
                        yDisp = j + (int)dispField->get_value(i,j,k,1);             // y Displacement
                        tDisp = k + (int)dispField->get_value(i,j,k,2);             // t displacement

                        ////if pure replacing of pixels
                        copy_pixel_values_nTuple_volume(imgVol, imgVol,xDisp, yDisp, tDisp, i, j, k); //copy pixel value of ANN to occluded pixel.
                        ///set_value_nTuple_volume(occVol,i,j,k,2,0);
                        ///set_value_nTuple_volume(imgVol,i,j,k,0,0);
                        continue;
                    }
                     
                    // init array weight and colours
                    for (ii=0;ii<(imgVol->patchSizeX)*(imgVol->patchSizeY)*(imgVol->patchSizeT); ii++) //all the pixel in a patch
					{
                        weights[ii] = (float)-1; //weight for all pixel in patch
                        colours[ii] = (float)-1; //red
                        colours[ii + nbNeighbours] = (float)-1;//green chanel
                        colours[ii + 2*nbNeighbours] = (float)-1; //blue chanel
					}

                    sumWeights = 0.0;
                    alpha = FLT_MAX;
                    correctInfo = 0;
                    avgColourR = 0.0;
                    avgColourG = 0.0;
                    avgColourB = 0.0;
                    
                    //boudary for each patch
                    iMin = max_int(i - hPatchSizeX,0);
                    iMax = min_int(i + hPatchSizeX,(imgVol->xSize)-1 );
                    jMin = max_int(j - hPatchSizeY,0);
                    jMax = min_int(j + hPatchSizeY,(imgVol->ySize)-1 );
                    kMin = max_int(k - hPatchSizeT,0);
                    kMax = min_int(k + hPatchSizeT,(imgVol->tSize)-1 );
                    
                    /*
                    MY_PRINTF("iMin : %d, iMax : %d\n",iMin,iMax);
                    MY_PRINTF("jMin : %d, jMax : %d\n",jMin,jMax);
                    MY_PRINTF("kMin : %d, kMax : %d\n",kMin,kMax);*/
                    /*first calculate the weights*/
                    //run for each patch, calculate the spatial temporal neighbor of the patch
                    for (kk=kMin; kk<=kMax;kk++) //t dimension
                        for (jj=jMin; jj<=jMax;jj++)//y dimension
                            for (ii=iMin; ii<=iMax;ii++) //x dimension
                            {
                                /*get ssd similarity*/
                                xDisp = ii + (int)dispField->get_value(ii,jj,kk,0); //get value dispX for pixel in the patch
                                yDisp = jj + (int)dispField->get_value(ii,jj,kk,1);
                                tDisp = kk + (int)dispField->get_value(ii,jj,kk,2);
                                /*(spatio-temporally) shifted values of the covering patches*/
                                xDispShift = xDisp - (ii-i); //i + dispField (dispField for pixel (ii,jj,kk))
                                yDispShift = yDisp - (jj-j);
                                tDispShift = tDisp - (kk-k);
                        
                                 if (useAllPatches == 1)    // if use all patches
                                 {
                                     
                                    alpha = (float)min_float(dispField->get_value(ii,jj,kk,3),alpha); //get min SSD and assign to alpha
                                    //transform matrice indice to array indice.
                                    weightInd = (int)((kk-kMin)*(imgVol->patchSizeX)*(imgVol->patchSizeY) + (jj-jMin)*(imgVol->patchSizeX) + ii-iMin);
                                    weights[weightInd] = dispField->get_value(ii,jj,kk,3); //weight is an array, SSD error of patch
                                    
                                    //colours is an array
                                    colours[weightInd] = (float)(imgVol->get_value(xDispShift,yDispShift,tDispShift,0)); //red
                                    colours[weightInd + nbNeighbours] = (float)(imgVol->get_value(xDispShift,yDispShift,tDispShift,1));//green
                                    colours[weightInd + 2*nbNeighbours] = (float)(imgVol->get_value(xDispShift,yDispShift,tDispShift,2)); //blue
                                    correctInfo = 1;
                                 }
                                 else   /*only use some of the patches*/
                                 {
                                     if (((occVol->get_value(ii,jj,kk,0)) == 0) || (occVol->get_value(ii,jj,kk,0) ==-1))
                                     {
                                        alpha = (float)min_float(dispField->get_value(ii,jj,kk,3),alpha); 
                                        weightInd = (int)((kk-kMin)*(imgVol->patchSizeX)*(imgVol->patchSizeY) + (jj-jMin)*(imgVol->patchSizeX) + ii-iMin);
                                        weights[weightInd] = dispField->get_value(ii,jj,kk,3);
                                        
                                        colours[weightInd] = (float)(imgVol->get_value(xDispShift,yDispShift,tDispShift,0));
                                        colours[weightInd + nbNeighbours] = (float)(imgVol->get_value(xDispShift,yDispShift,tDispShift,1));
                                        colours[weightInd + 2*nbNeighbours] = (float)(imgVol->get_value(xDispShift,yDispShift,tDispShift,2));
                                        correctInfo = 1;
                                     }
                                     else
                                     {
                                        weightInd = (int)((kk-kMin)*(imgVol->patchSizeX)*(imgVol->patchSizeY) + (jj-jMin)*(imgVol->patchSizeX) + ii-iMin);
                                        weights[weightInd] = -1;
                                         
                                        colours[weightInd] = -1;
                                        colours[weightInd + nbNeighbours] = -1;
                                        colours[weightInd + 2*nbNeighbours] = -1;
                                        continue;
                                     }
                                 }
                            }
                    
                    alpha = max_float(alpha,1);
                    if (correctInfo == 0)
                        continue;
                    
                    if (reconstructionType == BEST_PATCH_RECONSTRUCTION)
                    {
                        estimate_best_colour(imgVol,imgVol, weights, nbNeighbours, colours, sigmaColour, i, j, k);
                        continue;
                    }
                    //get the 75th percentile of the distances for setting the adaptive sigma
                    adaptiveSigma = get_adaptive_sigma(weights,(imgVol->patchSizeX)*(imgVol->patchSizeY)*(imgVol->patchSizeT),sigmaColour);
					adaptiveSigma = max_float(adaptiveSigma,(float)0.1);
                    
                    /* ///MY_PRINTF("alpha : %f\n",alpha);
                    //adjust the weights : note, the indices which are outside the image boundaries
                    //will have no influence on the final weights (they are initialised to 0)  */
                    for (kk=kMin; kk<=kMax;kk++)
                        for (jj=jMin; jj<=jMax;jj++)
                            for (ii=iMin; ii<=iMax;ii++)
                            {
                                if (useAllPatches)
                                {
                                    /*weights = exp( -weights/(2*sigma�*alpha))*/
                                    weightInd = (int)((kk-kMin)*(imgVol->patchSizeX)*(imgVol->patchSizeY) + (jj-jMin)*(imgVol->patchSizeX) + ii-iMin);
                                    weights[weightInd] = (float)(exp( - ((weights[weightInd])/(2*adaptiveSigma*adaptiveSigma)) ));/*exp( - ((weights[ii])/(2*sigmaColour*sigmaColour*alpha)) );*/
                                    //
                                    sumWeights = (float)(sumWeights+weights[weightInd]);
                                }
                                else   /*only use some of the patches*/
                                {
                                     if (((occVol->get_value(ii,jj,kk,0)) == 0) || (occVol->get_value(ii,jj,kk,0) ==-1))
                                     {
                                        /*weights = exp( -weights/(2*sigma�*alpha))*/
                                        weightInd = (int)((kk-kMin)*(imgVol->patchSizeX)*(imgVol->patchSizeY) + (jj-jMin)*(imgVol->patchSizeX) + ii-iMin);
                                        weights[weightInd] = (float)(exp( - ((weights[weightInd])/(2*adaptiveSigma*adaptiveSigma)) ));/*exp( - ((weights[ii])/(2*sigmaColour*sigmaColour*alpha)) );*/
                                        //
                                        sumWeights = (float)(sumWeights+weights[weightInd]);
                                     }
                                     else
                                         continue;
                                }
                            }

                    /*now calculate the pixel value(s)*/
                    for (kk=kMin; kk<=kMax;kk++)
                        for (jj=jMin; jj<=jMax;jj++)
                            for (ii=iMin; ii<=iMax;ii++)
                            {
                                if (useAllPatches)
                                {
                                    weightInd = (int)((kk-kMin)*(imgVol->patchSizeX)*(imgVol->patchSizeY) + (jj-jMin)*(imgVol->patchSizeX) + ii-iMin);
                                    /*get ssd similarity*/
                                    xDisp = ii + (int)dispField->get_value(ii,jj,kk,0);
                                    yDisp = jj + (int)dispField->get_value(ii,jj,kk,1);
                                    tDisp = kk + (int)dispField->get_value(ii,jj,kk,2);
                                    /*(spatio-temporally) shifted values of the covering patches*/
                                    xDispShift = xDisp - (ii-i);
                                    yDispShift = yDisp - (jj-j);
                                    tDispShift = tDisp - (kk-k);
                                    avgColourR = avgColourR + (float)(weights[weightInd])*(imgVol->get_value(xDispShift,yDispShift,tDispShift,0));
                                    avgColourG = avgColourG + (float)(weights[weightInd])*(imgVol->get_value(xDispShift,yDispShift,tDispShift,1));
                                    avgColourB = avgColourB + (float)(weights[weightInd])*(imgVol->get_value(xDispShift,yDispShift,tDispShift,2));
                                }
                                else
                                {
                                     if (((occVol->get_value(ii,jj,kk,0)) == 0) || (occVol->get_value(ii,jj,kk,0) ==-1))
                                     {
                                        weightInd = (int)((kk-kMin)*(imgVol->patchSizeX)*(imgVol->patchSizeY) + (jj-jMin)*(imgVol->patchSizeX) + ii-iMin);
                                        /*get ssd similarity*/
                                        xDisp = ii + (int)dispField->get_value(ii,jj,kk,0);
                                        yDisp = jj + (int)dispField->get_value(ii,jj,kk,1);
                                        tDisp = kk + (int)dispField->get_value(ii,jj,kk,2);
                                        /*(spatio-temporally) shifted values of the covering patches*/
                                        xDispShift = xDisp - (ii-i);
                                        yDispShift = yDisp - (jj-j);
                                        tDispShift = tDisp - (kk-k);
                                        avgColourR = avgColourR + (float)(weights[weightInd])*(imgVol->get_value(xDispShift,yDispShift,tDispShift,0));
                                        avgColourG = avgColourG + (float)(weights[weightInd])*(imgVol->get_value(xDispShift,yDispShift,tDispShift,1));
                                        avgColourB = avgColourB + (float)(weights[weightInd])*(imgVol->get_value(xDispShift,yDispShift,tDispShift,2));
                                     }
                                     else
                                         continue;
                                }
                            }
                         /*MY_PRINTF("SumWeights : %f\n",sumWeights);*/
                    imgVol->set_value(i,j,k,0,(T)(avgColourR/(sumWeights)));
                    imgVol->set_value(i,j,k,1,(T)(avgColourG/(sumWeights)));
                    imgVol->set_value(i,j,k,2,(T)(avgColourB/(sumWeights)));
                    /*set_value_nTuple_volume(occVol,i,j,k,0,0);*/
                }
            }

        free(weights);
        free(colours);
        return;
}
