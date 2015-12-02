//this function calculates the approximate nearest neighbour field for patch match
//in a volume, with multivalued pixels

#include "patch_match.h"
// T often a pointer of float data contain data of class nTupleVolume
template <class T>
// Wrap function for using PatchMatch in Matlab:
nTupleVolume<T> * wrapper_patch_Match_ANN(T *videoA, T *videoB, T *videoFirstGuess, T *videoOcc, T *videoMod,
                                    int xSizeA, int ySizeA, int tSizeA, 
									int xSizeB, int ySizeB, int tSizeB, int nTupleSize,
                                    int patchSizeX, int patchSizeY, int patchSizeT,
                                    const parameterStruct *params)
{
	//declarations
	nTupleVolume<T> *imgVolA, *imgVolB, *dispField, *firstGuessVol,*imgVolOcc, *imgVolMod;
    // reconstruct object of class nTupleVolume from input data
	imgVolA = new nTupleVolume<T>(nTupleSize, xSizeA, ySizeA, tSizeA, patchSizeX, patchSizeY, patchSizeT,1,videoA);
	imgVolB = new nTupleVolume<T>(nTupleSize, xSizeB, ySizeB, tSizeB, patchSizeX, patchSizeY, patchSizeT,1,videoB);

    if (videoFirstGuess == NULL)    //if no guess (first frame)
	{
		firstGuessVol = new nTupleVolume<T>();  //creat a guessVolume
	}
    else                            //use old guess data to creat guessVolume
	{
		firstGuessVol = new nTupleVolume<T>(NDIMS+1, xSizeA, ySizeA, tSizeA, patchSizeX, patchSizeY, patchSizeT,1,videoFirstGuess);
	}
    //get Data from firstGuessVol
	for (int k=0; k<tSizeA; k++)        //temp dimension
		for (int j=0; j<ySizeA; j++)    // y dimension
			for (int i=0; i<xSizeA; i++)// x dimension
				for (int p=0; p<4; p++) // n tuple (4 because shiftGuess have x,y,z vector and ssd error (use for shiftVol only)
				{
					float floatTemp = firstGuessVol->get_value(i,j,k,p);
				}

	//put value in imgVolOcc
	if (videoOcc == NULL)
	{
		imgVolOcc = new nTupleVolume<T>();
	}
	else
	{
		imgVolOcc = new nTupleVolume<T>(1, xSizeB, ySizeB, tSizeB, patchSizeX, patchSizeY, patchSizeT, 1, videoOcc);
	}
    
    	//put value in imgVolMod
	if (videoMod == NULL)
	{
		imgVolMod = new nTupleVolume<T>();
	}
	else
	{
		imgVolMod = new nTupleVolume<T>(1, xSizeA, ySizeA, tSizeA, patchSizeX, patchSizeY, patchSizeT,1 ,videoMod);
	}

    //start calculate dispField using patch_match_ANN function
   // mexPrintf("Before patch_match_ANN\n");
	dispField = patch_match_ANN(imgVolA, imgVolB, firstGuessVol, imgVolOcc,imgVolMod, params);
  //  mexPrintf("After patch_match_ANN\n");
    
    //if dissField ==NULL return NULL else return disField and free the memory
    if (dispField == NULL)
    {
        delete imgVolA;
        delete imgVolB;
        delete firstGuessVol;
        delete imgVolOcc;
        delete imgVolMod;

        return NULL;
    }
    //mexPrintf("After patch_match_ANN\n");
    //show_nTuple_volume(dispField);
	//put the displacement values in the output array

	//free memory
    //mexPrintf("Before destruction imgVolA\n");
	delete imgVolA;
    //mexPrintf("Before destruction imgVolB\n");
	delete imgVolB;
    //mexPrintf("Before destruction firstGuessVol\n");
    delete firstGuessVol;
    //mexPrintf("Before destruction imgVolOcc\n");
	delete imgVolOcc;
    //mexPrintf("Before destruction imgVolMod\n");
    delete imgVolMod;
    //mexPrintf("Before destruction dispField\n");
	return(dispField);
}

//this function calculates a nearest neighbour field, from imgVolA to imgVolB
template <class T>
nTupleVolume<T> * patch_match_ANN(nTupleVolume<T> *imgVolA, nTupleVolume<T> *imgVolB, 
        nTupleVolume<T> *firstGuessVol, nTupleVolume<T> *imgVolOcc,nTupleVolume<T> *imgVolMod,
        const parameterStruct *params)
{
	//decalarations
	int xSizeA, ySizeA, tSizeA, xSizeB, ySizeB, tSizeB, nTupleSize;
	int i, nbModified;
    clock_t startTime;
    double propagationTime,randomSearchTime;
	nTupleVolume<T> *dispField;

	//get image volume sizes
	xSizeA = imgVolA->xSize;
	ySizeA = imgVolA->ySize;
	tSizeA = imgVolA->tSize;
    nTupleSize = imgVolA->nTupleSize;       //nTuple size is 3 for color Image and 4 in shiftVol (the 4th contain SSD error)
	xSizeB = imgVolB->xSize;
	ySizeB = imgVolB->ySize;
	tSizeB = imgVolB->tSize;

	//check certain parameters
	if(nTupleSize != (imgVolB->nTupleSize) )
	{
		MY_PRINTF("Error in patch_match_ANN, the size of the vector associated to each pixel is different for the two image volumes.");
		return NULL;
	}
	if( (imgVolA->patchSizeX != (imgVolB->patchSizeX)) || (imgVolA->patchSizeY != (imgVolB->patchSizeY)) ||
		(imgVolA->patchSizeT != (imgVolB->patchSizeT))  )	//check that the patch sizes are equal
	{
		MY_PRINTF("Error in patch_match_ANN, the size of the patches are not equal in the two image volumes.");
		return NULL;
	}
	if ( ( imgVolA->patchSizeX > imgVolA->xSize) || ( imgVolA->patchSizeY > imgVolA->ySize) || ( imgVolA->patchSizeT > imgVolA->tSize) ||
		( imgVolA->patchSizeX > imgVolB->xSize) || ( imgVolA->patchSizeY > imgVolB->ySize) || ( imgVolA->patchSizeT > imgVolB->tSize)
		)	//check that the patch size is less or equal to each dimension in the images
	{
		MY_PRINTF("Error in patch_match_ANN, the patch size is to large for one or more of the dimensions of the image volumes.");
		return NULL;
	}

	//create the displacement field
    //mexPrintf("Before dispField creation\n");
    dispField = new nTupleVolume<T>(NDIMS+1, xSizeA, ySizeA, tSizeA, imgVolA->patchSizeX, imgVolA->patchSizeY, imgVolA->patchSizeT,1); //empty data dispField
    //mexPrintf("After dispField creation\n");
    //show_nTuple_volume(dispField);
	//randomly initialise the displacement field
    
    //#pragma omp parallel
    //mexPrintf("Hello\n");
    
    //cuda_full_search(dispField,imgVolA,imgVolB,imgVolOcc,imgVolMod);
    //mexPrintf("dispField[1] : %f\n",dispField->values[1]);
    //return(dispField);

    //Calculate time for PatchMatch
    propagationTime = 0.0;          //propagationTime
    randomSearchTime = 0.0;         //Random Search Time
    ssdTime = 0.0;                  // Calculate ssd
    randomGenerationTime = 0.0;
    if (params->fullSearch == 1)
    {
        //memset(dispField->values,0,(size_t)(dispField->xSize)*(dispField->ySize)*(dispField->tSize)*(dispField->nTupleSize)*sizeof(float));
        //Do initialise displace_ment_fiend, if have old dispField than use it, otherwise use random
        initialise_displacement_field(dispField, imgVolA,imgVolB, firstGuessVol, imgVolOcc,params);

        //memset(dispField->values,0,(size_t)(dispField->xSize)*(dispField->ySize)*(dispField->tSize)*(dispField->nTupleSize)*sizeof(float));
        startTime = clock();

        /*Start search (brut-force) search*/
        patch_match_full_search(dispField, imgVolA,imgVolB, imgVolOcc,imgVolMod,params);
        //-----------------//
        propagationTime = propagationTime + double(clock() - startTime);
    }
    else    //normal patchMatch
    {
        mexPrintf("Initialisation PatchMatch\n");
        initialise_displacement_field(dispField, imgVolA, imgVolB, firstGuessVol, imgVolOcc,params);
        //show_nTuple_volume(dispField);

       //check if the disField is valid dispField (occlusion, lie in inner boundary..)
        if (check_disp_field(dispField, imgVolA, imgVolB,imgVolOcc,params) == -1)
            return(dispField);
        //Start propagation:
        for (i=0; i<(params->nIters); i++)
       {
            startTime = clock();
            //Start propagation:
            nbModified = patch_match_propagation(dispField, imgVolA, imgVolB,imgVolOcc,imgVolMod,params,i);
            propagationTime = propagationTime + double(clock() - startTime);
            
            startTime = clock();
            //Start random search (with weight)
            nbModified = patch_match_random_search(dispField, imgVolA,imgVolB, imgVolOcc,imgVolMod,params);
            randomSearchTime = randomSearchTime + double(clock() - startTime);
        }
    }
    MY_PRINTF("Propagation time : %f s\n",propagationTime/CLOCKS_PER_SEC);
    MY_PRINTF("Random search time : %f s\n",randomSearchTime/CLOCKS_PER_SEC);

	return(dispField);
}
