#include "R.h"
#ifdef ES_DEV_ENV
	#include "..\..\..\IMat\RFunc.h"	 
#else
	#include "RFunc.h"	 
#endif


//#define ZERO_TOL 1e-16

	BOOL RestrictEigenValues (const IMatD &mEV, const IVecD & vClustSize, double dFact, IVecD &vTempNPp2, double dZeroTol) ;
	BOOL RestrictEigenValues_deter (const IMatD &mEV, const IVecD & vClustSize, double dFact, IVecD &vTempNPp2, double dZeroTol) ;
