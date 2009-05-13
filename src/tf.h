#ifdef _MSC_VER
	#include "..\..\IMat\RFunc.h"	 
	#include "..\..\IMat\IMat.h"
	#include "..\..\IMat\ITens.h"
#else
	#include "RFunc.h"	 
	#include "IMat.h"
	#include "ITens.h"
#endif

	class TClust
	{
	public:

		TClust (DWORD *pdwParIn, DWORD *pdwParOut, double *pdParIn, double *pdParOut, double *pdX, double *pdM, double *pdS, DWORD *pdwAssign, DWORD *pdwClustSize) ;
		BOOL tclust () ;
	protected:

		BOOL trimEval () ;

		void SetAllCovmatsIdent () ;

		void FindInitClustAssignment () ;
		void FindInitClustSize () ;
		void FindInitClustSize_R () ;

		void CalcDensity (const IMatD &mX, const IVecD &vDens, const IVecD &vCurM, const IVecD &vEVal, const IMatD &mEVec, const double dFact = 1) ;
		void FindNearestClust (const IVecD &vM, const IVecDW &vW) ;
		BOOL FindClustAssignment (BOOL bLastRun = FALSE) ;
		void CalcClustParams () ;
		DWORD ExistFlatClust () ;
		double CalcObjFunc () ;

		BOOL SumUpSigma () ;

		void CalcClusterSize () ;
		void FindOutliers () ;

		BOOL SaveCurResult (double dCurObj, DWORD dwCode = 0) ;

		DWORD	&m_n, &m_p, &m_K, &m_dwIter, &m_dwKSteps, &m_dwEqualWeights, &m_dwRestr, &m_dwTrace, 
				&m_dwConvCount, &m_dwIterSuccess, &m_dwCode ;//, &m_dwMinMode;

		double	&m_dAlpha, &m_dFactor, &m_dZeroTol, 
				&m_dBestObj ;

		const double  m_dDensFact ;
		const DWORD	m_dwNoTrim, m_dwTrim ;

		IVecDW	m_vInd, m_vIndBest, m_vIndOld, m_vCurInd, m_vClustSizeBest, m_vRank ;

		ITens (double, 3) m_tEVec, m_tCurS, m_tBestS, m_tBestEVec ;

		IMatD	m_mCurM, m_mBestM, m_mX, m_mLL, m_mEVal, m_mTempNP1, m_mTempNP2/*, m_mXTrim*/, m_mBestEVal ;
		IVecD	m_vWeights, m_vClustSize, m_vDisc, m_vDiscSorted, m_vTempN1, m_vTempN2, m_vTempNPp2 ;

		IMatArrayD m_amEVec, m_amCurS ;
		IVecArrayD m_avEVal, m_avCurM, m_avLL ;
	} ;
