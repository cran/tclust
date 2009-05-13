#include <math.h>
#include <R.h>

#include "tf.h"

#ifdef _MSC_VER
	#include "..\..\..\IMat\RFunc.h"	 
	#include "..\..\..\IMat\IMat.h"
#else
	#include "RFunc.h"	 
	#include "IMat.h"
#endif

#include <R_ext/Applic.h>
#include <Rmath.h>

#include "restr.h"

	void TClust::FindInitClustSize ()
	{

		if (m_dwEqualWeights)			//	for equal sized clusters
		{
			m_vClustSize.Reset (m_dwNoTrim / (double) m_K) ;
			m_vWeights.Reset (1.0 / (double) m_K) ;
		}
		else	
		{
			IVecD v (m_K) ;
			DWORD k, n ;

			for (k = m_K - 1; k != (DWORD) -1; k--)
				v(k) = runif (0, 1) ;		//	get array of 0-1 unif distr. numbers
			sort (v) ;						//	sort
			for (k = m_K - 2; k != (DWORD) -1; k--)	//	calc reverse cumsum
				v (k) += v(k + 1) ;

											//	sorting is good when m_dwNoTrim is huge.
											//	when checking which cluster a guess belongs to, the big clusters are checked earlier

			v /= v(0) ;						//	norm vector v (sum = 1 - i.e. first element is 1, as this these are cumulated sums)

			m_vClustSize.Reset (0) ;
			for (n = m_dwNoTrim; n > 0; n--)			//	
			{
				double dCur = runif (0, 1) ;	//	get rand. cluster assignment
														//		1- runif for R compatibility (gives similar results - still differences - why?)
				for (k = m_K - 1; k != (DWORD) -1; k--)	//	check in which group it can be found
					if (dCur <= v(k))
					{
						m_vClustSize (k) += 1 ;
						break ;
					}
			}

			m_vWeights = m_vClustSize / m_dwNoTrim ;
		}
	}

	void TClust::FindInitClustSize_R ()
	{	//	as FindInitClustSize, but returns the same results as the current R - code 2 be deleted!!

		if (m_dwEqualWeights)			//	for equal sized clusters
		{
			m_vClustSize.Reset (m_dwNoTrim / (double) m_K) ;
			m_vWeights.Reset (1.0 / (double) m_K) ;
		}
		else	
		{
			IVecD v (m_K) ;

			DWORD k, n ;
			
			//for (k = 0; k < m_K; k++)
			for (k = m_K - 1; k != (DWORD) -1; k--)
				v(k) = runif (0, 1) ;		//	calc cumulated sum of cluster propability (proportional cluster size)

			IVecDW Vorder (m_K) ;
			v.order (Vorder) ;
			sort (v) ;

			for (k = m_K - 2; k != (DWORD) -1; k--)
				v (k) += v(k + 1) ;
			v /= v(0) ;
//RtprintVec (v);										//	this is as in R (reverse order)
			m_vClustSize.Reset (0) ;
			for (n = m_dwNoTrim; n > 0; n--)			//	
			{
				double dCur = runif (0, 1) ;	//	get rand. cluster assignment
														//		1- runif for R compatibility (gives similar results - still differences - why?)
				for (k = m_K - 1; k != (DWORD) -1; k--)	//	check in which group it can be found
					if (dCur <= v(k))
					{
						m_vClustSize (k) += 1 ;
						break ;
					}
			}

			

			IVecD vtemp (m_K) ;
			vtemp << m_vClustSize ;

			for (k = m_K - 1; k != (DWORD) -1; k--)
				m_vClustSize (m_K - 1 - Vorder (k)) = vtemp (k) ;

			m_vWeights = m_vClustSize / m_dwNoTrim ;
if (int(m_dwTrace) >= 10)
{
	RtprintVec (m_vClustSize, "ClustSize:\t");
	RtprintVec (m_vWeights, "ClustWeights:\t");
}
		}
	}

	void TClust::FindInitClustAssignment ()
	{		//	finds initial cluster assignment
		DWORD k ;

		IVecD vP (m_K) ;

		for (k = 0; k < m_K; k++)
		{	//	for all clusters
			//	finds p+1 observations for forming the initial cluster
			SampleNoReplace (m_p + 1, m_mX.nrow (), m_vCurInd) ;

			ISub s (m_vCurInd) ;
			IMatD curX (m_mX (s/*ISub (m_vCurInd)*/, ISub ()))  ;

			//	calcs mean & cov struct. of selected observations.
			colMeans (curX, m_avCurM[k]) ;
			cov (curX, m_amCurS[k]) ;
		}

		m_tCurS *= (m_p) / (m_p  + 1.0) ;
	}

	BOOL TClust::SumUpSigma ()
	{
		DWORD i ;

		m_amCurS [0] *= m_vWeights (0) ;
		for (i = m_K - 1; i; i--)
			m_amCurS [0] += (m_amCurS [i] *= m_vWeights (i)) ;

		eigen_sqr (m_amCurS[0], m_avEVal [0], m_amEVec [0], FALSE) ;

		m_avEVal [0].Limit_L (0) ;

		for (i = m_K - 1; i; i--)
		{
			m_amCurS [i] << m_amCurS [0] ;
			m_avEVal [i] << m_avEVal [0] ;
			m_amEVec [i] << m_amEVec [0] ;
		}

		return m_avEVal[0].max () > 0 ;
	}

	void TClust::CalcClustParams ()
	{		//	calculates the mean and the cov matrix of the clusters currently defined by m_vIndTrim
		DWORD k ;

//		for (k = m_K - 1; k != (DWORD) -1; k--)
		for (k = 0; k < m_K; k++)
		{
			IMatD mCurX (m_mX (ISub (m_vInd == k), ISub ())) ;

			if (!m_dwEqualWeights)
				m_vWeights (k) = m_vClustSize (k) / (double) m_dwNoTrim ;

			if (m_vClustSize (k) > 0)
				colMeans (mCurX, m_avCurM [k]) ;

			if (m_vClustSize (k) > 1)
			{
				cov (mCurX, m_amCurS[k]) ;
				m_amCurS[k] *= ((m_vClustSize (k) - 1) / m_vClustSize (k)) ;

			}
			else
				m_amCurS[k].Reset (0) ;
		}
	}

	void TClust::CalcDensity (const IMatD &mX, const IVecD &vDens, const IVecD &vCurM, const IVecD &vEVal, const IMatD &mEVec, const double dFact)
	{
		m_mTempNP1.Reshape (mX.nrow (), mX.ncol ()) ;
		FC_ElOp<FC::FC_minus, DWORD>::OpMV_row (mX, vCurM, m_mTempNP1) ;	// TempNP1 = X centered

		m_mTempNP2.Reshape (mX.nrow (), mX.ncol ()) ;
		matmultmat (m_mTempNP1, mEVec, m_mTempNP2) ;						//	X scaled (Z)

		m_vTempN1.Reshape (m_p) ;

		FC_ElOp<FC::FC_pow, DWORD>::OpVE (vEVal, -0.5, m_vTempN1) ;			//	Gamma ^-0.5

		FC_ElOpAs<FC::FC_multiply>::OpMV_row (m_mTempNP2, m_vTempN1) ;		//	Z %*% Gamm ^-0.5

		FC_ElOpAs<FC::FC_sqr>::OpM (m_mTempNP2) ;							//	sqr (Z %*% Gamm ^-0.5)

		rowSums (m_mTempNP2, vDens) ;										//	 == mahalanobis

		FC_ElOpAs<FC::FC_multiply>::OpVE (vDens, -0.5) ;					//	mahalanobis * 0.5

		FC_ElOpAs<FC::FC_exp>::OpV (vDens) ;								//	exp (maha / 2)

		double dDet = prod (m_vTempN1) ;									//	||Sigma||
		vDens *= dDet * dFact * m_dDensFact ;								//
	}

	void TClust::FindNearestClust (const IVecD &vM, const IVecDW &vW)
	{
		DWORD r, c ;

		vW.Reset (0) ;

		for (r = m_mLL.nrow () - 1; r != (DWORD) -1; r--)
		{
			BOOL bFirst = TRUE ;
			for (c = m_mLL.ncol () - 1; c != (DWORD) -1; c--)
			{
				if (bFirst)
				{
					vM(r) = m_mLL (r, c) ;
					vW(r) = c ;
					bFirst = FALSE ;
				}
				else if (setmax (m_mLL (r, c), vM(r)))
					vW(r) = c ;
			}
		}
	}

	void TClust::FindOutliers ()
	{
		m_vDisc.order (m_vRank) ;	//	calculates the rank of m_vDisc and stores it in m_vRank
									//	setting the TrimIdx - array
		DWORD i ;

		for (i = m_dwTrim - 1; i!= (DWORD) -1; i--)
			m_vInd (m_vRank(i)) = -1 ;
	}

	void TClust::CalcClusterSize ()
	{
		m_vClustSize.Reset (0) ;

		DWORD k ;
		for (k = m_n - 1; k != (DWORD) -1; k--)
			if (m_vInd (k) != (DWORD) -1)
				m_vClustSize (m_vInd (k)) ++ ;

		m_vWeights << m_vClustSize ;
		m_vWeights /= m_dwNoTrim ;
	}

	BOOL TClust::FindClustAssignment (BOOL bLastRun) 
	{		//	calculates the vector m_vIndTrim by 
		DWORD k ;

		for (k = 0; k < m_K; k++)
			CalcDensity (m_mX, m_avLL[k], m_avCurM [k], m_avEVal[k], m_amEVec [k], m_vWeights (k)) ;

		FindNearestClust (m_vDisc, m_vInd) ;

		FindOutliers () ;

		if (m_vInd.Equal (m_vIndOld) && !bLastRun)	//	check the change of group assignment AFTER the calculation of outliers. 
			return FALSE ;							//	otherwise the alg will stop after the 1st iteration when k == 1, 
		m_vIndOld << m_vInd ;						//	since there wouldn't even be the change of the change of the cluster assignment

		if (!m_dwEqualWeights)
			CalcClusterSize () ;

		return TRUE ;
	}

	double TClust::CalcObjFunc ()
	{
		double dObj = 0 ;
		DWORD k ;

		for (k = 0; k < m_K; k++)
		{
			IMatD mCurX (m_mX (ISub (m_vInd == k), ISub ())) ;

			DWORD dwCurClustSize = mCurX.nrow () ;
			m_vClustSize (k) = dwCurClustSize ;

			if (!dwCurClustSize)
				continue ;

			m_vTempN2.Reshape (dwCurClustSize) ;
			CalcDensity (mCurX, m_vTempN2, m_avCurM [k], m_avEVal[k], m_amEVec [k]) ;

			FC_ElOpAs<FC::FC_log>::OpV (m_vTempN2) ;

			dObj += sum (m_vTempN2) ;

			if (!m_dwEqualWeights)	//	this can be avoided when equal-sized clusters are expected (only adding a fixed value to each object. function)
				dObj += dwCurClustSize * log (dwCurClustSize / (double) m_dwNoTrim) ;
		}

		return dObj ;
	}

	BOOL TClust::trimEval ()	//	should be called trimCov
	{
		if (m_dwRestr == 2)
			return SumUpSigma () ;

		DWORD k ;
		for (k = 0; k < m_K; k++)
			eigen_sqr (m_amCurS[k], m_avEVal[k], m_amEVec[k], FALSE) ;

		m_mEVal.Limit_L (0) ;

		if (m_dwRestr == 1 && m_p > 1)
		{
			if (!RestrictEigenValues_deter (m_mEVal, m_vClustSize, m_dFactor, m_vTempNPp2, m_dZeroTol))
				return FALSE ;
		}	
		else if (!RestrictEigenValues (m_mEVal, m_vClustSize, m_dFactor, m_vTempNPp2, m_dZeroTol))		//	trims the eigenvalues of the used clusters (given by m_vUsedK)
			return FALSE ;

			//	recalculates (recomposes) all covariance matrices

		for (k = 0; k < m_K; k++)
		{
			m_mTempNP1.Reshape (m_p, m_p) ;
			FC_ElOp<FC::FC_multiply, double>::OpMV_row (m_amEVec[k], m_avEVal[k], m_mTempNP1) ;	//	m_avEVal now is the squareroot of actually trimmed the EV

			matmultmat (m_mTempNP1, t (m_amEVec[k]), m_amCurS [k]) ;
		}
		return TRUE ;
	}

	BOOL TClust::SaveCurResult (double dCurObj, DWORD dwCode)
	{
		m_dBestObj = dCurObj ;
		m_mBestM << m_mCurM ;
		m_tBestS << m_tCurS ;
		FC_As::OpTT (m_vClustSizeBest, m_vClustSize) ;
		m_vIndBest << m_vInd ;
		m_dwCode = dwCode ;
		return TRUE ; 
	}

	void TClust::SetAllCovmatsIdent ()
	{
		//	setting all covariance matrices to ident matrix and correspondingly changes the eigenvalues & vectors

		DWORD k ;
		for (k = 0; k < m_K; k++)
		{
			m_amCurS[k].setdiag () ;
			m_amEVec[k].setdiag2 () ;
		}
		m_mEVal.Reset (1) ;
	}

	BOOL TClust::tclust ()
	{	 
		DWORD &i = m_dwIterSuccess = 0, j ;

		m_vIndOld.Reset ((DWORD) -1) ;

		double dLastObj = 0 ;

		for (i = 0; i < m_dwIter; i++)
		{
			FindInitClustAssignment () ;


			FindInitClustSize_R () ;

			for (j = 0; 1; j++)
			{

				if (!trimEval ())
					if (i)
						return !SaveCurResult (CalcObjFunc (), 2) ;
					else
						SetAllCovmatsIdent () ;
				if (!FindClustAssignment () ||		//	returns false, if the cluster assignment has not changed -> break 
					j == m_dwKSteps)				//	max number of iterations reached -> break
					break ;

				if (int(m_dwTrace) >= 2)
				{
					double dCurObj = CalcObjFunc () ;

					if (j && dLastObj > dCurObj)
						Rprintf ("Objective function dropped from %.10f to %.10f in (%d/%d)\r\n", dLastObj, dCurObj, i, j) ;
					else
						Rprintf ("Objective function %.10f in (%d/%d)\r\n", dCurObj, i, j) ;

					dLastObj = dCurObj ;
				}

				CalcClustParams () ;
			}

			if (j < m_dwKSteps)						//	did this iteration converge?
				m_dwConvCount ++ ;

			double dCurObj = CalcObjFunc () ;

			if (!i || dCurObj > m_dBestObj)			//	store if the current solution is the best solution so far..
				SaveCurResult (dCurObj, j >= m_dwKSteps) ;
		}

		return TRUE ;
	}

	TClust::TClust (DWORD *pdwParIn, DWORD *pdwParOut, double *pdParIn, double *pdParOut, double *pdX, double *pdM, double *pdS, DWORD *pdwAssign, DWORD *pdwClustSize) :
				//	parameters
		m_n (pdwParIn[0]), m_p (pdwParIn[1]), m_K (pdwParIn[2]), m_dwIter (pdwParIn[3]), m_dwKSteps (pdwParIn[4]), m_dwEqualWeights (pdwParIn[5]), m_dwRestr (pdwParIn[6]), m_dwTrace (pdwParIn[7]), 
		m_dwConvCount (pdwParOut[0] = 0), m_dwIterSuccess (pdwParOut[1] = 0), m_dwCode (pdwParOut[2]),

		m_dAlpha (pdParIn[0]), m_dFactor (pdParIn[1]), m_dZeroTol (pdParIn[2]),
		m_dBestObj (pdParOut[0]),

			//	vals 2b calculated
		m_dDensFact (pow (2 * M_PI, m_p / -2.0)), m_dwNoTrim ((DWORD) floor (m_n * (1-m_dAlpha))),m_dwTrim (m_n-m_dwNoTrim ), 

			//	matrix stuff...
		m_vInd (m_n), m_vIndBest (m_n, pdwAssign), m_vIndOld (m_n), m_vCurInd (m_p + 1), 

		m_vClustSizeBest (m_K, pdwClustSize), m_vRank (m_n), 
		m_tEVec (m_p, m_p, m_K), m_tCurS (m_p, m_p, m_K),
		m_tBestS (IIter<3> (m_p, m_p, m_K), pdS), m_tBestEVec (m_p, m_p, m_K), 
		m_mCurM	(m_p, m_K), m_mBestM (m_p, m_K, pdM), m_mX (m_n, m_p, pdX), m_mLL (m_n, m_K), m_mEVal (m_p, m_K),
		m_mTempNP1 (m_n, m_p), m_mTempNP2 (m_n, m_p), m_mBestEVal (m_p, m_K), m_vWeights (m_K), m_vClustSize (m_K),
		m_vDisc (m_n), m_vDiscSorted (m_n), m_vTempN1(m_n), m_vTempN2 (m_n), m_vTempNPp2 (2 * m_n * m_p + 2)/*,  m_vTrimIdx(m_n), m_vTrimIdxBest(m_n)*/
	{
			//	splitting matrices and tensors into vector - and matrix arrays
		m_tEVec.GetIMatArray (0, 1, m_amEVec) ;
		m_tCurS.GetIMatArray (0, 1, m_amCurS) ;
		m_mCurM.GetIVecArray (0, m_avCurM) ;
		m_mEVal.GetIVecArray (0, m_avEVal) ;
		m_mLL.GetIVecArray	(0, m_avLL) ;

		tclust () ;
	}

	EXPORT void tclust (DWORD *pdwParIn, DWORD *pdwParOut, double *pdParIn, double *pdParOut, double *pdX, double *pdM, double *pdS, DWORD *pdwAssign, DWORD *pdwClustSize)
	{
		TClust (pdwParIn, pdwParOut, pdParIn, pdParOut, pdX, pdM, pdS, pdwAssign, pdwClustSize) ;
	}

	EXPORT void RestrictEigenValues (DWORD *pdwPar, double *pdPar, double *pdEV, double *pdwClustSize)
	{
		IMatD mEV (pdwPar[0], pdwPar[1], pdEV) ;
		IVecD vClustSize (pdwPar[1], pdwClustSize) ;

		IVecD vTemp ;

		RestrictEigenValues (mEV, vClustSize, pdPar[0], vTemp, pdPar[1]) ;
	}

	EXPORT void RestrictEigenValues_deter (DWORD *pdwPar, double *pdPar, double *pdEV, double *pdwClustSize)
	{
		IMatD mEV (pdwPar[0], pdwPar[1], pdEV) ;
		IVecD vClustSize (pdwPar[1], pdwClustSize) ;

		IVecD vTemp ;

		pdwPar [2] = RestrictEigenValues_deter (mEV, vClustSize, pdPar[0], vTemp, pdPar[1]) ;
	}
