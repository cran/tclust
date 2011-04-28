#include "tkmeans.h"

	CTKMeans::CTKMeans (t_size n, t_size p, t_size k, double dAlpha, double dZeroTol, double *pdX, int *pnAssign, double *pdClustSize, double *pdClustWeights, int nEqualWeights, int nTrace, double *pdM)
		: CClust (n, p, k, dAlpha, dZeroTol, pdX, pnAssign, pdClustSize, pdClustWeights, nEqualWeights, nTrace)
		, CClust_CM (pdM)
	{
		
	}

	void CTKMeans::CalcDensity (const SCMatD &mX, const SVecD &vDens, t_size k, const double dFact)
	{
		const SVecD &vCurM = m_mCurM.GetColRef (k) ;

		vDens.Reset (0) ;

		EO<UOP::Apa_sqr_BsC>::VMcVct_NC (*vDens, mX, vCurM) ;				//	vDens <- rowSums (sqr (mX - vCurM))

		EO<SOP::a_neg>::V (*vDens) ;										//	vDens <- -vDens

//		EO<UOP::Aa_Bm_exp_Adm2>::VSc (*vDens, dFact * m_dDensFact) ;
	}

	double CTKMeans::CalcObjFunc ()
	{
		double dObj = 0 ;

		t_size k ;

		SVecD vDensity (m_aTemp[3], 0) ;
		SMatD mCurX (m_aTemp[4], m_n, m_p) ;

//		double *pdClustSize = m_vClustSize ;

		for (k = m_K - 1; k != NAI; --k)
		//for (k = 0; k < m_K; ++k)
		{
			LoadCluster (mCurX, k) ;				//	2do: is this necessary? only summing up the according cells would be sufficient.. -> it doesn't cost too much..

			t_size dwClustSize = mCurX.nrow () ;

			if (!dwClustSize)
				continue ;

			vDensity.Require (dwClustSize) ;

			CalcDensity (mCurX, vDensity, k) ;

			EO<SOP::a_plus>::SVc (dObj, vDensity) ;

			//EO<UOP::Apa_logB>::SVc (dObj, vDensity) ;
		}

		return dObj ;
	}