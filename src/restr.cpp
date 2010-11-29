#include "restr.h"

	BOOL CheckClusterSingularity (const IMatD &mEV, const IVecD & vClustSize, double dZeroTol) ;

	double CalcDiff_log (const IMatD &ev, const IVecD &ks, double dLower, double dUpper)
	{
		if (dLower > dUpper)
			swap (&dLower, &dUpper) ;

		double dRetC = 0, dRet = 0 ;

		DWORD rm1 = ev.nrow () - 1 ;
		DWORD c, r ;

		for (c = ev.ncol () - 1; c != (DWORD) -1; c--)	//	iterates the clusters (coloumns)
		{
			dRetC = 0 ;
			for (r = rm1; r != (DWORD) -1; r--)			//	iterates the EVs of the current cluster
			{
				const double &dCur = ev(r, c) ;

				if (dCur < dLower)
					dRetC += log (dLower) + dCur / dLower ;
				else if (dCur > dUpper)
					dRetC += log (dUpper) + dCur / dUpper ;
				else
					dRetC += log (dCur) + 1 ;
			}

			dRet += ks(c) * dRetC ;
		}

		return dRet ;
	}

	void CalcRST (const IVecD &vec, double dL, double dU, DWORD &dwR, double &dS, double &dT)
	{
		DWORD c;
		dwR = 0 ;
		dS = dT = 0 ;

		for (c = vec.size () - 1; c != (DWORD) -1; c--)
		{
			double &dCur = vec (c) ;

			if (dCur < dL)
				dS += dCur ;
			if (dCur > dU)
				dT += dCur ;
			if (dCur < dL || dCur > dU)
				dwR ++ ;
		}
	}

	BOOL ZeroGroupsGetMeanEigenvalues (const IMatD &mEV, const IVecD & vClustSize)
	{
			// --> just care about groups without observations. set their eigenvalues to the mean of all other EVs

		DWORD dwK = vClustSize.size () ;
		DWORD k, dwCount = 0 ;
		double dSum = 0 ;

		for (k = dwK - 1; k != (DWORD) -1; k--)
			if (vClustSize (k))
			{
				dwCount ++ ;
				dSum += mEV.GetCol (k).sum () ;
			}

		dSum /= dwCount * mEV.nrow () ;

		for (k = dwK - 1; k != (DWORD) -1; k--)
			if (!vClustSize (k))
				mEV.GetCol (k).Reset (dSum) ;
		return TRUE ;
	}

	void GetCheckArray (const IMatD &mEV, const IVecD & vClustSize, double dFact, IVecD &vdCheckEV, double dMax)
	{
		vdCheckEV.Reshape (mEV.size () * 2 + 2) ;

		vdCheckEV(0) = 0 ;
		vdCheckEV(1) = dMax ;

		DWORD dwCur = 2 ;

		DWORD rm1 = mEV.nrow () - 1 ;
		DWORD c, r ;

		for (c = mEV.ncol () - 1; c != (DWORD) -1; c--)	//	iterates the clusters (coloumns)
			for (r = rm1; r != (DWORD) -1; r--)			//	iterates the EVs of the current cluster
			{
				const double &dCur = mEV(r, c) ;
				vdCheckEV (dwCur++) = dCur ;
				vdCheckEV (dwCur++) = dCur / dFact;
			}

		sort (vdCheckEV) ;
			//	calcs the average of all the neigbours of vdCheckEV
		for (c = 1; c < dwCur; c++)
			vdCheckEV(c - 1) = (vdCheckEV(c - 1) + vdCheckEV(c)) / 2 ;
		vdCheckEV.Reshape (--dwCur) ;
	}

	BOOL RestrictEigenValues (const IMatD &mEV, const IVecD & vClustSize, double dFact, IVecD &vTempNPp2, double dZeroTol)
	{
		if (!CheckClusterSingularity (mEV, vClustSize, dZeroTol))
			return FALSE ;	//	all eigenvalues of all clusters with at least one observation are 0

		if (dFact <0)
			dFact = -dFact ;

		if (dFact < 1)
			dFact = 1 / dFact ;

		double	dMin = 0, dMax = 0 ;
		DWORD i ;

//		mEV.MinMax (dMin, dMax) ;	//	bug: only clusters with vClustSize > 0 shall be considered here...

		BOOL bFoundOne = FALSE ;	//	XXXC 20101108	calculating the min / max of Clustersizes only for clusters with a size > 0
		for (i = mEV.ncol () - 1; i != (DWORD) -1; i--)
			if (vClustSize (i) > dZeroTol)
			{
				mEV.GetCol (i).MinMax (dMin, dMax, !bFoundOne) ;
				bFoundOne = TRUE ;
			}

		DWORD dwK = vClustSize.size () ;

		if (//dMin <= dZeroTol ||
			dMax / dMin > dFact)	//	min (ev) < max (ev) / dFact --> so we have to restrict the eigenvalues
		{
			IVecD &vdCheckEV = vTempNPp2 ;

			GetCheckArray (mEV, vClustSize, dFact, vdCheckEV, dMax) ;

			double dMinVal = 0, dMinM = 0 ;

			double dSumU = 0, dSumL = 0 ;
			
			DWORD c, r, dwCur = vdCheckEV.size () ;
			for (c = 0; c < dwCur; c++)
			//for (c = vdCheckEV.size () - 1; c != (DWORD) -1; c--)
			{
				dSumU = dSumL = 0 ;
				double &dCurL = vdCheckEV (c), dCurU = dCurL * dFact ;

				for (r = dwK - 1; r != (DWORD) -1; r--)
				{
					DWORD dwR ;
					double dS, dT ;
					CalcRST (mEV.GetCol (r), dCurL, dCurU,dwR, dS, dT) ; 

					dSumU += vClustSize (r) * (dS + dT/ dFact) ;
					dSumL += vClustSize (r) * dwR ;
				}

				dSumU /= dSumL ;

				double dCurVal = CalcDiff_log (mEV, vClustSize, dSumU, dSumU * dFact) ;

				if (!c || dMinVal > dCurVal)
				{
					dMinVal = dCurVal ;
					dMinM = dSumU ;
				}
			}

			mEV.Limit (dMinM, dMinM * dFact) ;
		}
		else	//	we don't have to restrict eigenvalues.
			ZeroGroupsGetMeanEigenvalues (mEV, vClustSize) ;

		return CheckClusterSingularity (mEV, vClustSize, dZeroTol) ;
	}

	BOOL CheckClusterSingularity (const IMatD &mEV, const IVecD & vClustSize, double dThreshold)
	{
		DWORD i, j;

		for (i = mEV.ncol () - 1; i != (DWORD) -1; i--)
			if (vClustSize (i) > 0)
				for (j = mEV.nrow () - 1; j != (DWORD) -1; j--)
					if (mEV (j, i) > dThreshold)
						return TRUE ;
		return FALSE ;
	}

	BOOL RestrictEigenValues_deter (const IMatD &mEV, const IVecD & vClustSize, double dFact, IVecD &vTempNPp2, double dZeroTol)
	{
		DWORD p = mEV.nrow (), k = mEV.ncol () ;

		IMatD mDeter (1, k) ;
		IVecD vDeter (mDeter.GetRow (0)), vMins (k);

		const double dpInv = 1.0/p ;

		vDeter.Reshape (k) ;
		colProds (mEV, vDeter) ;

		if (!CheckClusterSingularity (mDeter, vClustSize, dZeroTol))
			return FALSE ;

		DWORD i, j ;
						//////
						//	calculating the min of the determinant vector, setting all values < dZeroTol to dZeroTol
		for (i = k - 1; i != (DWORD) -1; i--)			//	for all coloumns (clusters)
		{
			IVecD curCol = mEV.GetCol (i) ;

			double &dCurMin = vMins (i) ;

			for (j = 0; j < p; j++)						//	for each eigenvalue ov the current cluster
			{
				double &dCurVal = curCol (j) ;

				if (dCurVal <= dZeroTol)				//	cur value is smaller than zero tol
					dCurVal = dCurMin = dZeroTol ;
				else if (!j || dCurVal < dCurMin)		//	or smaller than the current minimum (or the first value added)
					dCurMin = dCurVal ;
			}

			curCol.Limit_U (dCurMin / dZeroTol) ;		//	war dCurMin * 1e15
			curCol /= pow (prod (curCol), dpInv) ;
		}

		double dMin = 0, dMax = 0 ;
		vDeter.MinMax (dMin, dMax) ;

		mDeter ^= dpInv ;
		if (dMax / dMin <= dFact)
			ZeroGroupsGetMeanEigenvalues (mDeter, vClustSize) ;
		else
			RestrictEigenValues (mDeter, vClustSize, pow (dFact, dpInv), vTempNPp2, dZeroTol) ;

		mEV.byrow () *= vDeter ;

		return CheckClusterSingularity (mEV, vClustSize, dZeroTol) ;
	}
