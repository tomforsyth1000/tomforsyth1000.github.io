

// Cone creation/intersection library.
// Written by Tom Forsyth, somewhere in 2004.

// No warranties implied, expressed, impressed or explied.
// So don't go using this stuff for life support machines, now, you hear?
// But otherwise, do whatever you like with it.

// Sorry about the slightly odd math library. It's inherited from StarTopia,
// and is rather more "organic" than I'd ideally like. But it works.
// class M31 is a vector of three floats. I think everything else should be fairly obviously named.



// A view cone originates at the origin, has its axis along vAxis (always normalised). All points on it are defined by
// (vAxis * P)/P.length() = fCosAngle.
struct ViewCone
{
	M31 vAxis;
	float fCosAngle;

	void Check ( void ) const
	{
		ASSERT ( fabsf ( vAxis.LengthSquared() - 1.0f ) < 0.0001f );
	}
};


// Make a cone from the union of two unit-length vectors
ViewCone ConeMake ( const M31 &v1, const M31 &v2 )
{
	ASSERT ( fabsf ( v1.LengthSquared() - 1.0f ) < 0.01f );
	ASSERT ( fabsf ( v2.LengthSquared() - 1.0f ) < 0.01f );
	ViewCone res;
	res.vAxis = ( v1 + v2 ).Normalise();
	res.fCosAngle = DOT_PRODUCT ( v1, res.vAxis );
	return res;
}

// Make a cone from a bounding sphere.
ViewCone ConeMake ( const M31 &vCentre, float fRadius )
{
	ViewCone res;
	float fOneOverLengthSq = 1.0f / vCentre.LengthSquared();
	float fOneOverLength = sqrtf ( fOneOverLengthSq );
	res.vAxis = vCentre * fOneOverLength;

	if ( 1.0f <= ( fRadius * fOneOverLength ) )
	{
		// Point is inside the sphere.
		res.fCosAngle = -1.0f;
		return res;
	}

#if 0
	res.fCosAngle = cosf ( asin ( fRadius / fLength ) );
#else
	// AdjacentLength = sqrt ( fLength * fLength - fRadius * fRadius )
	// CosAngle = AdjacentLength / fLength
	//  = sqrt ( fLength * fLength - fRadius * fRadius ) / fLength
	//  = sqrt ( 1.0f - ( fRadius * fRadius ) / ( fLength * fLength ) )
	res.fCosAngle = sqrtf ( 1.0f - ( fRadius * fRadius ) * fOneOverLengthSq );
	ASSERT ( fabsf ( res.fCosAngle - cosf ( asinf ( fRadius * fOneOverLength ) ) ) < 0.00001f );
#endif
	return res;
}



bool ConeIntersect ( const ViewCone &vc1, const ViewCone &vc2 )
{
	if ( ( vc1.fCosAngle < 0.0f ) || ( vc2.fCosAngle < 0.0f ) )
	{
		return true;
	}

	vc1.Check();
	vc2.Check();

	// Ugh - lots of acosfs. Can't see a better way to do it.
	float fThetaTotal = acosf ( DOT_PRODUCT ( vc1.vAxis, vc2.vAxis ) );
	float fTheta1 = acosf ( vc1.fCosAngle );
	float fTheta2 = acosf ( vc2.fCosAngle );

	return ( fTheta1 + fTheta2 > fThetaTotal );
}



enum eConeUnionResult
{
	CUR_NORMAL,			// A new bounding cone was found.
	CUR_1ENCLOSES2,		// Cone 1 encloses cone 2 already.
	CUR_2ENCLOSES1,		// Cone 2 encloses cone 1 already.
	CUR_NOBOUND,		// There is no way to union these - the cone would subtend more than 180 degrees.
};

// Returns the cone that bounds both the input cones.
eConeUnionResult ConeUnion ( ViewCone *presult, const ViewCone &cone1, const ViewCone &cone2, bool bAlwaysSetResult = false )
{
	ASSERT ( presult != NULL );

	const float fVeryCloseEnough = 0.00001f;
	const float fCloseEnough = 0.001f;

	// Just check if they share axis.
	float fOneDotTwo = DOT_PRODUCT ( cone2.vAxis, cone1.vAxis );
	if ( fOneDotTwo > ( 1.0f - fVeryCloseEnough ) )
	{
		// Yep. OK, the test is really simple - which is bigger?
		if ( cone1.fCosAngle < cone2.fCosAngle )
		{
			if ( bAlwaysSetResult )
			{
				presult->vAxis = cone1.vAxis;
				presult->fCosAngle = cone1.fCosAngle - fVeryCloseEnough;
			}
			return CUR_1ENCLOSES2;
		}
		else
		{
			if ( bAlwaysSetResult )
			{
				presult->vAxis = cone2.vAxis;
				presult->fCosAngle = cone2.fCosAngle - fVeryCloseEnough;
			}
			return CUR_2ENCLOSES1;
		}
	}
	else if ( fOneDotTwo < ( -1.0f + fVeryCloseEnough ) )
	{
		// They point in completely opposite directions.
		if ( bAlwaysSetResult )
		{
			presult->fCosAngle = -1.0f;
		}
		return CUR_NOBOUND;
	}


	// Find the plane that includes both axis - this is the plane that the final cone axis will lie on as well.
	M31 vPlaneNormal = ( CROSS_PRODUCT ( cone2.vAxis, cone1.vAxis ) ).Normalise();

	// Now for each cone, find the "outer vector", which is the vector along the cone that lies in the plane,
	// furthest from the other cone's axis.
	// So define the vector vP = ( vAxis ^ vPlaneNormal ).norm()
	// So vP.vAxis = 0 and vP.vP=1.
	// Define:
	// vOuter = vAxis + lambda * ( ( vAxis ^ vPlaneNormal ).norm() )
	// and also:
	// ( vOuter * vAxis ) / vOuter.length() = fCosAngle
	// thus:
	// lambda = +/- sqrt ( ( 1 - fCosAngle^2 ) / ( fCosAngle^2 ) )
	//
	// For cone1, use +ve lambda, for cone2, use -ve.

	M31 vP1 = ( CROSS_PRODUCT ( vPlaneNormal, cone1.vAxis ) ).Normalise();
	float fCosAngleSquared1 = cone1.fCosAngle * cone1.fCosAngle;
	M31 vOuter1 = cone1.vAxis + vP1 * sqrtf ( ( 1.0f - fCosAngleSquared1 ) / fCosAngleSquared1 );
	vOuter1 = vOuter1.Normalise();

	M31 vP2 = ( CROSS_PRODUCT ( vPlaneNormal, cone2.vAxis ) ).Normalise();
	float fCosAngleSquared2 = cone2.fCosAngle * cone2.fCosAngle;
	M31 vOuter2 = cone2.vAxis - vP2 * sqrtf ( ( 1.0f - fCosAngleSquared2 ) / fCosAngleSquared2 );
	vOuter2 = vOuter2.Normalise();

	// Check to see if either outer vector is actually inside the other cone.
	// If it is, then that cone completely encloses the other.
#ifdef _DEBUG
	float fDebug1 = DOT_PRODUCT ( vOuter2, cone1.vAxis );
#endif
	if ( DOT_PRODUCT ( vOuter2, cone1.vAxis ) + fVeryCloseEnough >= cone1.fCosAngle )
	{
		// Obviously this means cone1 should be the wider one.
		ASSERT ( cone1.fCosAngle <= cone2.fCosAngle + fCloseEnough );
		// And cone2's axis should be inside it as well.
		ASSERT ( DOT_PRODUCT ( cone1.vAxis, cone2.vAxis ) + fCloseEnough >= cone1.fCosAngle );
		if ( bAlwaysSetResult )
		{
			*presult = cone1;
		}
		return CUR_1ENCLOSES2;
	}
#ifdef _DEBUG
	float fDebug2 = DOT_PRODUCT ( vOuter2, cone1.vAxis );
#endif
	if ( DOT_PRODUCT ( vOuter1, cone2.vAxis ) + fVeryCloseEnough >= cone2.fCosAngle )
	{
		// Obviously this means cone1 should be the wider one.
		ASSERT ( cone2.fCosAngle <= cone1.fCosAngle + fCloseEnough );
		// And cone2's axis should be inside it as well.
		ASSERT ( DOT_PRODUCT ( cone2.vAxis, cone1.vAxis ) + fCloseEnough >= cone2.fCosAngle );
		if ( bAlwaysSetResult )
		{
			*presult = cone2;
		}
		return CUR_2ENCLOSES1;
	}

	// And if the cross-product of the two points the opposite direction to vPlaneNormal, then the bounding
	// cone has turned "inside out" and covers an angle more than 180 degrees.
	if ( DOT_PRODUCT ( CROSS_PRODUCT ( vOuter1, vOuter2 ), vPlaneNormal ) >= 0.0f )
	{
		if ( bAlwaysSetResult )
		{
			presult->fCosAngle = -1.0f;
		}
		return CUR_NOBOUND;
	}

#ifdef _DEBUG
	// Taking copies coz presult might be one of the sources.
	float fCosAngle1 = cone1.fCosAngle;
	float fCosAngle2 = cone2.fCosAngle;
#endif

	// And now find the average of those two - that is the centre of the new cone.
	presult->vAxis = ( vOuter1 + vOuter2 ).Normalise();
	// And the angle.
	presult->fCosAngle = DOT_PRODUCT ( presult->vAxis, vOuter1 );
	// Of course it shouldn't matter which you use to find the angle.
#ifdef _DEBUG
	float f = DOT_PRODUCT ( presult->vAxis, vOuter2 );
	// The smaller fCosAngle, the bigger the allowable error.
	float fAllowedError = ( 1.0f - fabsf ( presult->fCosAngle ) ) * 0.2f + 0.001f;
	ASSERT ( fabsf ( presult->fCosAngle - f ) < fAllowedError );
#endif

	// And obviously the resulting cone can't be narrower than either of the two sources.
	// remember that narrower = higher value)
	ASSERT ( presult->fCosAngle <= fCosAngle1 + fCloseEnough );
	ASSERT ( presult->fCosAngle <= fCosAngle2 + fCloseEnough );

	// All done.
	return CUR_NORMAL;
}



// This is a version of ConeUnion, mangled so that all it does is decide if one cone encloses the other.
bool DoesCone2EncloseCone1 ( const ViewCone &cone1, const ViewCone &cone2 )
{
	const float fVeryCloseEnough = 0.00001f;
	const float fCloseEnough = 0.001f;

#ifdef _DEBUG
	ViewCone DebugResCone;
	eConeUnionResult DebugRes = ConeUnion ( &DebugResCone, cone1, cone2, false );
#endif

	// Just check if they share axis.
	float fOneDotTwo = DOT_PRODUCT ( cone2.vAxis, cone1.vAxis );
	if ( fOneDotTwo > ( 1.0f - fVeryCloseEnough ) )
	{
		// Yep. OK, the test is really simple - which is bigger?
		if ( cone1.fCosAngle < cone2.fCosAngle )
		{
			//return CUR_1ENCLOSES2;
			ASSERT ( DebugRes == CUR_1ENCLOSES2 );
			return false;
		}
		else
		{
			//return CUR_2ENCLOSES1;
			ASSERT ( DebugRes == CUR_2ENCLOSES1 );
			return true;
		}
	}
	else if ( fOneDotTwo < ( -1.0f + fVeryCloseEnough ) )
	{
		// They point in completely opposite directions.
		//return CUR_NOBOUND;
		ASSERT ( DebugRes == CUR_NOBOUND );
		return false;
	}


	// Find the plane that includes both axis - this is the plane that the final cone axis will lie on as well.
	M31 vPlaneNormal = ( CROSS_PRODUCT ( cone2.vAxis, cone1.vAxis ) ).Normalise();

	// Now for each cone, find the "outer vector", which is the vector along the cone that lies in the plane,
	// furthest from the other cone's axis.
	// So define the vector vP = ( vAxis ^ vPlaneNormal ).norm()
	// So vP.vAxis = 0 and vP.vP=1.
	// Define:
	// vOuter = vAxis + lambda * ( ( vAxis ^ vPlaneNormal ).norm() )
	// and also:
	// ( vOuter * vAxis ) / vOuter.length() = fCosAngle
	// thus:
	// lambda = +/- sqrt ( ( 1 - fCosAngle^2 ) / ( fCosAngle^2 ) )
	//
	// For cone1, use +ve lambda, for cone2, use -ve.

	M31 vP1 = ( CROSS_PRODUCT ( vPlaneNormal, cone1.vAxis ) ).Normalise();
	float fCosAngleSquared1 = cone1.fCosAngle * cone1.fCosAngle;
	M31 vOuter1 = cone1.vAxis + vP1 * sqrtf ( ( 1.0f - fCosAngleSquared1 ) / fCosAngleSquared1 );
	vOuter1 = vOuter1.Normalise();

#ifdef _DEBUG
	M31 vP2 = ( CROSS_PRODUCT ( vPlaneNormal, cone2.vAxis ) ).Normalise();
	float fCosAngleSquared2 = cone2.fCosAngle * cone2.fCosAngle;
	M31 vOuter2 = cone2.vAxis - vP2 * sqrtf ( ( 1.0f - fCosAngleSquared2 ) / fCosAngleSquared2 );
	vOuter2 = vOuter2.Normalise();
#endif

	// Check to see if either outer vector is actually inside the other cone.
	// If it is, then that cone completely encloses the other.
#ifdef _DEBUG
	float fDebug2 = DOT_PRODUCT ( vOuter2, cone1.vAxis );
#endif
	if ( DOT_PRODUCT ( vOuter1, cone2.vAxis ) + fVeryCloseEnough >= cone2.fCosAngle )
	{
		// Obviously this means cone1 should be the wider one.
		ASSERT ( cone2.fCosAngle <= cone1.fCosAngle + fCloseEnough );
		// And cone2's axis should be inside it as well.
		ASSERT ( DOT_PRODUCT ( cone2.vAxis, cone1.vAxis ) + fCloseEnough >= cone2.fCosAngle );
		//return CUR_2ENCLOSES1;
		ASSERT ( DebugRes == CUR_2ENCLOSES1 );
		return true;
	}

#if 0
	float fDebug1 = DOT_PRODUCT ( vOuter2, cone1.vAxis );
	if ( DOT_PRODUCT ( vOuter2, cone1.vAxis ) + fVeryCloseEnough >= cone1.fCosAngle )
	{
		// Obviously this means cone1 should be the wider one.
		ASSERT ( cone1.fCosAngle <= cone2.fCosAngle + fCloseEnough );
		// And cone2's axis should be inside it as well.
		ASSERT ( DOT_PRODUCT ( cone1.vAxis, cone2.vAxis ) + fCloseEnough >= cone1.fCosAngle );
		//return CUR_1ENCLOSES2;
		return false;
	}

	// And if the cross-product of the two points the opposite direction to vPlaneNormal, then the bounding
	// cone has turned "inside out" and covers an angle more than 180 degrees.
	if ( DOT_PRODUCT ( CROSS_PRODUCT ( vOuter1, vOuter2 ), vPlaneNormal ) >= 0.0f )
	{
		//return CUR_NOBOUND;
		return false;
	}
#endif

	//return CUR_NORMAL;

	ASSERT ( DebugRes != CUR_2ENCLOSES1 );

	return false;
}



float FindTanFromCos ( float fVal )
{
	float fTanTheta = ( sqrtf ( 1 - fVal * fVal ) ) / fVal;
	ASSERT ( fabsf ( fTanTheta - tanf(acosf(fVal)) ) < 0.001f );
	return fTanTheta;
}


// Calculate the size of a texture required for a given view cone and its pixel density.
float ConeSizeOfTexture ( const ViewCone &vc, float fPixelDensity )
{
	float fTanTheta = FindTanFromCos ( vc.fCosAngle );
	return fPixelDensity * fTanTheta;
}


