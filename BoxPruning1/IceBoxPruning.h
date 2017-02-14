///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for box pruning.
 *	\file		IceBoxPruning.h
 *	\author		Pierre Terdiman
 *	\date		January, 29, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef __ICEBOXPRUNING_H__
#define __ICEBOXPRUNING_H__

	struct MESHMERIZER_API Axes
	{
		udword	Axis0;
		udword	Axis1;
		udword	Axis2;
	};

	// Optimized versions
	FUNCTION MESHMERIZER_API bool CompleteBoxPruning(udword nb, const AABB** list, Container& pairs, const Axes& axes);
	FUNCTION MESHMERIZER_API bool BipartiteBoxPruning(udword nb0, const AABB** list0, udword nb1, const AABB** list1, Container& pairs, const Axes& axes);

#endif // __ICEBOXPRUNING_H__
