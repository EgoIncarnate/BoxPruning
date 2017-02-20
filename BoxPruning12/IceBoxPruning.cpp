///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for the "box pruning revisited" project.
 *	\file		IceBoxPruning.cpp
 *	\author		Pierre Terdiman
 *	\date		February 2017
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Precompiled Header
#include "Stdafx.h"

using namespace Meshmerizer;

// InsertionSort has better coherence, RadixSort is better for one-shot queries.
#define PRUNING_SORTER	RadixSort
//#define PRUNING_SORTER	InsertionSort

static __forceinline int intersects2D(const AABB& a, const AABB& b)
{
       if(    b.mMax.y < a.mMin.y || a.mMax.y < b.mMin.y
       ||     b.mMax.z < a.mMin.z || a.mMax.z < b.mMin.z)
              return 0;
       return 1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Bipartite box pruning. Returns a list of overlapping pairs of boxes, each box of the pair belongs to a different set.
 *	\param		nb0		[in] number of boxes in the first set
 *	\param		list0	[in] list of boxes for the first set
 *	\param		nb1		[in] number of boxes in the second set
 *	\param		list1	[in] list of boxes for the second set
 *	\param		pairs	[out] list of overlapping pairs
 *	\param		axes	[in] projection order (0,2,1 is often best)
 *	\return		true if success.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Meshmerizer::BipartiteBoxPruning(udword nb0, const AABB* list0, udword nb1, const AABB* list1, Container& pairs)
{
	// Checkings
	if(!nb0 || !list0 || !nb1 || !list1)
		return false;

	AABB* BoxList0 = new AABB[nb0+1];
	AABB* BoxList1 = new AABB[nb1+1];
	udword* Remap0;
	udword* Remap1;
	{
		// Allocate some temporary data
		float* PosList0 = new float[nb0+1];
		float* PosList1 = new float[nb1+1];

		// 1) Build main lists using the primary axis
		for(udword i=0;i<nb0;i++)
			PosList0[i] = list0[i].mMin.x;
		PosList0[nb0] = FLT_MAX;
		for(udword i=0;i<nb1;i++)
			PosList1[i] = list1[i].mMin.x;
		PosList1[nb1] = FLT_MAX;

		// 2) Sort the lists
		static PRUNING_SORTER RS0, RS1;	// Static for coherence.
		Remap0 = RS0.Sort(PosList0, nb0+1).GetRanks();
		Remap1 = RS1.Sort(PosList1, nb1+1).GetRanks();

		for(udword i=0;i<nb0;i++)
			BoxList0[i] = list0[Remap0[i]];
		BoxList0[nb0].mMin.x = FLT_MAX;

		for(udword i=0;i<nb1;i++)
			BoxList1[i] = list1[Remap1[i]];
		BoxList1[nb1].mMin.x = FLT_MAX;

		DELETEARRAY(PosList1);
		DELETEARRAY(PosList0);
	}

	// 3) Prune the lists
	udword Index0 = 0;
	udword RunningAddress1 = 0;
	while(RunningAddress1<nb1 && Index0<nb0)
	{
		const AABB& Box0 = BoxList0[Index0];

		const float MinLimit = Box0.mMin.x;
		while(BoxList1[RunningAddress1].mMin.x<MinLimit)
			RunningAddress1++;

		const udword RIndex0 = Remap0[Index0];
		udword Index1 = RunningAddress1;
		const float MaxLimit = Box0.mMax.x;
		while(BoxList1[Index1].mMin.x<=MaxLimit)
		{
			if(intersects2D(Box0, BoxList1[Index1]))
				pairs.Add(RIndex0).Add(Remap1[Index1]);

			Index1++;
		}
		Index0++;
	}

	////

	Index0 = 0;
	udword RunningAddress0 = 0;
	while(RunningAddress0<nb0 && Index0<nb1)
	{
		const AABB& Box1 = BoxList1[Index0];

		const float MinLimit = Box1.mMin.x;
		while(BoxList0[RunningAddress0].mMin.x<=MinLimit)
			RunningAddress0++;

		const udword RIndex1 = Remap1[Index0];
		udword Index1 = RunningAddress0;
		const float MaxLimit = Box1.mMax.x;
		while(BoxList0[Index1].mMin.x<=MaxLimit)
		{
			if(intersects2D(BoxList0[Index1], Box1))
				pairs.Add(Remap0[Index1]).Add(RIndex1);

			Index1++;
		}
		Index0++;
	}

	DELETEARRAY(BoxList1);
	DELETEARRAY(BoxList0);

	return true;
}


// Munge the float bits to return produce an unsigned order-preserving
// ranking of floating-point numbers.
// (Old trick: http://stereopsis.com/radix.html FloatFlip, with a new
// spin to get rid of -0.0f)

// In /fp:precise, we can just calc "x + 0.0f" and get what we need.
// But fast math optimizes it away. Could use #pragma float_control,
// but that prohibits inlining of MungeFloat. So do this silly thing
// instead.
float g_global_this_always_zero = 0.0f;

union FloatOrInt32
{
	float f;
	sdword s;
};

static inline udword MungeFloat(float f)
{
	FloatOrInt32 u;
    u.f = f + g_global_this_always_zero;  // NOT a nop! Canonicalizes -0.0f to +0.0f
    udword toggle = (u.s >> 31) & ~(1u << 31);
    return u.s ^ toggle;
}
//#pragma float_control(pop)

struct SIMD_AABB_X
{
	__forceinline	SIMD_AABB_X()	{}
	__forceinline	~SIMD_AABB_X()	{}

	void	InitFrom(const AABB& b)
	{
		mMinX	= MungeFloat(b.mMin.x);
		mMaxX	= MungeFloat(b.mMax.x);
	}

    sdword mMinX;
    sdword mMaxX;
};

struct SIMD_AABB_YZ
{
	__forceinline	SIMD_AABB_YZ()	{}
	__forceinline	~SIMD_AABB_YZ()	{}

	void	InitFrom(const AABB& b)
	{
		mMinY	= b.mMin.y;
		mMinZ	= b.mMin.z;
		mMaxY	= b.mMax.y;
		mMaxZ	= b.mMax.z;
	}

    float mMinY;
    float mMinZ;
    float mMaxY;
    float mMaxZ;
};

/*static __forceinline int intersects2D(const SIMD_AABB_YZ& a, const SIMD_AABB_YZ& b)
{
       if(    b.mMaxY <= a.mMinY || a.mMaxY < b.mMinY
       ||     b.mMaxZ <= a.mMinZ || a.mMaxZ < b.mMinZ)
              return 0;
       return 1;
}*/

	#include <xmmintrin.h>
	#include <emmintrin.h>

/*#define SIMD_OVERLAP_INIT(box)	\
       const __m128 b = _mm_shuffle_ps(_mm_load_ps(&box.mMinY), _mm_load_ps(&box.mMinY), 78);

#define SIMD_OVERLAP_TEST(box)						\
       const __m128 a = _mm_load_ps(&box.mMinY);	\
       const __m128 d = _mm_cmpgt_ps(a, b);			\
       if(_mm_movemask_ps(d)==12)*/

static void __cdecl outputPair(udword id0, udword id1, Container& pairs, const udword* remap)
{
	pairs.Add(remap[id0]).Add(remap[id1]);
}

static void __cdecl outputPair2(udword id0, udword id1, Container& pairs)
{
	pairs.Add(id0).Add(id1);
}

/*static void debug(udword id0, udword id1, Container& pairs, const udword* remap)
{
	if(id0==9 && remap[id1]==6)
		int FoundIt=0;
}*/

// Count trailing zeroes
static inline udword Ctz32(udword x)
{
	unsigned long idx;
	_BitScanForward(&idx, x);
	return idx;
}

// Reports a bunch of intersections as specified by a base index and a bit mask.
static void ReportIntersections(Container& pairs, udword remap_id0, const udword* remap, uword base_id1, udword mask)
{
	while (mask)
	{
		pairs.Add(remap_id0).Add(remap[base_id1 + Ctz32(mask)]);
		mask &= mask - 1;
	}
}

static void Error()
{
	printf("ERROR!\n");
}

template<typename T>
static inline T* PtrAddBytes(T* ptr, ptrdiff_t bytes)
{
	return (T*)((char*)ptr + bytes);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Complete box pruning. Returns a list of overlapping pairs of boxes, each box of the pair belongs to the same set.
 *	\param		nb		[in] number of boxes
 *	\param		list	[in] list of boxes
 *	\param		pairs	[out] list of overlapping pairs
 *	\param		axes	[in] projection order (0,2,1 is often best)
 *	\return		true if success.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Meshmerizer::CompleteBoxPruning(udword nb, const AABB* list, Container& pairs)
{
	// Checkings
	if(!nb || !list)
		return false;

	udword nbpad = nb+5;
	ptrdiff_t BoxBytesP = nbpad*sizeof(FloatOrInt32);
	ptrdiff_t BoxBytesN = -BoxBytesP;
	ptrdiff_t BoxBytes3N = 3*BoxBytesN;

	// BoxSOA: in order, arrays for MaxX,MinX (int), MaxY,MinY,MaxZ,MinZ (float).
	FloatOrInt32* BoxSOA = new FloatOrInt32[nbpad*6];

	// Our default origin is actually pointing at array number 3 (MinY).
	FloatOrInt32* BoxBase = PtrAddBytes(BoxSOA, 3*BoxBytesP);
	FloatOrInt32* BoxEnd = BoxBase + nb;

	udword* Remap;
	//{
		// Allocate some temporary data
		float* PosList = new float[nb+1];

		// 1) Build main list using the primary axis
		for(udword i=0;i<nb;i++)
			PosList[i] = list[i].mMin.x;
		PosList[nb] = FLT_MAX;

		// 2) Sort the list
		static PRUNING_SORTER RS;	// Static for coherence
		Remap = RS.Sort(PosList, nb+1).GetRanks();

		// 3) Prepare the SoA box array
		for(udword i=0;i<nb;i++)
		{
			const AABB& Box = list[Remap[i]];
			FloatOrInt32 *OutBoxI = &BoxBase[i];
			PtrAddBytes(OutBoxI,  BoxBytes3N)->s = MungeFloat(Box.mMax.x);
			PtrAddBytes(OutBoxI, 2*BoxBytesN)->s = MungeFloat(Box.mMin.x);
			PtrAddBytes(OutBoxI, 1*BoxBytesN)->f = Box.mMax.y;
			PtrAddBytes(OutBoxI, 0*BoxBytesP)->f = Box.mMin.y;
			PtrAddBytes(OutBoxI, 1*BoxBytesP)->f = Box.mMax.z;
			PtrAddBytes(OutBoxI, 2*BoxBytesP)->f = Box.mMin.z;
		}
		for(udword i=nb;i<nbpad;i++)
		{
			FloatOrInt32 *OutBoxI = &BoxBase[i];
			PtrAddBytes(OutBoxI,  BoxBytes3N)->s = -0x80000000;
			PtrAddBytes(OutBoxI, 2*BoxBytesN)->s = 0x7fffffff;
			PtrAddBytes(OutBoxI, 1*BoxBytesN)->f = -FLT_MAX;
			PtrAddBytes(OutBoxI, 0*BoxBytesP)->f = FLT_MAX;
			PtrAddBytes(OutBoxI, 1*BoxBytesP)->f = -FLT_MAX;
			PtrAddBytes(OutBoxI, 2*BoxBytesP)->f = FLT_MAX;
		}
		DELETEARRAY(PosList);
	//}

	// 4) Prune the list
	FloatOrInt32 *Box0Ptr = BoxBase; // corresponds to Index0
	FloatOrInt32 *RunningPtr = BoxBase; // corresponds to RunningAddress
	while(Box0Ptr < BoxEnd && RunningPtr < BoxEnd)
	{
		const sdword MinLimit = PtrAddBytes(Box0Ptr, 2*BoxBytesN)->s;
		while (PtrAddBytes(RunningPtr++, 2*BoxBytesN)->s < MinLimit);

		const FloatOrInt32 *MaxLimitPtr = PtrAddBytes(PtrAddBytes(Box0Ptr, 2*BoxBytesN), BoxBytesN);
		const sdword MaxLimit = MaxLimitPtr->s;
		udword RIndex0 = Remap[Box0Ptr - BoxBase];

		__m128i MaxLimitVec = _mm_set1_epi32(MaxLimitPtr->s);
		__m128 Box0MaxY = _mm_set1_ps(PtrAddBytes(Box0Ptr, 1*BoxBytesN)->f);
		__m128 Box0MinY = _mm_set1_ps(PtrAddBytes(Box0Ptr, 0*BoxBytesP)->f);
		__m128 Box0MaxZ = _mm_set1_ps(PtrAddBytes(Box0Ptr, 1*BoxBytesP)->f);
		__m128 Box0MinZ = _mm_set1_ps(PtrAddBytes(Box0Ptr, 2*BoxBytesP)->f);

		// Main loop
		const FloatOrInt32 *Box1Ptr = RunningPtr;
		while (PtrAddBytes(Box1Ptr + 3, 2*BoxBytesN)->s <= MaxLimit) // while Box[Index1+3].mMinX <= MaxLimit
		{
			// Load 4 boxes worth
			__m128 Box1MaxY = _mm_loadu_ps(&PtrAddBytes(Box1Ptr, 1*BoxBytesN)->f);
			__m128 Box1MinY = _mm_loadu_ps(&PtrAddBytes(Box1Ptr, 0*BoxBytesP)->f);
			__m128 Box1MaxZ = _mm_loadu_ps(&PtrAddBytes(Box1Ptr, 1*BoxBytesP)->f);
			__m128 Box1MinZ = _mm_loadu_ps(&PtrAddBytes(Box1Ptr, 2*BoxBytesP)->f);
			Box1Ptr += 4;

			// Intersection test:
			//  !(b.MaxY <= a.MinY) && (b.MinY < a.MaxY) && !(b.MaxZ <= a.MinZ) && (b.MinZ <= a.MaxZ)
			__m128 Cmp;
			Cmp = _mm_cmpnle_ps(Box1MaxY, Box0MinY);
			Cmp = _mm_and_ps(Cmp, _mm_cmplt_ps(Box1MinY, Box0MaxY));
			Cmp = _mm_and_ps(Cmp, _mm_cmpnle_ps(Box1MaxZ, Box0MinZ));
			Cmp = _mm_and_ps(Cmp, _mm_cmplt_ps(Box1MinZ, Box0MaxZ));

			int Mask = _mm_movemask_ps(Cmp);
			if (Mask)
				ReportIntersections(pairs, RIndex0, Remap, (Box1Ptr - 4) - BoxBase, Mask);
		}

		// Tail group: first box is in, but one or more boxes with mMinX past MaxLimit inside.
		if (PtrAddBytes(Box1Ptr, 2*BoxBytesN)->s <= MaxLimit)
		{
			__m128i Box1MinX = _mm_loadu_si128((const __m128i *)&PtrAddBytes(Box1Ptr, 2*BoxBytesN)->s);
			__m128 OutsideMask = _mm_castsi128_ps(_mm_cmpgt_epi32(Box1MinX, MaxLimitVec));

			// Load 4 boxes worth
			__m128 Box1MaxY = _mm_loadu_ps(&PtrAddBytes(Box1Ptr, 1*BoxBytesN)->f);
			__m128 Box1MinY = _mm_loadu_ps(&PtrAddBytes(Box1Ptr, 0*BoxBytesP)->f);
			__m128 Box1MaxZ = _mm_loadu_ps(&PtrAddBytes(Box1Ptr, 1*BoxBytesP)->f);
			__m128 Box1MinZ = _mm_loadu_ps(&PtrAddBytes(Box1Ptr, 2*BoxBytesP)->f);

			// Intersection test:
			//  !(b.MaxY <= a.MinY) && (b.MinY < a.MaxY) && !(b.MaxZ <= a.MinZ) && (b.MinZ <= a.MaxZ)
			__m128 Cmp;
			Cmp = _mm_andnot_ps(OutsideMask, _mm_cmpnle_ps(Box1MaxY, Box0MinY));
			Cmp = _mm_and_ps(Cmp, _mm_cmplt_ps(Box1MinY, Box0MaxY));
			Cmp = _mm_and_ps(Cmp, _mm_cmpnle_ps(Box1MaxZ, Box0MinZ));
			Cmp = _mm_and_ps(Cmp, _mm_cmplt_ps(Box1MinZ, Box0MaxZ));

			int Mask = _mm_movemask_ps(Cmp);
			if (Mask)
				ReportIntersections(pairs, RIndex0, Remap, Box1Ptr - BoxBase, Mask);
		}
		Box0Ptr++;
	}

	DELETEARRAY(BoxSOA);
	return true;
}

