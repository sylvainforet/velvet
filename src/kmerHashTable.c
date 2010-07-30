/*
Copyright 2009 Sylvain Foret (sylvain.foret@anu.edu.au) 

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

/**
 * Hash table with direct addressing and quadratic probing
 * Loosely based from Glib's GHash table
 * http://www.gtk.org
 * This version has `in-line' nodes, no tombstones and resizes in place
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "globals.h"
#include "kmerHashTable.h"
#include "kmer.h"
#include "utility.h"


#define HASH_TABLE_MIN_SHIFT 16  /* 1 << 3 == 8 buckets */

static const uint64_t primeMod[] =
{
  1ul,          /* For 1 << 0 */
  2ul,
  3ul,
  7ul,
  13ul,
  31ul,
  61ul,
  127ul,
  251ul,
  509ul,
  1021ul,
  2039ul,
  4093ul,
  8191ul,
  16381ul,
  32749ul,
  65521ul,      /* For 1 << 16 */
  131071ul,
  262139ul,
  524287ul,
  1048573ul,
  2097143ul,
  4194301ul,
  8388593ul,
  16777213ul,
  33554393ul,
  67108859ul,
  134217689ul,
  268435399ul,
  536870909ul,
  1073741789ul,
  2147483647ul,  /* For 1 << 31 */
  4294967291ul,
  8589934583ul,
  17179869143ul,
  34359738337ul,
  68719476731ul,
  137438953447ul,
  274877906899ul,
  549755813881ul,
  1099511627689ul
};

typedef struct kmerHashNode_st KmerHashNode;

struct kmerHashNode_st
{
	KmerKey keyHash;
	IDnum position;
	IDnum seqID;
	Kmer kmer;
} ATTRIBUTE_PACKED;

struct kmerHashTable_st
{
	KmerHashNode *nodes;

	uint64_t size;
	uint64_t mod;
	uint64_t mask;

	Coordinate nNodes;
};

static void
hashTableSetShift (KmerHashTable *hTable, int shift)
{
	int i;
	uint64_t mask = 0;

	hTable->size = 1 << shift;
	hTable->mod  = primeMod[shift];

	for (i = 0; i < shift; i++)
	{
		mask <<= 1;
		mask |= 1;
	}

	hTable->mask = mask;
}

static int
hashTableFindClosestShift (Coordinate n)
{
	int i;

	for (i = 0; n; i++)
		n >>= 1;

	return i;
}

/* Thomas Wang (64 bits) and Robert Juenkins (32 bits) bit mixing functions
 * I have a feeling that the 64 bit version is not quite right ...
 *
 * 1 is OR'ed to the key to denote that it is not an empty node
 */
static inline KmerKey
kmerHashFunc(Kmer * kmer)
{
#if KMER_LONGLONGS
	KmerKey key = kmer->longlongs[0];

/* These `#if' are ugly, but they should cater for the large majority of cases
 * and allow to avoid loops */
#if KMER_LONGLONGS > 1
	key ^= kmer->longlongs[1];
#endif
#if KMER_LONGLONGS > 2
	key ^= kmer->longlongs[2];
#endif

	key = (~key) + (key << 21);
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8);
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4);
	key = key ^ (key >> 28);
	key = key + (key << 31);

	return key | 1;
#elif KMER_LONGS
	KmerKey key = kmer->longs[0];

	key += ~(key << 15);
	key ^=  (key >> 10);
	key +=  (key << 3);
	key ^=  (key >> 6);
	key += ~(key << 11);
	key ^=  (key >> 16);

	return key | 1;

#elif KMER_INTS
	return kmer->ints | 1;
#elif KMER_CHARS
	return kmer->chars | 1;
#endif
}

static void
hashTableSetShiftFromSize (KmerHashTable *hTable, Coordinate size)
{
	int shift;

	shift = hashTableFindClosestShift (size);
	if (HASH_TABLE_MIN_SHIFT > shift)
		shift = HASH_TABLE_MIN_SHIFT;
	hashTableSetShift(hTable, shift);
}

static void
hashTableResize (KmerHashTable *hTable)
{
	const Coordinate oldSize = hTable->size;
	Coordinate i;
	boolean *touched;

#define IS_TOUCHED(i)  ((touched[(i) / 8] >> ((i) % 8)) & 1)
#define SET_TOUCHED(i) (touched[(i) / 8] |= (1 << ((i) % 8)))

	hashTableSetShiftFromSize(hTable, hTable->nNodes * 2);
	hTable->nodes = reallocOrExit(hTable->nodes, hTable->size, KmerHashNode);
	memset(hTable->nodes + oldSize, 0, (hTable->size - oldSize) * sizeof(KmerHashNode));
	touched = callocOrExit(oldSize / 8 + 1, boolean);

	for (i = 0; i < oldSize; i++)
	{
		KmerHashNode *node;
		KmerHashNode *nextNode;
		KmerHashNode tmpNode;
		KmerKey hashVal;
		Coordinate step;

		if (IS_TOUCHED(i))
			continue;
		node = &hTable->nodes[i];
		if (node->keyHash == 0)
			continue;
		hashVal = node->keyHash % hTable->mod;
		if (hashVal == i)
		{
			SET_TOUCHED(i);
			continue;
		}
		tmpNode = *node;
		nextNode = &hTable->nodes[hashVal];
		node->keyHash = 0;
		step = 0;
		while (nextNode->keyHash)
		{
			/* Hit an untouched element */
			if (hashVal < oldSize && !IS_TOUCHED(hashVal))
			{
				KmerHashNode tmpp;

				tmpp = tmpNode;
				tmpNode = *nextNode;
				*nextNode = tmpp;
				SET_TOUCHED(hashVal);
				hashVal = tmpNode.keyHash % hTable->mod;
				step = 0;
			}
			else
			{
				step++;
				hashVal += step;
				hashVal &= hTable->mask;
			}
			nextNode = &hTable->nodes[hashVal];
		}
		*nextNode = tmpNode;
	}
	free(touched);
#undef IS_TOUCHED
#undef SET_TOUCHED
}

static inline void hashTableMaybeResize (KmerHashTable *hTable)
{
	Coordinate nOccupied = hTable->nNodes;
	Coordinate size = hTable->size;

	if ((size > hTable->nNodes * 4 && size > 1 << HASH_TABLE_MIN_SHIFT) ||
			(size <= nOccupied + (nOccupied / 16)))
		hashTableResize(hTable);
}

KmerHashTable *newKmerHashTable(void)
{
	KmerHashTable *hTable;

	hTable = mallocOrExit(1, KmerHashTable);
	hTable->nNodes = 0;
	hashTableSetShift(hTable, HASH_TABLE_MIN_SHIFT);
	hTable->nodes = callocOrExit(hTable->size, KmerHashNode);

	return hTable;
}

void destroyKmerHashTable(KmerHashTable * hTable)
{
	if (hTable)
	{
		if (hTable->nodes)
			free(hTable->nodes);
		free(hTable);
	}
}

static inline KmerKey
hashTableLookupNodeForInsertion (KmerHashTable *hTable,
				 Kmer *kmer, KmerKey *hashReturn)
{
	KmerHashNode *node;
	Coordinate nodeIndex;
	Coordinate step = 0;
	KmerKey hashValue;

	hashValue = kmerHashFunc(kmer);

	*hashReturn = hashValue;
	nodeIndex = hashValue % hTable->mod;
	node = &hTable->nodes[nodeIndex];

	while (node->keyHash)
	{
		/*  We first check if our full hash values
		 *  are equal so we can avoid calling the full-blown
		 *  key equality function in most cases.
		 */
		if (node->keyHash == hashValue)
		{
			if (!compareKmers(&node->kmer, kmer))
				return nodeIndex;
		}
		step++;
		nodeIndex += step;
		nodeIndex &= hTable->mask;
		node = &hTable->nodes[nodeIndex];
	}

	return nodeIndex;
}
boolean findOrInsertOccurenceInKmerHashTable(KmerHashTable *hTable, Kmer *kmer,
					     IDnum *seqID, Coordinate *position)
{
	KmerHashNode *node;
	Coordinate nodeIndex;
	KmerKey keyHash;
	KmerKey oldHash;

	nodeIndex = hashTableLookupNodeForInsertion(hTable, kmer, &keyHash);
	node = &hTable->nodes[nodeIndex];
	oldHash = node->keyHash;

	if (oldHash > 1)
	{
		*seqID = node->seqID;
		*position = node->position;
		return true;
	}

	copyKmers(&node->kmer, kmer);
	node->keyHash = keyHash;
	hTable->nNodes++;

	node->seqID = *seqID;
	node->position = *position;
	hashTableMaybeResize(hTable);

	return false;
}
