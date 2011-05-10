/*
Copyright 2009 Sylvain Foret (sylvain.foret@anu.edu.au) 

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "hirschberg.h"
#include "recycleBin.h"
#include "tightString.h"
#include "utility.h"


// TODO make this configurable?
#define GAP_OPEN     -12
#define GAP_EXTEND    -4
#define SUB_MAT_MATCH  5
#define SUB_MAT_MIS   -4

static const IDnum SUB_MAT[4][4] = {
	{SUB_MAT_MATCH, SUB_MAT_MIS  , SUB_MAT_MIS  , SUB_MAT_MIS  },
	{SUB_MAT_MIS  , SUB_MAT_MATCH, SUB_MAT_MIS  , SUB_MAT_MIS  },
	{SUB_MAT_MIS  , SUB_MAT_MIS  , SUB_MAT_MATCH, SUB_MAT_MIS  },
	{SUB_MAT_MIS  , SUB_MAT_MIS  , SUB_MAT_MIS  , SUB_MAT_MATCH}
};


typedef struct HirschbergOp_st HirschbergOp;

struct HirschbergOp_st {
	HirschbergOp *next;
	IDnum         i1;
	IDnum         i2;
	IDnum         j1;
	IDnum         j2;
	IDnum         tb;
	IDnum         te;
	char          type;
};

struct Hirschberg_st {
	unsigned char *seq1;
	unsigned char *seq2;

	IDnum         *CC;
	IDnum         *DD;
	IDnum         *SS;
	IDnum         *RR;

	AlignCode     *aln;

	HirschbergOp  *stack;

	IDnum          l1;
	IDnum          l2;

	IDnum          score;

	IDnum 	       alnLen;

	boolean        flipped;
};

#define HIRSCHBERG_OP_BLOCK_SIZE 10000

static RecycleBin *hirschbergOpMemory = NULL;


static HirschbergOp*
newHirschbergOp(IDnum i1,
		IDnum i2,
		IDnum j1,
		IDnum j2,
		IDnum tb,
		IDnum te,
		char type)
{
	HirschbergOp *op;

	if (hirschbergOpMemory == NULL)
		hirschbergOpMemory = newRecycleBin(sizeof(HirschbergOp), HIRSCHBERG_OP_BLOCK_SIZE);
	op = allocatePointer(hirschbergOpMemory);

	op->next = NULL;
	op->i1 = i1;
	op->i2 = i2;
	op->j1 = j1;
	op->j2 = j2;
	op->tb = tb;
	op->te = te;
	op->type = type;

	return op;
}

static void
destroyHirschbergOp(HirschbergOp *op)
{
	deallocatePointer(hirschbergOpMemory, op);
}

void
cleanupHirschbergMemory(void)
{
	if (hirschbergOpMemory != NULL) {
		destroyRecycleBin(hirschbergOpMemory);
		hirschbergOpMemory = NULL;
	}
}

Hirschberg*
newHirschberg(TightString *seq1,
	      TightString *seq2)
{
	Hirschberg  *aln;
	TightString *s1;
	TightString *s2;
	IDnum l1;
	IDnum l2;
	IDnum i;

	aln = callocOrExit(1, Hirschberg);

	l1 = getLength(seq1);
	l2 = getLength(seq2);
	if (l1 > l2) {
		aln->flipped = false;
		aln->l1 = l1;
		aln->l2 = l2;
		s1 = seq1;
		s2 = seq2;
	}
	else {
		aln->flipped = true;
		aln->l1 = l2;
		aln->l2 = l1;
		s1 = seq2;
		s2 = seq1;
	}

	aln->seq1 = mallocOrExit(l1, unsigned char);
	aln->seq2 = mallocOrExit(l2, unsigned char);

	for (i = aln->l1 - 1; i >= 0; i--)
		aln->seq1[i] = getNucleotide(i, s1);
	for (i = aln->l2 - 1; i >= 0; i--)
		aln->seq2[i] = getNucleotide(i, s2);

	aln->CC = callocOrExit(aln->l2 + 1, IDnum);
	aln->DD = callocOrExit(aln->l2 + 1, IDnum);
	aln->RR = callocOrExit(aln->l2 + 1, IDnum);
	aln->SS = callocOrExit(aln->l2 + 1, IDnum);
	aln->aln = callocOrExit(aln->l1 + aln->l2 + 1, AlignCode);

	return aln;
}

void
destroyHirschberg(Hirschberg *aln)
{
	if (aln != NULL) {
		if (aln->seq1)
			free(aln->seq1);
		if (aln->seq2)
			free(aln->seq2);
		if (aln->CC)
			free(aln->CC);
		if (aln->DD)
			free(aln->DD);
		if (aln->RR)
			free(aln->RR);
		if (aln->SS)
			free(aln->SS);
		if (aln->aln)
			free(aln->aln);
		free(aln);
	}
}

static void
hirschbergFlipBack(Hirschberg *aln)
{
	unsigned char *tmpSeq;
	IDnum tmpL;
	IDnum i;

	if (!aln->flipped)
		return;

	for (i = 0; i < aln->alnLen; i++) {
		if (aln->aln[i] == ALIGN_CODE_INSERT)
			aln->aln[i] = ALIGN_CODE_DELETE;
		else if (aln->aln[i] == ALIGN_CODE_DELETE)
			aln->aln[i] = ALIGN_CODE_INSERT;
	}
	tmpSeq = aln->seq1;
	aln->seq1 = aln->seq2;
	aln->seq2 = tmpSeq;
	tmpL = aln->l1;
	aln->l1 = aln->l2;
	aln->l2 = tmpL;
	aln->flipped = false;
}

inline static IDnum
gapScore(IDnum k) {
	if (k <= 0)
		return 0;
	return GAP_OPEN + k * GAP_EXTEND;
}

static void
hirschbergAlignOneVsN(Hirschberg *aln,
		      AlignCode **ptr,
		      IDnum i1,
		      IDnum j1,
		      IDnum j2,
		      IDnum tb,
		      IDnum te)
{
	IDnum sMax;
	IDnum jMax;
	IDnum j;
	AlignCode *current = *ptr;
	const unsigned char ci = aln->seq1[i1];

	if (tb < te)
		tb = te;

	sMax = tb + GAP_EXTEND + gapScore(j2 - j1);
	jMax = -1;

	for (j = j1; j < j2; j++) {
		const unsigned char cj = aln->seq2[j];
		const IDnum s = gapScore(j - j1) + SUB_MAT[ci][cj] + gapScore(j2 -j -1);
		if (s > sMax) {
			sMax = s;
			jMax = j;
		}
	}
	if (jMax == -1) {
		for (j = j1; j < j2; j++) {
			*current++ = ALIGN_CODE_INSERT;
		}
		*current++ = ALIGN_CODE_DELETE;
	}
	else {
		for (j = j1; j < jMax; j++) {
			*current++ = ALIGN_CODE_INSERT;
		}
		if (ci == aln->seq2[jMax])
			*current++ = ALIGN_CODE_MATCH;
		else
			*current++ = ALIGN_CODE_MISMATCH;
		for (j = jMax + 1; j < j2; j++) {
			*current++ = ALIGN_CODE_INSERT;
		}
	}
	*ptr = current;
}

static void
hirschbergForward(Hirschberg *aln,
		  IDnum i1,
		  IDnum i2,
		  IDnum j1,
		  IDnum j2,
		  IDnum tb)
{
	IDnum i;
	IDnum j;
	IDnum t;

	aln->CC[j1] = 0;
	t = GAP_OPEN;
	for (j = j1 + 1; j <= j2; j++) {
		t += GAP_EXTEND;
		aln->CC[j] = t;
		aln->DD[j] = t + GAP_OPEN;
	}
	t = tb;
	for (i = i1 + 1; i <= i2; i++) {
		IDnum j;
		IDnum c;
		IDnum e;
		IDnum s;
		const unsigned char ci = aln->seq1[i - 1];

		s = aln->CC[j1];
		t += GAP_EXTEND;
		c = t;
		aln->CC[j1] = c;
		e = t + GAP_OPEN;
		for (j = j1 + 1; j <= j2; j++) {
			const unsigned char cj = aln->seq2[j - 1];

			if (e < c + GAP_OPEN)
				e = c + GAP_OPEN;
			e += GAP_EXTEND;
			if (aln->DD[j] < aln->CC[j] + GAP_OPEN)
				aln->DD[j] = aln->CC[j] + GAP_OPEN;
			aln->DD[j] += GAP_EXTEND;
			c = aln->DD[j];
			if (c < e)
				c = e;
			if (c < s + SUB_MAT[ci][cj])
				c = s + SUB_MAT[ci][cj];
			s = aln->CC[j];
			aln->CC[j] = c;
		}
	}
	aln->DD[j1] = aln->CC[j1];
}

static void
hirschbergReverse(Hirschberg *aln,
		  IDnum i1,
		  IDnum i2,
		  IDnum j1,
		  IDnum j2,
		  IDnum te)
{
	IDnum i;
	IDnum j;
	IDnum t;

	aln->RR[j2] = 0;
	t = GAP_OPEN;
	for (j = j2 - 1; j >= j1; j--) {
		t += GAP_EXTEND;
		aln->RR[j] = t;
		aln->SS[j] = t + GAP_OPEN;
	}
	t = te;
	for (i = i2 - 1; i >= i1; i--) {
		IDnum j;
		IDnum c;
		IDnum e;
		IDnum s;
		const unsigned char ci = aln->seq1[i];

		s = aln->RR[j2];
		t += GAP_EXTEND;
		c = t;
		aln->RR[j2] = c;
		e = t + GAP_OPEN;
		for (j = j2 - 1; j >= j1; j--) {
			const unsigned char cj = aln->seq2[j - 1];

			if (e < c + GAP_OPEN)
				e = c + GAP_OPEN;
			e += GAP_EXTEND;
			if (aln->SS[j] < aln->RR[j] + GAP_OPEN)
				aln->SS[j] = aln->RR[j] + GAP_OPEN;
			aln->SS[j] += GAP_EXTEND;
			c = aln->SS[j];
			if (c < e)
				c = e;
			if (c < s + SUB_MAT[ci][cj])
				c = s + SUB_MAT[ci][cj];
			s = aln->RR[j];
			aln->RR[j] = c;
		}
	}
	aln->SS[j2] = aln->RR[j2];
}

void
hirschbergAlign(Hirschberg *aln)
{
	HirschbergOp *op;
	AlignCode *current = aln->aln;
	IDnum i;
	IDnum j;

	// Start with free end gap penalties
	op = newHirschbergOp(0, aln->l1, 0, aln->l2, 0 , 0, 1);
	aln->stack = op;

	while (aln->stack != NULL) {
		IDnum i1;
		IDnum i2;
		IDnum j1;
		IDnum j2;
		IDnum tb;
		IDnum te;
		char type;

		// Pop the next op
		op = aln->stack;
		aln->stack = op->next;
		i1 = op->i1;
		i2 = op->i2;
		j1 = op->j1;
		j2 = op->j2;
		tb = op->tb;
		te = op->te;
		type = op->type;
		destroyHirschbergOp(op);

		// Handle type 2 special case
		if (type == 2) {
			*current++ = ALIGN_CODE_DELETE;
			*current++ = ALIGN_CODE_DELETE;
			continue;
		}
		// Boundary conditions
		if (i2 <= i1) {
			for (j = j1; j < j2; j++)
				*current++ = ALIGN_CODE_INSERT;
		}
		else if (j2 <= j1) {
			for (i = i1; i1 < i2; i++)
				*current++ = ALIGN_CODE_DELETE;
		}
		else if (i1 + 1 == i2) {
			hirschbergAlignOneVsN(aln, &current, i1, j1, j2, tb, te);
		}
		// Core dynamic programming algo
		else {
			const IDnum mid = (i1 + i2) / 2;
			IDnum midc;
			IDnum midj;

			hirschbergForward(aln, i1, mid, j1, j2, tb);
			hirschbergReverse(aln, mid, i2, j1, j2, te);

			type = 1;
			midc = aln->CC[j1] + aln->RR[j1];
			midj = j1;
			for (j = j1; j <= j2; j++) {
				const IDnum c = aln->CC[j] + aln->RR[j];
				if (c >= midc) {
					if (c > midc
					    || (aln->CC[j] != aln->DD[j]
						&& aln->RR[j] != aln->SS[j])) {
						midc = c;
						midj = j;
					}
				}
			}
			for (j = j2; j >= j1; j--) {
				const IDnum c = aln->DD[j] + aln->SS[j] - GAP_OPEN;
				if (c > midc) {
					midc = c;
					midj = j;
					type = 2;
				}
			}
			if (type == 1) {
				op = newHirschbergOp(mid, i2, midj, j2, GAP_OPEN, te, 1);
				op->next = aln->stack;
				aln->stack = op;
				op = newHirschbergOp(i1, mid, j1, midj, tb, GAP_OPEN, 1);
				op->next = aln->stack;
				aln->stack = op;
			}
			else {
				op = newHirschbergOp(mid + 1, i2, midj, j2, 0, te, 1);
				op->next = aln->stack;
				aln->stack = op;
				op = newHirschbergOp(mid - 1, mid, 0, 0, 0, 0, 2);
				op->next = aln->stack;
				aln->stack = op;
				op = newHirschbergOp(i1, mid - 1, j1, midj, tb, 0, 1);
				op->next = aln->stack;
				aln->stack = op;
			}
		}
	}
	aln->alnLen = current - aln->aln;
	hirschbergFlipBack(aln);
}

boolean
hirschbergCompare(Hirschberg *aln,
		  IDnum window,
		  Time divergence,
		  Time gaps)
{
	IDnum alnIdx;
	IDnum i = 0;
	IDnum j = 0;
	IDnum nMismatches = 0;
	IDnum nGaps = 0;
	IDnum maxMismatches;
	IDnum maxGaps;

	// Resize window to the smallest sequence
	if (aln->l1 < window || aln->l2 < window) {
		window = aln->l1;
		if (window < aln->l2)
			window = aln->l2;
	}

	maxMismatches = round(divergence * window);
	maxGaps = round(gaps * window);

	for (alnIdx = 0; alnIdx < aln->alnLen; alnIdx++) {
		AlignCode c = aln->aln[alnIdx];

		if (c == ALIGN_CODE_MATCH) {
			i++;
			j++;
		}
		else if (c == ALIGN_CODE_MISMATCH) {
			nMismatches++;
			i++;
			j++;
		}
		else if (c == ALIGN_CODE_INSERT) {
			nGaps++;
			j++;
		}
		else if (c == ALIGN_CODE_DELETE) {
			nGaps++;
			i++;
		}
		else {
			velvetLog("Error: unknown align code\n");
			fflush(stdout);
			abort();
		}
		if (alnIdx < window - 1)
			continue;

		if (nMismatches <= maxMismatches && nGaps <= maxGaps)
			return true;

		c = aln->aln[alnIdx - window + 1];

		if (c == ALIGN_CODE_INSERT || c == ALIGN_CODE_DELETE)
			nGaps--;
		else if (c == ALIGN_CODE_MISMATCH)
			nMismatches--;
	}
	return false;
}

void
hirschbergMapSlowOntoFast(Hirschberg *aln,
			  Coordinate *fastToSlowMapping,
			  Coordinate *slowToFastMapping)
{
	IDnum alnIdx;
	IDnum i = 0;
	IDnum j = 0;

	for (alnIdx = 0; alnIdx < aln->alnLen; alnIdx++) {
		const AlignCode c = aln->aln[alnIdx];

		if (c == ALIGN_CODE_MATCH || c == ALIGN_CODE_MISMATCH) {
			fastToSlowMapping[i] = j;
			slowToFastMapping[j] = i;
			i++;
			j++;
		}
		else if (c == ALIGN_CODE_INSERT) {
			slowToFastMapping[j] = i - 1;
			j++;
		}
		else if (c == ALIGN_CODE_DELETE) {
			fastToSlowMapping[i] = j - 1;
			i++;
		}
		else {
			velvetLog("Error: unknown align code\n");
			fflush(stdout);
			abort();
		}
	}
	fastToSlowMapping[aln->l1] = aln->l2;
	slowToFastMapping[aln->l2] = aln->l1;
}
