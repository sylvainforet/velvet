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

#ifndef _HIRSCHBERG_H_
#define _HIRSCHBERG_H_

#include "globals.h"

typedef enum {
	ALIGN_CODE_NONE = 0,
	ALIGN_CODE_MATCH,
	ALIGN_CODE_MISMATCH,
	ALIGN_CODE_DELETE,
	ALIGN_CODE_INSERT
}
AlignCode;

typedef struct Hirschberg_st Hirschberg;


Hirschberg* newHirschberg(TightString *seq1,
			  TightString *seq2);

void destroyHirschberg(Hirschberg *aln);

void cleanupHirschbergMemory(void);

void hirschbergAlign(Hirschberg *aln);

boolean hirschbergCompare(Hirschberg *aln,
			  IDnum window,
			  Time divergence,
			  Time gaps);

void hirschbergMapSlowOntoFast(Hirschberg *aln,
			       Coordinate *fastToSlowMapping,
			       Coordinate *slowToFastMapping);

void hirschbergPrint(Hirschberg *aln);

#endif /* _HIRSCHBERG_H_ */
