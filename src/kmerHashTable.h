/**
 *
 */

#ifndef _HASH_TABLE_DAQP_H_
#define _HASH_TABLE_DAQP_H_

#include "globals.h"


typedef struct kmerHashTable_st KmerHashTable;

KmerHashTable* newKmerHashTable(void);
void destroyKmerHashTable(KmerHashTable *hTable);
boolean findOrInsertOccurenceInKmerHashTable(KmerHashTable *hTable, Kmer *kmer,
					     IDnum *seqID, Coordinate *position);

#endif /* _HASH_TABLE_DAQP_H_ */
