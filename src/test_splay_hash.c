/*
 *
 */

#include <dirent.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>

#include "globals.h"
#include "readSet.h"
#include "splay.h"
#include "kmerHashTable.h"
#include "tightString.h"
#include "crc.h"
#include "utility.h"
#include "kmer.h"

struct kmerTable_st {
	KmerHashTable *hTable;
	SplayTree **splayTrees;
	IDnum lastIndex;
	int WORDLENGTH;
	boolean useSplay;
};

KmerTable *newKmerTable(int WORDLENGTH, boolean useSplay)
{
	KmerTable *kmerTable = mallocOrExit(1, KmerTable);
	kmerTable->WORDLENGTH = WORDLENGTH;
	kmerTable->useSplay = useSplay;
	if (useSplay)
		kmerTable->splayTrees = callocOrExit(CRC_HASH_BUCKETS, SplayTree *);
	else
		kmerTable->hTable = newKmerHashTable();
	kmerTable->lastIndex = 0;
	return kmerTable;
}

void destroyKmerTable(KmerTable * kmerTable)
{
	velvetLog("Destroying kmer table\n");

	if (kmerTable->useSplay)
	{
		destroyAllSplayTrees();
		free(kmerTable->splayTrees);
	}
	else
		destroyKmerHashTable(kmerTable->hTable);
	free(kmerTable);

	velvetLog("Kmer table destroyed\n");
}

static int hash_kmer(Kmer * kmer)
{
	//return crc32_v((char *) kmer, KMER_BYTE_SIZE);
#define CRC_HASH_MASK    0x0000000000ffffffL
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

	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);

	return key & CRC_HASH_MASK;
#elif KMER_LONGS
	KmerKey key = kmer->longs[0];

	key += ~(key << 15);
	key ^=  (key >> 10);
	key +=  (key << 3);
	key ^=  (key >> 6);
	key += ~(key << 11);
	key ^=  (key >> 16);

	return key & CRC_HASH_MASK;

#elif KMER_INTS
	return kmer->ints & CRC_HASH_MASK;
#elif KMER_CHARS
	return kmer->chars & CRC_HASH_MASK;
#endif
}

static boolean findOrInsertOccurenceInKmerTable(Kmer * kmer, IDnum * seqID,
						 Coordinate * position,
						 KmerTable * kmerTable)
{
	if (kmerTable->useSplay)
		return findOrInsertOccurenceInSplayTree(kmer, seqID, position, &kmerTable->splayTrees[hash_kmer(kmer)]);
	return findOrInsertOccurenceInKmerHashTable(kmerTable->hTable, kmer, seqID, position);
}

static void inputSequenceIntoKmerTable(TightString * tString, KmerTable * kmerTable) 
{
	IDnum currentIndex;
	Coordinate readNucleotideIndex = 0;
	Coordinate writeNucleotideIndex = 0;
	Kmer word;
	Kmer antiWord;
	IDnum sequenceID;
	Coordinate coord;
	Nucleotide nucleotide;

	clearKmer(&word);
	clearKmer(&antiWord);

	kmerTable->lastIndex++;

	currentIndex = kmerTable->lastIndex;

	// Neglect any string shorter than WORDLENGTH :
	if (GETlENGTH(TsTRING) < KMERtABLE->wordlength) {
		return;
	}
	// Fill in the initial word :
	for (readNucleotideIndex = 0;
	     readNucleotideIndex < kmerTable->WORDLENGTH - 1;
	     readNucleotideIndex++) { 
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);
		reversePushNucleotide(&antiWord, 3 - nucleotide);
	}

	while (readNucleotideIndex < getLength(tString)) {
		// Shift word:
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);
		reversePushNucleotide(&antiWord, 3 - nucleotide);

		sequenceID = currentIndex;
		coord = writeNucleotideIndex;

		if (compareKmers(&word, &antiWord) <= 0) {
			findOrInsertOccurenceInKmerTable(&word,
					&sequenceID,
					&coord,
					kmerTable);
		} else {
			sequenceID = -sequenceID;
			findOrInsertOccurenceInKmerTable(&antiWord,
					&sequenceID,
					&coord,
					kmerTable);
			sequenceID = -sequenceID;
		}
		writeNucleotideIndex++;
		readNucleotideIndex++;
	}
	return;
}

static
void inputSequenceArrayIntoKmerTableAndArchive(ReadSet * reads,
						KmerTable * kmerTable)
{
	IDnum index;
	IDnum sequenceCount = reads->readCount;
	TightString *array;

	velvetLog("Inputting sequences...\n");
	array = reads->tSequences;
	for (index = 0; index < sequenceCount; index++)
	{
		/*
		if (index % 1000000 == 0)
			velvetLog("Inputting sequence %d / %d\n", index, sequenceCount);
		*/
		inputSequenceIntoKmerTable(getTightStringInArray(array, index), kmerTable);
	}
	velvetLog("Done inputting sequences\n");
}

static void printUsage()
{
	puts("Usage:");
	puts("./test_splay_hash directory hash_length {[-file_format][-read_type] filename}");
	puts("");
	puts("\tdirectory\t\t: directory name for output files");
	printf("\thash_length\t\t: odd integer (if even, it will be decremented) <= %i (if above, will be reduced)\n", MAXKMERLENGTH);
	puts("\tfilename\t\t: path to sequence file or - for standard input");	
	puts("");
	puts("File format options:");
	puts("\t-fasta");
	puts("\t-fastq");
	puts("\t-raw");
	puts("\t-fasta.gz");
	puts("\t-fastq.gz");
	puts("\t-raw.gz");
	puts("\t-sam");
	puts("\t-bam");
	puts("\t-eland");
	puts("\t-gerald");
	puts("");
	puts("Read type options:");
	puts("\t-short");
	puts("\t-shortPaired");
	puts("\t-short2");
	puts("\t-shortPaired2");
	puts("");
	puts("Output:");
	puts("\tdirectory/Sequences");
	puts("\t\t[Both files are picked up by graph, so please leave them there]");
}

int main(int argc, char **argv)
{
	ReadSet *allSequences;
	KmerTable *kmerTable;
	int hashLength;
	char *directory, *filename, *seqFilename, *buf;
	DIR *dir;
	FILE *seqFile;
	boolean double_strand;
	struct timeval start, end, diff;
	int i;

	setProgramName("test_splay_hash");

	if (argc < 4) {
		printUsage();
		return 0;
	}

	directory = argv[1];
	filename = mallocOrExit(strlen(directory) + 100, char);
	seqFilename = mallocOrExit(strlen(directory) + 100, char);
	buf = mallocOrExit(strlen(directory) + 100, char);

	hashLength = atoi(argv[2]);

	if (hashLength > MAXKMERLENGTH) {
		printf
		    ("Velvet can't handle k-mers as long as %i! We'll stick to %i if you don't mind.\n",
		     hashLength, MAXKMERLENGTH);
		hashLength = MAXKMERLENGTH;
	} else if (hashLength <= 0) {
		printf("Invalid hash length: %s\n", argv[2]);
		printUsage();
		return 0;
	} 

	if (hashLength % 2 == 0) {
		printf
		    ("Velvet can't work with even length k-mers, such as %i. We'll use %i instead, if you don't mind.\n",
		     hashLength, hashLength - 1);
		hashLength--;
	}
	resetWordFilter(hashLength);

	dir = opendir(directory);

	if (dir == NULL)
		mkdir(directory, 0777);

	strcpy(seqFilename, directory);
	strcat(seqFilename, "/Sequences");
	seqFile = fopen(seqFilename, "r");
	if (seqFile == NULL)
		parseDataAndReadFiles(seqFilename, argc - 2, &(argv[2]), &double_strand);
	else
		fclose(seqFile);

	allSequences = importReadSet(seqFilename);
	velvetLog("%i sequences in total.\n", allSequences->readCount);
	convertSequences(allSequences);

	for (i = 0; i < 10; i++)
	{
		gettimeofday(&start, NULL);
		kmerTable = newKmerTable(hashLength, true);
		inputSequenceArrayIntoKmerTableAndArchive(allSequences, kmerTable);
		destroyKmerTable(kmerTable);
		gettimeofday(&end, NULL);
		timersub(&end, &start, &diff);
		velvetLog("Iterated with splay trees in %ld.%06ld s\n", diff.tv_sec, diff.tv_usec);
	}
	for (i = 0; i < 10; i++)
	{
		gettimeofday(&start, NULL);
		kmerTable = newKmerTable(hashLength, false);
		inputSequenceArrayIntoKmerTableAndArchive(allSequences, kmerTable);
		destroyKmerTable(kmerTable);
		gettimeofday(&end, NULL);
		timersub(&end, &start, &diff);
		velvetLog("Iterated with hash tables in %ld.%06ld s\n", diff.tv_sec, diff.tv_usec);
	}
	free(allSequences->tSequences);
	allSequences->tSequences = NULL;
	destroyReadSet(allSequences);

	if (dir)
		closedir(dir);
	free(filename);
	free(seqFilename);
	free(buf);

	return 0;
}
