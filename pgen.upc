#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include <upc_io.h>
#include <stdbool.h>

#include "packingDNAseq.h"
#include "kmer_hash.h"

typedef struct{
	char kmer[KMER_PACKED_LENGTH];
	char l_ext;	
	char r_ext;
	int64_t next; // stores the index of the next kmer in hashtable
} kmer_upc_t;

typedef struct{	
	shared[1] int64_t* list; // stores the starting kmers' indices
	int64_t size;
} start_list_upc_t;

typedef struct{
	shared[1] kmer_upc_t* heap; // stores the kmers
} memory_heap_upc_t;

typedef struct{
	shared[1] int64_t* tableHead; // stores the index of the head in each bucket
} hash_table_upc_t;

void create_all_upc(int64_t nKmers, int64_t nKmersPerThread, int64_t nBuckets, 
	shared memory_heap_upc_t* memoryHeap, shared hash_table_upc_t* hashTable,
	shared start_list_upc_t* startList)
{
	//Init memory heap
	memoryHeap->heap = (shared[1] kmer_upc_t*) upc_all_alloc(THREADS * nKmersPerThread, sizeof(kmer_upc_t));
	if(memoryHeap->heap == NULL) {
		fprintf(stderr, "ERROR: Could not allocate memory for the heap!\n");
      	exit(1);
	}
	
	//Init hashtable
	hashTable->tableHead = (shared[1] int64_t*) upc_all_alloc(nBuckets, sizeof(int64_t));
	if(hashTable->tableHead == NULL) {
		fprintf(stderr, "ERROR: Could not allocate memory for the hash table: %lld buckets of %lu bytes\n", nBuckets, sizeof(int64_t));
      	exit(1);
	}

	//********************************** weird!!!!!!!!
	shared int64_t* bucketList = hashTable->tableHead; // Why don't we directly use hashTable??????

	upc_forall(int64_t i = 0; i < nBuckets; i++; &bucketList[i])
	{
		bucketList[i] = -1;
	}

	//Init start list
	startList->list = (shared[1] int64_t*) upc_all_alloc(nKmers, sizeof(int64_t));
	if (startList->list == NULL) {
		fprintf(stderr, "ERROR: Could not allocate memory for the start list!\n");
		exit(1);
	}
	if (MYTHREAD == 0)
		startList->size = 0;

}

void add_kmer_to_start_list_upc(start_list_upc_t* startList, int64_t kmerIdx)
{
   int64_t* size_ptr = &(startList->size); //??????????????????????????????
   int64_t index = bupc_atomicU64_fetchadd_relaxed(size_ptr, 1);
   startList->list[index] = kmerIdx;
}

void add_kmer_upc(hash_table_upc_t* hashTable, memory_heap_upc_t* memoryHeap, const unsigned char *kmer, 
			char left_ext, char right_ext, int64_t kmerIdx, int64_t nBuckets)
{
	/* Pack a k-mer sequence appropriately */
	char packedKmer[KMER_PACKED_LENGTH];
	packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
	int64_t hashVal = hashkmer(nBuckets, (char*) packedKmer);

	/* Add the contents to the appropriate kmer struct in the heap */
	kmer_upc_t new_kmer;
	new_kmer.l_ext = left_ext;
	new_kmer.r_ext = right_ext;
	new_kmer.next = -1;
	memcpy(new_kmer.kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));



	//memcpy((memoryHeap+kmerIdx)->kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));
	//(memoryHeap+kmerIdx)->l_ext = left_ext;
	//(memoryHeap+kmerIdx)->r_ext = right_ext;

	/*Add the contents to the appropriate bucket in the hashtable */
	shared[1] kmer_upc_t* kmer_ptr = memoryHeap->heap + kmerIdx;
	*kmer_ptr = new_kmer;

	shared[1] int64_t* bucket_ptr = hashTable->tableHead + hashVal;

	int64_t old = -1;
	int64_t bucket_index = old;
	do
	{
		kmer_ptr->next = bucket_index;
		old = bucket_index;
		bucket_index = bupc_atomicI64_cswap_relaxed(bucket_ptr, old, kmerIdx);
	} while(bucket_index != old);
	///////////////check you pointer!!!!!!!!
	//shared[1] int64_t curBucketVal = *(hashTable + hashVal); //Assume the current bucket value is -1 (no kmers hashed into the bucket yet)
	//int64_t oldVal;

	//do{
	//	(memoryHeap+kmerIdx)->next = curBucketVal;
		//oldVal = curBucketVal;
	//	newBucketVal = bupc_atomicI64_cswap_relaxed(hashTable+hashVal, curBucketVal, kmerIdx);
	//} while (newBucketVal != curBucketVal)

}

kmer_upc_t lookup_kmer_upc(memory_heap_upc_t* memoryHeap, hash_table_upc_t* hashTable, 
							const unsigned char *kmer, int64_t nBuckets)
{
   char packedKmer[KMER_PACKED_LENGTH];
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   int64_t hashVal = hashkmer(nBuckets, (char*) packedKmer);
   int64_t kmerIdx = hashTable->tableHead[hashVal];
   kmer_upc_t result = memoryHeap->heap[kmerIdx];
   
   while (kmerIdx != -1) {
      if ( memcmp(packedKmer, result.kmer, KMER_PACKED_LENGTH * sizeof(char)) == 0 ) {
         break;
      }
      kmerIdx = result.next;
      result = memoryHeap->heap[kmerIdx];
   }
   return result;
}


int dealloc_memoryheap_upc(shared memory_heap_upc_t* memoryHeap)
{
    if(MYTHREAD == 0)
    	upc_free(memoryHeap->heap);
    return 0;
}

int dealloc_hashtable_upc(shared hash_table_upc_t* hashTable)
{
   	if(MYTHREAD == 0)
   		upc_free(hashTable->tableHead);
   	return 0;
}

int dealloc_startlist_upc(shared start_list_upc_t* startList)
{
	if(MYTHREAD == 0)
		upc_free(startList->list);
	return 0;
}

int main(int argc, char *argv[]){

	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;


	/** Read input **/
	upc_barrier;
	inputTime -= gettime();
	///////////////////////////////////////////
	// Your code for input file reading here //
	///////////////////////////////////////////
	char* input_UFX_name = argv[1];
	static shared int64_t nKmers; // shared by all threads
	if(MYTHREAD == 0) nKmers = getNumKmersInUFX(input_UFX_name); // once is enough

	upc_barrier;

	int64_t nKmersPerThread = (nKmers + THREADS - 1) / THREADS;
	int64_t nKmersThisThread = min(nKmersPerThread * (MYTHREAD + 1), nKmers) - nKmersPerThread * MYTHREAD;
	int64_t nBuckets = nKmers * LOAD_FACTOR;
	nBuckets = (nBuckets + THREADS - 1) / THREADS * THREADS;

	unsigned char* working_buffer;
	int64_t total_chars_to_read = nKmersThisThread * LINE_SIZE;
   	working_buffer = (unsigned char*) malloc(total_chars_to_read * sizeof(unsigned char));
   	FILE* inputFile = fopen(input_UFX_name, "r");
   	fseek(inputFile, nKmersPerThread * MYTHREAD * LINE_SIZE * sizeof(unsigned char), 0);
   	int64_t cur_chars_read = fread(working_buffer, sizeof(unsigned char), total_chars_to_read , inputFile);

   	fclose(inputFile);


   	upc_barrier;
	inputTime += gettime();

	/** Graph construction **/
	constrTime -= gettime();
	///////////////////////////////////////////
	// Your code for graph construction here //
	///////////////////////////////////////////
	init_LookupTable();
	static shared memory_heap_upc_t memoryHeap;
	static shared hash_table_upc_t hashTable;
	static shared start_list_upc_t startList;

	create_all_upc(nKmers, nKmersPerThread, nBuckets, &memoryHeap, &hashTable, &startList);

	int64_t startIdx = nKmersPerThread * MYTHREAD;
	//put_kmers_in_hashtable(hashtable, memory_heap, &start_list, aux, buffer, buffer_size, start_index);
	int64_t ptr = 0;
	int64_t kmerIdx = startIdx;
	char left_ext, right_ext;
	while(ptr < cur_chars_read) {
		left_ext = (char)working_buffer[ptr+KMER_LENGTH+1];
		right_ext = (char)working_buffer[ptr+KMER_LENGTH+2];
		add_kmer_upc(&hashTable, &memoryHeap, &working_buffer[ptr], left_ext, right_ext, kmerIdx, nBuckets);
		if(left_ext == "F")
			add_kmer_to_start_list_upc(&startList, kmerIdx);
		ptr += LINE_SIZE;
		kmerIdx++;
	}
	free(working_buffer);     

	upc_barrier;
	constrTime += gettime();

	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////
	char* filename = (char*) malloc((strlen("output") + 5 + 64) * sizeof(char));
	sprintf(filename, "%s_%d.out", "output", MYTHREAD);
	FILE* upcOutputFile = fopen(filename, "w");
	free(filename);
	///////////////
	/* Pick start nodes from the startKmersList */
    char unpackedKmer[KMER_LENGTH+1];
    unpackedKmer[KMER_LENGTH] = "\0";
    char cur_contig[MAXIMUM_CONTIG_SIZE];


    shared int64_t* list_ptr = startList.list;
    upc_forall(int64_t i = 0; i < startList.size; i++; &list_ptr[i]) //??????????????????????
    {
    	int64_t startKmerIdx = startList.list[i];
    	kmer_upc_t curKmer = memoryHeap.heap[startKmerIdx];
    	unpackSequence((unsigned char*) curKmer.kmer, unpackedKmer, KMER_LENGTH);
    	memcpy(cur_contig, unpackedKmer, KMER_LENGTH * sizeof(char));
    	int64_t posInContig = KMER_LENGTH;
    	char right_ext = curKmer.r_ext;

    	while(right_ext != 'F')
    	{
    		cur_contig[posInContig] = right_ext;
    		posInContig++;
    		curKmer = lookup_kmer_upc(&memoryHeap, &hashTable, (const unsigned char *)&cur_contig[posInContig-KMER_LENGTH], nBuckets);
    		right_ext = curKmer.r_ext;
    	}

    	cur_contig[posInContig] = '\0';
    	fprintf(upcOutputFile, "%s\n", cur_contig);
    }

	fclose(upcOutputFile);
	////////////////////////

	upc_barrier;
	traversalTime += gettime();

	dealloc_startlist_upc(&startList);
	dealloc_memoryheap_upc(&memoryHeap);
	dealloc_hashtable_upc(&hashTable);

	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}
