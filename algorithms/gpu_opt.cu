#include <iostream>
#include <cmath>
#include <string.h>
#include <thrust/scan.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>

# if !defined USED_READS_SIZE
# define USED_READS_SIZE 2048
# endif

# if !defined ETH
# define ETH 10
# endif

# if !defined KMER
# define KMER 12
# endif

# if !defined DIVIDE_DATA_BY
# define DIVIDE_DATA_BY 1
# endif

# if !defined UCHAR4_OPTIMIZATION
# define UCHAR4_OPTIMIZATION true
# endif

# if !defined FORCE_UCHAR4_OPTIMIZATION
# define FORCE_UCHAR4_OPTIMIZATION false
# endif

#define NEW_VAL(s_val, t_val) min(diag + (s[i - 1] != t[j - 1]), min(s_val + (s[i - 1] != 'P'),  t_val + (t[j - 1] != 'P')));

const int THREADS_PER_BLOCK = 32;

struct Histogram{
    unsigned char numA;
    unsigned char numT;
    unsigned char numC;
    unsigned char numG;
};

struct IndexTable{
    int index;
    int count;
};

template <unsigned int n, typename T>
/**
 * A cyclic buffer of size n
 */
struct CyclicBuffer{
    /**
     * The buffer itself
     */
    T buffer[n];
    /**
     * The buffer start location
     */
    unsigned int start;
    /**
     * The buffer length
     */
    unsigned int length;
    /**
     * Resets the buffer
     */
    __device__ void reset(){
        start = 0; length = 0;
    }
    /**
     * Returns the xth element in the buffer
     * @param x
     * @return
     */
    __device__ T get(unsigned int x){
        return buffer[(start + x) % n];
    }
    /**
     * Sets the xth element in the buffer
     * @param x
     * @param data
     */
    __device__ void set(unsigned int x, T data){
        buffer[(start + x) % n] = data;
    }
    /**
     * Enlarges the buffer by numElements. Should only be called by thread 0.
     * @param numElements
     */
    __device__ void push(unsigned int numElements){
        if (n <= length + numElements){
            printf("length is bigger then n\n");
            return;
        }
        if(threadIdx.x == 0) length += numElements;
    }
    /**
     * Shrinks the buffer by numElements. Should only be called by thread 0.
     * @param numElements
     */
    __device__ void pop(unsigned int numElements){
        if(threadIdx.x == 0){
            start = (start + numElements) % n;
            if (length < numElements){
                printf("length is smaller then numElements\n");
//                exit(1);
            }
            length -= numElements;
        }
    }
    /**
     * Returns true if buffer contains data.
     * @param data
     */
    __device__ bool contains(T data){
        for(int i = 0; i < length; i++){
            if (get(i) == data){
                return true;
            }
        }
        return false;
    }
};

__device__ int editDistance(const char* s, const char* t){
#if (READ_LENGTH <= 115 || !UCHAR4_OPTIMIZATION) && !FORCE_UCHAR4_OPTIMIZATION
    // The last row
    int arr[READ_LENGTH + 1];

    // Initialize arr1 to be the first row of the DP matrix
    for(int j = 0; j <= READ_LENGTH; j++){
        arr[j] = j;
    }

    // Fill the remaining rows
    for(int i = 1; i <= READ_LENGTH; i++) {

        int diag = arr[0];
        arr[0] = i;

        for (int j = 1; j <= READ_LENGTH; j++) {

            int new_val = min(diag + (s[i - 1] != t[j - 1]),
                              min(arr[j] + (s[i - 1] != 'P'), arr[j - 1] + (t[j - 1] != 'P')));
            diag = arr[j];
            arr[j] = new_val;

        }

    }
    return arr[READ_LENGTH];
#else
    uchar4 arr[READ_LENGTH/4 + 1];

    // Initialize arr to be the first row of the DP matrix
    for (int j = 0; j <= READ_LENGTH; j+=4){
        arr[j/4].x = j;
        arr[j/4].y = j+1;
        arr[j/4].z = j+2;
        arr[j/4].w = j+3;
    }

    // Fill the remaining rows
    for(int i = 1; i <= READ_LENGTH; i++) {

        int diag = arr[0].x;
        arr[0].x = i;

        for (int j = 1; j <= READ_LENGTH; j++) {
            switch (j%4){
                case 0:
                {
                    // x
                    int new_val = NEW_VAL(arr[j/4].x, arr[j/4 - 1].w);
                    diag = arr[j/4].x;
                    arr[j/4].x = new_val;
                }
                    break;
                case 1:
                {
                    // y
                    int new_val = NEW_VAL(arr[j/4].y, arr[j/4].x);
                    diag = arr[j/4].y;
                    arr[j/4].y = new_val;
                }
                    break;
                case 2:
                {
                    // z
                    int new_val = NEW_VAL(arr[j/4].z, arr[j/4].y);
                    diag = arr[j/4].z;
                    arr[j/4].z = new_val;
                }
                    break;
                case 3:
                {
                    // w
                    int new_val = NEW_VAL(arr[j/4].w, arr[j/4].z);
                    diag = arr[j/4].w;
                    arr[j/4].w = new_val;
                }
                    break;
            }
        }
    }

    int last =  READ_LENGTH/4;
    switch (READ_LENGTH % 4){
        case 0: return arr[last].x;
        case 1: return arr[last].y;
        case 2: return arr[last].z;
        case 3: return arr[last].w;
    }

#endif
}


__constant__ int RHO[26] = {0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0};//possibilities

/**
 * Calculates the decimal value of a DNA sequence.
 * @param dnaSeq
 */
__device__ int dnaToDecimal(const char* dnaSeq) {
    int decimalNum = 0;
    for (int i = KMER-1; i >= 0; i--) {
        char base = dnaSeq[i];
        if (base < 'A' || base > 'Z'){
            printf("illegal base: %c", base);
            return 0;
        }
        int nucleotideValue = RHO[base - 'A'];
        decimalNum = decimalNum * 4 + nucleotideValue;
    }
    return decimalNum;
}

/**
 * Calculates the base counts for each read.
 * @param reads
 * @param histograms
 * @param num_reads
 */
__global__ void computeHistogram(const char *reads, Histogram *histograms, int num_reads){
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= num_reads) return;
    histograms[index].numA = 0;
    histograms[index].numT = 0;
    histograms[index].numG = 0;
    histograms[index].numC = 0;
    for(int i = 0; i < READ_LENGTH; i++){
        histograms[index].numA += (reads[index * READ_LENGTH + i] == 'A');
        histograms[index].numT += (reads[index * READ_LENGTH + i] == 'T');
        histograms[index].numG += (reads[index * READ_LENGTH + i] == 'G');
        histograms[index].numC += (reads[index * READ_LENGTH + i] == 'C');
    }
}


/**
 * Calculates the number of reads containing each KMER
 * @param reads
 * @param read_counts
 * @param num_reads
 */
__global__ void computeReadCounts(const char *reads, int *read_counts, int num_reads){
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= num_reads) return;
    for (int i = 0; i <= READ_LENGTH - KMER; i++){
        char kmer[KMER];
        for (int j = 0; j < KMER; j++) kmer[j] = reads[(READ_LENGTH * index)+j+i];
        atomicAdd(&read_counts[dnaToDecimal(kmer)], 1);
    }
}


/**
 * initializes an array of IndexTable instances that for each KMER contains the first read in the read chunk and the number of reads containing the KMER
 * @param reads
 * @param read_counts
 * @param tmp_read_counts
 * @param read_chunk
 * @param prefix
 * @param index_table
 * @param num_reads
 */
__global__ void computeIndexTable(const char *reads, int *read_counts, int *tmp_read_counts, int *read_chunk, int *prefix, IndexTable *index_table, int num_reads){
    int read = blockIdx.x * blockDim.x + threadIdx.x;
    if (read >= num_reads) return;
    for (int i = 0; i <= READ_LENGTH - KMER; i++){
        char kmer[KMER];
        for (int j = 0; j < KMER; j++) kmer[j] = reads[(READ_LENGTH * read)+j+i];
        int kmerNum = dnaToDecimal(kmer);
        int value = atomicAdd(&tmp_read_counts[kmerNum], 1);
        int index = prefix[kmerNum] + value;
        read_chunk[index] = read;
        index_table[kmerNum].index = prefix[kmerNum];
        index_table[kmerNum].count = read_counts[kmerNum];
    }
}

__global__ void findClosest(const char *reads, int *min_num, int *min_index, Histogram *histograms, IndexTable *index_table, int *read_chunk, unsigned long long *clock, long *counter_sum){
    unsigned long long start = clock64();
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    long counter = 0;

    // will contain all the potential closest reads
    __shared__ CyclicBuffer <2 * THREADS_PER_BLOCK, int> samplesBuffer;
    samplesBuffer.reset();

    // will contain the first USED_READS_SIZE(number) that were checked with edit distance
    __shared__ CyclicBuffer <USED_READS_SIZE, int> usedReads;
    usedReads.reset();

    __shared__ int min_distance;
    __shared__ int minIdx;
    __shared__ int count;
    __shared__ char read[READ_LENGTH];


    //creates the main read in the memory:
    if (threadIdx.x == 0) {
        min_distance = READ_LENGTH;
        minIdx = -1;
        for (int i = 0; i < READ_LENGTH; i++) read[i] = reads[(READ_LENGTH * blockIdx.x) + i];
    }
    __syncthreads();

    // iterate over each KMER in the block's read
    for (int i = 0; i <= READ_LENGTH - KMER; i++) {
        // copy the KMER to array
        char kmer[KMER];
        for (int j = 0; j < KMER; j++) kmer[j] = read[j + i];

        IndexTable kmer_reads = index_table[dnaToDecimal(kmer)];

        // iterate over each read containing the current KMER
        for (int j = 0; j < kmer_reads.count; j += THREADS_PER_BLOCK) {
            if (j + threadIdx.x < kmer_reads.count) {

                int comp_read = read_chunk[kmer_reads.index + j + threadIdx.x];
                count = 0;
                if (blockIdx.x != comp_read && !usedReads.contains(comp_read)) {
                    int samplesLength = samplesBuffer.length;

                    // calculates base counts diff between block read and another read.
                    int diff = abs(histograms[blockIdx.x].numA - histograms[comp_read].numA) +
                               abs(histograms[blockIdx.x].numT - histograms[comp_read].numT) +
                               abs(histograms[blockIdx.x].numG - histograms[comp_read].numG) +
                               abs(histograms[blockIdx.x].numC - histograms[comp_read].numC);

                    // if diff is smaller the 2*ETH (this is the filter) add it to the buffer
                    if (diff < 2 * ETH) {
                        int readBufferIdx = atomicAdd(&count, 1);
                        samplesBuffer.set(readBufferIdx + samplesLength, comp_read);
                    }
                }
            }
            __syncthreads();

            if (threadIdx.x == 0) {
                samplesBuffer.push(count); //enlarges the length of the buffer
            }
            __syncthreads();

            // if buffer is big enough, check 32 reads from the buffer
            if (samplesBuffer.length >= THREADS_PER_BLOCK) {
                counter++;
                // calculate edit distance
                int edit_distance = editDistance(read, (reads + READ_LENGTH * (samplesBuffer.get(threadIdx.x))));
                if (edit_distance < min_distance) {
                    int previous = atomicMin(&min_distance, edit_distance);
                    // make sure that this thread has the actual smallest edit distance
                    if (edit_distance < previous && edit_distance == min_distance) {
                        minIdx = samplesBuffer.get(threadIdx.x);
                    }
                }

                // removes checked reads and adds them the used reads buffer
                if (threadIdx.x == 0){
                    if (usedReads.length + THREADS_PER_BLOCK < USED_READS_SIZE){
                        usedReads.push(THREADS_PER_BLOCK);
                    }else{
                        usedReads.push(USED_READS_SIZE - usedReads.length -1);
                    }
                }
                int threadLoc = usedReads.length + threadIdx.x;
                __syncthreads();
                if (threadLoc < USED_READS_SIZE){
                    usedReads.set(threadLoc, samplesBuffer.get(threadIdx.x));
                }
                if (threadIdx.x == 0) {
                    samplesBuffer.pop(THREADS_PER_BLOCK);
                }
                __syncthreads();
            }

        }
    }
    // empty samplesBuffer one last time
    if(samplesBuffer.length > 0) {
        if (threadIdx.x < samplesBuffer.length){
            int edit_distance = editDistance(read, (reads + READ_LENGTH * (samplesBuffer.get(threadIdx.x))));
            counter++;

            if (edit_distance < min_distance) {
                int previous = atomicMin(&min_distance, edit_distance);
                if (edit_distance < previous && edit_distance == min_distance) {
                    minIdx = samplesBuffer.get(threadIdx.x);
                }
            }
        }
        __syncthreads();

    }

    if(threadIdx.x == 0){
        min_num[blockIdx.x] = min_distance;
        min_index[blockIdx.x] = minIdx;

    }
    __syncthreads();
    unsigned long long end = clock64();
    atomicAdd(clock, (end - start));
    counter_sum[index] = counter;
}

int main(){

    // receive dataset as input
    std::string readsStr;
    std::cin >> readsStr;
    const char *reads = readsStr.c_str();
    char *d_reads; cudaMalloc(&d_reads,reads_length * sizeof(char));
    cudaMemcpy(d_reads, reads, reads_length * sizeof(char), cudaMemcpyHostToDevice);

    // get file length
    int reads_length = readsStr.length();

    // calculate the number of reads in the dataset
    int num_reads = reads_length / READ_LENGTH;

    // create histogram in the GPU memory
    Histogram *d_histograms; cudaMalloc(&d_histograms, num_reads * sizeof(Histogram));

    // initialize clock in CPU and GPU memory
    unsigned long long clock = 0;
    unsigned long long *d_clock; cudaMalloc(&d_clock, sizeof(unsigned long long));

    // initialize a counter that counts every edit distance calculated (used for testing purposes only)
    long *counts = (long*) malloc(num_reads * THREADS_PER_BLOCK * sizeof(long));
    long *d_counts; cudaMalloc(&d_counts, num_reads * THREADS_PER_BLOCK * sizeof(long));

    // initialize an array of the minimum edit distance for each read
    int *min_num = (int*) malloc(num_reads * sizeof(int));
    int *d_min_num; cudaMalloc(&d_min_num, num_reads * sizeof(int));

    // initialize an array of the index of the minimum edit distance for each read
    int *min_index = (int*) malloc(num_reads * sizeof(int));
    int *d_min_index; cudaMalloc(&d_min_index, num_reads * sizeof(int));

    // initialize an array of the counts of each KMER
    int *d_read_counts; cudaMalloc(&d_read_counts, std::pow(4,KMER) * sizeof(int)); cudaMemset(d_read_counts, 0, std::pow(4,KMER) * sizeof(int));

    // initialize an array that sum of all the values before each value in the counts array.
    // will be used to locate each KMER in the chunk array
    int *read_prefix = (int*) malloc(std::pow(4,KMER) * sizeof(int));
    int *d_read_prefix; cudaMalloc(&d_read_prefix, std::pow(4,KMER) * sizeof(int));

    // creating a device(GPU) pointer enables to use reduce in the thrust library
    thrust::device_ptr<int> dev_ptr(d_read_counts);

    // calculate the base counts for each read
    computeHistogram<<<num_reads/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(d_reads,d_histograms, num_reads);

    // calculates the number of reads containing each KMER
    computeReadCounts<<<num_reads/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(d_reads, d_read_counts, num_reads);

    // sums all the values in read counts
    int read_counts_sum = thrust::reduce(thrust::device, dev_ptr, dev_ptr + std::pow(4,KMER),0);

    // calculate the prefix
    thrust::exclusive_scan(thrust::device, dev_ptr, dev_ptr + std::pow(4,KMER), d_read_prefix);

    // initialize an array that will contain the actual reads sorted by KMER
    int *d_read_chunk; cudaMalloc(&d_read_chunk, read_counts_sum * sizeof(int));
    int *read_chunk = (int*) malloc(read_counts_sum * sizeof(int));

    // an array that will contain the IndexTable instances that for each KMER contains an index to the first read in read chunk and the amount of reads
    IndexTable *d_index_table; cudaMalloc(&d_index_table, std::pow(4,KMER) * sizeof(IndexTable));cudaMemset(d_index_table, 0, std::pow(4,KMER) * sizeof(IndexTable));
    IndexTable *index_table = (IndexTable*) malloc(std::pow(4,KMER) * sizeof(IndexTable));
    int *d_tmp_read_counts; cudaMalloc(&d_tmp_read_counts, std::pow(4,KMER) * sizeof(int)); cudaMemset(d_tmp_read_counts, 0, std::pow(4,KMER) * sizeof(int));

    // calculate index table
    computeIndexTable<<<num_reads/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(d_reads, d_read_counts, d_tmp_read_counts, d_read_chunk, d_read_prefix, d_index_table, num_reads);

    // delete no longer used arrays
    cudaFree(d_read_counts);
    cudaFree(d_tmp_read_counts);

    // for each read finds the closest read and its distance
    findClosest<<<num_reads/DIVIDE_DATA_BY,THREADS_PER_BLOCK>>>(d_reads,d_min_num, d_min_index, d_histograms,d_index_table,d_read_chunk, d_clock,d_counts);

    // clock that calculates average thread runtime, it was used for testing and debugging
    cudaMemcpy(&clock, d_clock, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
    cudaMemcpy(counts, d_counts, num_reads * THREADS_PER_BLOCK * sizeof(long), cudaMemcpyDeviceToHost);
    thrust::device_ptr<long> dev_counts_ptr(d_counts);

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    int clockRate = prop.clockRate * 1000;// its in KHz so * 1000 to Hz

    // copy results to CPU memory
    cudaMemcpy(min_num, d_min_num, num_reads * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(min_index, d_min_index, num_reads * sizeof(int), cudaMemcpyDeviceToHost);

    // print info
    std::cout<< "settings: " <<ETH << "," <<KMER << "," << DIVIDE_DATA_BY << "," << USED_READS_SIZE << "," <<static_cast<double> (clock) / (num_reads*THREADS_PER_BLOCK) / clockRate << "," <<static_cast<double> (thrust::reduce(thrust::device, dev_counts_ptr, dev_counts_ptr + num_reads * THREADS_PER_BLOCK,0l)) / (num_reads*THREADS_PER_BLOCK) << std::endl;

    // print results
    std::cout << "read index" << ","<< "closest read" << ","<< "edit distance" << std::endl;
    for(int i = 0; i < num_reads; i++){
        std::cout << i << ","<< min_index[i]<< ","<< min_num[i] << std::endl;
    }
    return 0;
}