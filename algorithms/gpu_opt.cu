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
//const int READ_LENGTH = 115;

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

__global__ void computeReadCounts(const char *reads, int *read_counts, int num_reads){
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= num_reads) return;
    for (int i = 0; i <= READ_LENGTH - KMER; i++){
        char kmer[KMER];
        for (int j = 0; j < KMER; j++) kmer[j] = reads[(READ_LENGTH * index)+j+i];
        atomicAdd(&read_counts[dnaToDecimal(kmer)], 1);
//        printf("%d: %d\n",dnaToDecimal(kmer),(value+1));

    }
}

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
    // initializing variables
    __shared__ CyclicBuffer <2 * THREADS_PER_BLOCK, int> samplesBuffer;
    samplesBuffer.reset();  // resets the buffer
//    extern __shared__ int usedReads[];
    __shared__ CyclicBuffer <USED_READS_SIZE, int> usedReads;
    usedReads.reset();
    __shared__ int min_distance;
    __shared__ int minIdx;
    __shared__ int count;
    __shared__ char read[READ_LENGTH];


    //creates the main read in the memory:
    if (threadIdx.x == 0) {

//        usedReads[blockIdx.x] = 1;
        min_distance = READ_LENGTH;
        minIdx = -1;
        for (int i = 0; i < READ_LENGTH; i++) read[i] = reads[(READ_LENGTH * blockIdx.x) + i];
    }
    __syncthreads();

    for (int i = 0; i <= READ_LENGTH - KMER; i++) {
        char kmer[KMER];
        for (int j = 0; j < KMER; j++) kmer[j] = read[j + i];
        IndexTable kmer_reads = index_table[dnaToDecimal(kmer)];

        for (int j = 0; j < kmer_reads.count; j += THREADS_PER_BLOCK) {
//            printf("%d\n",kmer_reads.count);
            if (j + threadIdx.x < kmer_reads.count) {

                int comp_read = read_chunk[kmer_reads.index + j + threadIdx.x];
                count = 0;
                if (blockIdx.x != comp_read && !usedReads.contains(comp_read)) {

                    int samplesLength = samplesBuffer.length;  //enables to know the current length
                    // Compute histogram diff between our read and read comp_read

                    int diff = abs(histograms[blockIdx.x].numA - histograms[comp_read].numA) +
                               abs(histograms[blockIdx.x].numT - histograms[comp_read].numT) +
                               abs(histograms[blockIdx.x].numG - histograms[comp_read].numG) +
                               abs(histograms[blockIdx.x].numC - histograms[comp_read].numC);
                    // If diff < 2*ETH, add to samples buffer
                    if (diff < 2 * ETH) {
                        int readBufferIdx = atomicAdd(&count, 1);  //if I did normal add it would be overwritten, returns count before add
                        samplesBuffer.set(readBufferIdx + samplesLength, comp_read);
                    }
                }
            }
            __syncthreads();

            if (threadIdx.x == 0) {
                samplesBuffer.push(count);//enlarges the length of the buffer
//                printf("%d, %d\n",samplesBuffer.length,count);
            }
            __syncthreads();

            if (samplesBuffer.length >= THREADS_PER_BLOCK) {
                counter++;
                int edit_distance = editDistance(read, (reads + READ_LENGTH * (samplesBuffer.get(threadIdx.x)))); // samplesBuffer[threadIdx.x]
                if (edit_distance < min_distance) {
                    int previous = atomicMin(&min_distance, edit_distance);
                    if (edit_distance < previous && edit_distance == min_distance) {
                        minIdx = samplesBuffer.get(threadIdx.x);
                    }
                }
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
                    samplesBuffer.pop(THREADS_PER_BLOCK); //removes all of the checked reads
                }
                __syncthreads();
            }

        }
    }
    // empty samplesBuffer one last time
    if(samplesBuffer.length > 0) {
        if (threadIdx.x < samplesBuffer.length){
            int edit_distance = editDistance(read,
                                             (reads + READ_LENGTH * (samplesBuffer.get(threadIdx.x)))); // samplesBuffer[threadIdx.x]
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
//    printf("The minimum_num of thread %d block %d is: %d\n",threadIdx.x,blockIdx.x,minimum_num);

    if(threadIdx.x == 0){
        min_num[blockIdx.x] = min_distance;
        min_index[blockIdx.x] = minIdx;
//        free(usedReads);

    }
    __syncthreads();
    unsigned long long end = clock64();
    atomicAdd(clock, (end - start));
    counter_sum[index] = counter;
}

int main(){

    std::string readsStr;
    std::cin >> readsStr;
    int reads_length = readsStr.length();

    int num_reads = reads_length / READ_LENGTH;
    const char *reads = readsStr.c_str();
    char *d_reads; cudaMalloc(&d_reads,reads_length * sizeof(char));
    cudaMemcpy(d_reads, reads, reads_length * sizeof(char), cudaMemcpyHostToDevice);

//    Histogram *histograms = (Histogram*) malloc(num_reads * sizeof(Histogram));
    Histogram *d_histograms; cudaMalloc(&d_histograms, num_reads * sizeof(Histogram));

    unsigned long long clock = 0;
    unsigned long long *d_clock; cudaMalloc(&d_clock, sizeof(unsigned long long));

    long *counts = (long*) malloc(num_reads * THREADS_PER_BLOCK * sizeof(long));
    long *d_counts; cudaMalloc(&d_counts, num_reads * THREADS_PER_BLOCK * sizeof(long));

    // will be a list of the minimum edit distance for each read
    int *min_num = (int*) malloc(num_reads * sizeof(int));
    int *d_min_num; cudaMalloc(&d_min_num, num_reads * sizeof(int));

    // will be a list of the index of the minimum edit distance for each read
    int *min_index = (int*) malloc(num_reads * sizeof(int));
    int *d_min_index; cudaMalloc(&d_min_index, num_reads * sizeof(int));

    int *d_read_counts; cudaMalloc(&d_read_counts, std::pow(4,KMER) * sizeof(int)); cudaMemset(d_read_counts, 0, std::pow(4,KMER) * sizeof(int));
    int *read_prefix = (int*) malloc(std::pow(4,KMER) * sizeof(int));
    int *d_read_prefix; cudaMalloc(&d_read_prefix, std::pow(4,KMER) * sizeof(int));

    thrust::device_ptr<int> dev_ptr(d_read_counts);

    computeHistogram<<<num_reads/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(d_reads,d_histograms, num_reads);

    computeReadCounts<<<num_reads/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(d_reads, d_read_counts, num_reads);

    int read_counts_sum = thrust::reduce(thrust::device, dev_ptr, dev_ptr + std::pow(4,KMER),0);// sum read_counts
    thrust::exclusive_scan(thrust::device, dev_ptr, dev_ptr + std::pow(4,KMER), d_read_prefix);// create a prefix of read_counts

    int *d_read_chunk; cudaMalloc(&d_read_chunk, read_counts_sum * sizeof(int));// create index table with the length of the sum of read_counts
    int *read_chunk = (int*) malloc(read_counts_sum * sizeof(int));
    IndexTable *d_index_table; cudaMalloc(&d_index_table, std::pow(4,KMER) * sizeof(IndexTable));cudaMemset(d_index_table, 0, std::pow(4,KMER) * sizeof(IndexTable));
    IndexTable *index_table = (IndexTable*) malloc(std::pow(4,KMER) * sizeof(IndexTable));
    int *d_tmp_read_counts; cudaMalloc(&d_tmp_read_counts, std::pow(4,KMER) * sizeof(int)); cudaMemset(d_tmp_read_counts, 0, std::pow(4,KMER) * sizeof(int));

    computeIndexTable<<<num_reads/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(d_reads, d_read_counts, d_tmp_read_counts, d_read_chunk, d_read_prefix, d_index_table, num_reads);
    findClosest<<<num_reads/DIVIDE_DATA_BY,THREADS_PER_BLOCK>>>(d_reads,d_min_num, d_min_index, d_histograms,d_index_table,d_read_chunk, d_clock,d_counts);

    cudaMemcpy(&clock, d_clock, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
    cudaMemcpy(counts, d_counts, num_reads * THREADS_PER_BLOCK * sizeof(long), cudaMemcpyDeviceToHost);
    thrust::device_ptr<long> dev_counts_ptr(d_counts);

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    int clockRate = prop.clockRate * 1000;// its in KHz so * 1000 to Hz

    std::cout<< "settings: " <<ETH << "," <<KMER << "," << DIVIDE_DATA_BY << "," << USED_READS_SIZE << "," <<static_cast<double> (clock) / (num_reads*THREADS_PER_BLOCK) / clockRate << "," <<static_cast<double> (thrust::reduce(thrust::device, dev_counts_ptr, dev_counts_ptr + num_reads * THREADS_PER_BLOCK,0l)) / (num_reads*THREADS_PER_BLOCK) << std::endl;

    cudaMemcpy(min_num, d_min_num, num_reads * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(min_index, d_min_index, num_reads * sizeof(int), cudaMemcpyDeviceToHost);

    std::cout << "read index" << ","<< "closest read" << ","<< "edit distance" << std::endl;
    for(int i = 0; i < num_reads; i++){
        std::cout << i << ","<< min_index[i]<< ","<< min_num[i] << std::endl;
    }
    return 0;
}