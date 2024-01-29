#include <iostream>
#include <cmath>
#include <string.h>
#include <thrust/scan.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>

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
        if (n < length + numElements){
            printf("length is bigger then n\n");
//            exit(1);
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
};

const int N = 115;
const int THREADS_PER_BLOCK = 32;
const int ETH = 200;
const int K = 12;
const int K_CLOSEST = 3;

int main(){

    std::string readsStr;
    std::cin >> readsStr;
    int reads_length = readsStr.length();

    int num_reads = reads_length / N;
    int *d_num_reads; cudaMalloc(&d_num_reads, sizeof(int));
    cudaMemcpy(d_num_reads, &num_reads, sizeof(int), cudaMemcpyHostToDevice);

    const char *reads = readsStr.c_str();
    char *d_reads; cudaMalloc(&d_reads,reads_length * sizeof(char));
    cudaMemcpy(d_reads, reads, reads_length * sizeof(char), cudaMemcpyHostToDevice);

//    Histogram *histograms = (Histogram*) malloc(num_reads * sizeof(Histogram));
    Histogram *d_histograms; cudaMalloc(&d_histograms, num_reads * sizeof(Histogram));

    // will be a list of the minimum edit distance for each read
    int *min_num = (int*) malloc(num_reads * sizeof(int));
    int *d_min_num; cudaMalloc(&d_min_num, num_reads * sizeof(int));

    // will be a list of the index of the minimum edit distance for each read
    int *min_index = (int*) malloc(num_reads * sizeof(int));
    int *d_min_index; cudaMalloc(&d_min_index, num_reads * sizeof(int));

    int *d_read_counts; cudaMalloc(&d_read_counts, std::pow(4,K) * sizeof(int)); cudaMemset(d_read_counts, 0, std::pow(4,K) * sizeof(int));
    int *read_prefix = (int*) malloc(std::pow(4,K) * sizeof(int));
    int *d_read_prefix; cudaMalloc(&d_read_prefix, std::pow(4,K) * sizeof(int));

    thrust::device_ptr<int> dev_ptr(d_read_counts);

    computeHistogram<<<num_reads/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(d_reads,d_histograms, d_num_reads);

    computeReadCounts<<<num_reads/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(d_reads, d_read_counts, d_num_reads);

    int read_counts_sum = thrust::reduce(thrust::device, dev_ptr, dev_ptr + std::pow(4,K),0);// sum read_counts
    thrust::exclusive_scan(thrust::device, dev_ptr, dev_ptr + std::pow(4,K), d_read_prefix);// create a prefix of read_counts

    int *d_read_chunk; cudaMalloc(&d_read_chunk, read_counts_sum * sizeof(int));// create index table with the length of the sum of read_counts
    int *read_chunk = (int*) malloc(read_counts_sum * sizeof(int));
    IndexTable *d_index_table; cudaMalloc(&d_index_table, std::pow(4,K) * sizeof(IndexTable));cudaMemset(d_index_table, 0, std::pow(4,K) * sizeof(IndexTable));
    IndexTable *index_table = (IndexTable*) malloc(std::pow(4,K) * sizeof(IndexTable));
    int *d_tmp_read_counts; cudaMalloc(&d_tmp_read_counts, std::pow(4,K) * sizeof(int)); cudaMemset(d_tmp_read_counts, 0, std::pow(4,K) * sizeof(int));

    computeIndexTable<<<num_reads/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(d_reads, d_read_counts, d_tmp_read_counts, d_read_chunk, d_read_prefix, d_index_table, d_num_reads);
    cudaFree(d_read_counts);
    cudaFree(d_tmp_read_counts);
    findClosest<<<num_reads,THREADS_PER_BLOCK>>>(d_reads,d_min_num, d_min_index, d_histograms,d_index_table,d_read_chunk);

    cudaMemcpy(min_num, d_min_num, num_reads * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(min_index, d_min_index, num_reads * sizeof(int), cudaMemcpyDeviceToHost);

    std::cout << "read index" << ","<< "closest read" << ","<< "edit distance" << std::endl;
    for(int i = 0; i < num_reads; i++){
        std::cout << i << ","<< min_index[i]<< ","<< min_num[i] << std::endl;
    }

    return 0;
}

__global__ void findClosest(const char *reads, int *min_num, int *min_index, Histogram *histograms, IndexTable *index_table, int *read_chunk, int *minClosestDists, int *minClosestDistsIdx){
    __shared__ CyclicBuffer <2 * THREADS_PER_BLOCK, int> samplesBuffer;
    samplesBuffer.reset();  // resets the buffer

    __shared__ int minDistsMaxIdx;
    __shared__ int minClosestDists[K_CLOSEST];
    __shared__ int minClosestDistsIdx[K_CLOSEST];
    __shared__ int distThread[THREADS_PER_BLOCK];
    __shared__ int distThreadIdx[THREADS_PER_BLOCK];
    __shared__ int distances[THREADS_PER_BLOCK]
    __shared__ int count;
    __shared__ char read[N];


    //creates the main read in the memory:
    if (threadIdx.x == 0) {
        minDistsMaxIdx = 0;
        for (int i = 0; i < N; i++) read[i] = reads[(N * blockIdx.x) + i];
    }
    __syncthreads();

    for (int i = 0; i <= N - K; i++) {
        char kmer[K];
        for (int j = 0; j < K; j++) kmer[j] = read[j + i];
        IndexTable kmer_reads = index_table[dnaToDecimal(kmer)];

        for (int j = 0; j < kmer_reads.count; j += THREADS_PER_BLOCK) {
//            printf("%d\n",kmer_reads.count);
            if (j + threadIdx.x < kmer_reads.count) {
                int comp_read = read_chunk[kmer_reads.index + j + threadIdx.x];
                count = 0;
                if (blockIdx.x != comp_read) { //TODO: instead have a set of compared read, put blockIdx.x in it

                    int length = samplesBuffer.length;  //enables to know the current length
                    // Compute histogram diff between our read and read comp_read

                    int diff = abs(histograms[blockIdx.x].numA - histograms[comp_read].numA) +
                               abs(histograms[blockIdx.x].numT - histograms[comp_read].numT) +
                               abs(histograms[blockIdx.x].numG - histograms[comp_read].numG) +
                               abs(histograms[blockIdx.x].numC - histograms[comp_read].numC);
                    // If diff < 2*ETH, add to samples buffer
                    if (diff < 2 * ETH) {
                        int readBufferIdx = atomicAdd(&count,
                                                      1);  //if I did normal add it would be overwritten, returns count before add
                        samplesBuffer.set(readBufferIdx + length, comp_read);
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
                int edit_distance = editDistance(read, (reads + N * (samplesBuffer.get(
                        threadIdx.x)))); // samplesBuffer[threadIdx.x]
                distances[thread.Idx] = (edit_distance < ) ? edit_distance : -1 { // TODO: serialize
                __syncthreads();
                for(int i = 0; i < THREADS_PER_BLOCK; i++){
                    if (distances[i] == -1) { continue; }
                    if (distances[i] < minClosestDists[minDistsMaxIdx]){
                        minClosestDistsIdx[minDistsMaxIdx] = samplesBuffer.get(i);
                        minClosestDists[minDistsMaxIdx] = distances[i];

                        distThread[thread.Idx] = 0;
                        distThreadIdx[thread.Idx] = -1;
                        for(int i = thread.Idx; i < K_CLOSEST; i += THREADS_PER_BLOCK){
                            if (i >= K_CLOSEST){continue;}
                            distThreadIdx[thread.Idx] = (minClosestDists[i] > distThread[thread.Idx]) ? i : distThreadIdx[thread.Idx];
                            distThread[thread.Idx] = max(distThread[thread.Idx], minClosestDists[i]);
                        }
                        __syncthreads();
                        for(int s = THREADS_PER_BLOCK / 2; s > 0; s >>= 1) {
                            if (thread < s) {
                                distThreadIdx[thread] = (distThread[thread + s] > distThread[thread]) ?
                                                        distThreadIdx[thread + s] : distThreadIdx[thread];
                                distThread[thread] = max(distThread[thread], distThread[thread + s]);
                            }
                        }
                        __syncthreads();
                        if(thread == 0) {
                            minDistsMaxIdx = distThreadIdx[0];
                        }
                    }
                }

                __syncthreads();
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
                                             (reads + N * (samplesBuffer.get(threadIdx.x)))); // samplesBuffer[threadIdx.x]
            if (edit_distance < minDistsMax) {
                int previous = atomicMin(&minDistsMax, edit_distance);
                if (edit_distance < previous && edit_distance == minDistsMax) {
                    minDistsMaxIdx = samplesBuffer.get(threadIdx.x);
                }
            }
        }
        __syncthreads();

    }
// TODO: change so we give minClosestDists and minClosestIdx
    if(threadIdx.x == 0){
        min_num[blockIdx.x] = minDistsMax;
        min_index[blockIdx.x] = minDistsMaxIdx;
    }
    __syncthreads();
}




__device__ int editDistance(const char* s, const char* t){

    // The last row
    int arr[N + 1];

    // Initialize arr to be the first row of the DP matrix
    for(int j = 0; j <= N; j++){
        arr[j] = j;
    }

    // Fill the remaining rows
    for(int i = 1; i <= N; i++){
    
        int diag = arr[0];
        arr[0] = i;

        for(int j = 1; j <= N; j++){

            int new_val = min(diag + (s[i - 1] != t[j - 1]), min(arr[j] + (s[i-1] != 'P'), arr[j - 1] + (t[i-1] != 'P')));
            diag = arr[j];
            arr[j] = new_val;

        }

    }

    return arr[N];

}

__global__ void computeHistogram(const char *reads,Histogram *histograms, int *num_reads){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= *num_reads) return;
    histograms[index].numA = 0;
    histograms[index].numT = 0;
    histograms[index].numG = 0;
    histograms[index].numC = 0;
    for(int i = 0; i < N; i++){
        histograms[index].numA += (reads[index * N + i] == 'A');
        histograms[index].numT += (reads[index * N + i] == 'T');
        histograms[index].numG += (reads[index * N + i] == 'G');
        histograms[index].numC += (reads[index * N + i] == 'C');
    }
    
}

__global__ void computeReadCounts(const char *reads, int *read_counts, int *num_reads){
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= *num_reads) return;
    for (int i = 0; i <= N - K; i++){
        char kmer[K];
        for (int j = 0; j < K; j++) kmer[j] = reads[(N * index)+j+i];
        atomicAdd(&read_counts[dnaToDecimal(kmer)], 1);
    }
}

__global__ void computeIndexTable(const char *reads, int *read_counts, int *tmp_read_counts, int *read_chunk, int *prefix, IndexTable *index_table, int *num_reads){
    int read = blockIdx.x * blockDim.x + threadIdx.x;
    if (read >= *num_reads) return;
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

__constant__ int RHO[26] = {0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0};//possibilities
__device__ int dnaToDecimal(const char* dnaSeq) {
    int decimalNum = 0;
    for (int i = K - 1; i >= 0; i--) {
        int nucleotideValue = RHO[dnaSeq[i] - 'A'];
        decimalNum = decimalNum * 4 + nucleotideValue;
    }
    return decimalNum;
}