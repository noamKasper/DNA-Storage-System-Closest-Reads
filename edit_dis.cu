#include <iostream>
#include <cmath>
#include <string.h>

struct Histogram{
    unsigned char numA;
    unsigned char numT;
    unsigned char numC;
    unsigned char numG;
};

struct DynamicVectorBlock{
    /**
     * List of elemneקts
     */
    int x[32];
    /**
     * A pointer to the next node
     */
    DynamicVectorBlock *next;
};

struct DynamicVector{
    /**
     * A pointer to the first node of the linked vector
     */
    DynamicVectorBlock *first_node;
    /**
     * A pointer to the last node of the linked vector
     */
    DynamicVectorBlock *last_node;
    /**
     * The number of elements in the last node
     */
    int n;
    /**
     * If a thread is using the push method it would be true
     */
    int lock;

    /**
     * resets the vector
     */
    __device__ void init(){

        first_node = (DynamicVectorBlock*) malloc(sizeof(DynamicVectorBlock));
        first_node->next = nullptr;
        last_node = first_node;
        n = 0;
        lock = 0;
    }
    /**
     * Add an element to the end of the vector
     * @param element
     */
    __device__ void push(int element){
        while (atomicExch(&lock,1)){}// wait until a thread is done
        last_node->x[atomicAdd(&n, 1)] = element;
        if (n == 32){
            last_node->next = (DynamicVectorBlock*) malloc(sizeof(DynamicVectorBlock));
            last_node = last_node->next;
            last_node->next = nullptr;
        }
        atomicExch(&lock,0);

    }
};

const int NUM_READS = 256;
const int N = 256;
const int THREADS_PER_BLOCK = 32;
const int ETH = 200;
const int K = 12;

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

            int new_val = min(diag + (s[i - 1] != t[j - 1]), min(arr[j], arr[j - 1]) + 1);
            diag = arr[j];
            arr[j] = new_val;

        }

    }

    return arr[N];

}

//__device__ const char* decimalToDna(int num) {
//
//    // Stores DNA representation of number.
//    char DnaNum[K];
//    char DNA[5] = "ACGT";
//
//    for (int i = 0 ;i < K;i++ ){
//        DnaNum[i] = DNA[num % 4];
//        num /= 4;
//    }
//    return DnaNum;
//}
__constant__ int RHO[26] = {0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0};//possibilities
__device__ int dnaToDecimal(const char* dnaSeq) {
    int decimalNum = 0;
    for (int i = K; i >= 0; i--) {
        int nucleotideValue = RHO[dnaSeq[i] - 'A'];
        decimalNum = decimalNum * 4 + nucleotideValue;
    }
    return decimalNum;
}


__global__ void computeHistogram(const char *reads,Histogram *histograms){
    int index = blockIdx.x * blockDim.x + threadIdx.x;
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
            printf("length is bigger then numElements");
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
                printf("length is smaller then numElements");
//                exit(1);
            }
            length -= numElements;
        }
    }
};

__global__ void resetIndexTable(DynamicVector *index_table){
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    printf("hi\n");
    if (index < std::pow(4,K)){
        index_table[index].init();
    }
}

__global__ void computeIndexTable(const char *reads, DynamicVector *index_table){// also might be a situation in which there will be the same read twice in one index in the index_table
    __shared__ char read[N];

//    index_table.init();// how? perhaps make a check in push add a kerenel to init
    if (threadIdx.x == 0){
        for (int i = 0; i < N; i++) read[i] = reads[(N * blockIdx.x) + i];
    }
    for(int i = threadIdx.x; i <= N - K; i +=THREADS_PER_BLOCK){
        char kmer[K];
        for (int j = 0; j < K; j++) kmer[j] = read[j+i];
        index_table[dnaToDecimal(kmer)].push(blockIdx.x);
    }
    __syncthreads();
}

__global__ void findClosest(const char *reads, int *min_num, int *min_index, Histogram *histograms){

    // initializing variables
    __shared__ CyclicBuffer <2 * THREADS_PER_BLOCK, int> samplesBuffer;
    samplesBuffer.reset();  // resets the buffer
    
    __shared__ int min_distance;
    __shared__ int minIdx;
    __shared__ int count;
    __shared__ char read[N];


    //creates the main read in the memory:
    if (threadIdx.x == 0) {
        min_distance = N;
        minIdx = -1;
        for (int i = 0; i < N; i++) read[i] = reads[(N * blockIdx.x) + i];
    }
    __syncthreads();


    for (int i = 0; i < NUM_READS; i += THREADS_PER_BLOCK){
        if (blockIdx.x != (i + threadIdx.x)) {
            count = 0;
            int length = samplesBuffer.length;  //enables to know the current length
            // Compute histogram diff between our read and read i + threadIdx.x

            int diff =  abs(histograms[blockIdx.x].numA - histograms[i + threadIdx.x].numA) +
                    abs(histograms[blockIdx.x].numT - histograms[i + threadIdx.x].numT) +
                    abs(histograms[blockIdx.x].numG - histograms[i + threadIdx.x].numG) +
                    abs(histograms[blockIdx.x].numC - histograms[i + threadIdx.x].numC);
            // If diff < 2*ETH, add to samples buffer
            if (diff < 2*ETH){
                int readBufferIdx = atomicAdd(&count, 1);  //if I did normal add it would be overwritten, returns count before add
                samplesBuffer.set(readBufferIdx + length, i + threadIdx.x);
            }
        }

        __syncthreads();
        if (threadIdx.x == 0){
            samplesBuffer.push(count);//enlarges the length of the buffer
        }
        __syncthreads();

        if(samplesBuffer.length >= THREADS_PER_BLOCK){
            int edit_distance = editDistance(read, (reads + N * (samplesBuffer.get(threadIdx.x)))); // samplesBuffer[threadIdx.x]
            if(edit_distance < min_distance){
                int previous = atomicMin(&min_distance, edit_distance);
                if(edit_distance < previous && edit_distance == min_distance){
                    minIdx = samplesBuffer.get(threadIdx.x);
                }
            }
            __syncthreads();
            if(threadIdx.x == 0){
                samplesBuffer.pop(THREADS_PER_BLOCK); //removes all of the checked reads
            }
            __syncthreads();
        }


    }
    // empty samplesBuffer one last time
    if(samplesBuffer.length > 0) {
        if (threadIdx.x < samplesBuffer.length){
            int edit_distance = editDistance(read,
                                             (reads + N * (samplesBuffer.get(threadIdx.x)))); // samplesBuffer[threadIdx.x]
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
    }
    __syncthreads();

}

int main(){

    std::string readsStr;
    std::cin >> readsStr;
    const char *reads = readsStr.c_str();

    char *d_reads; cudaMalloc(&d_reads,NUM_READS*sizeof(char)*N);
    cudaMemcpy(d_reads, reads, NUM_READS * sizeof(char)*N, cudaMemcpyHostToDevice);

//    Histogram *histograms = (Histogram*) malloc(NUM_READS * sizeof(Histogram));
    Histogram *d_histograms; cudaMalloc(&d_histograms, NUM_READS * sizeof(Histogram));

    DynamicVector *index_table = (DynamicVector*) malloc(std::pow(4,K) * sizeof(DynamicVector));
    DynamicVector *d_index_table; cudaMalloc(&d_index_table, std::pow(4,K) * sizeof(DynamicVector));
    // will be a list of the minimum edit distance for each read
    int *min_num = (int*) malloc(NUM_READS * sizeof(int));
    int *d_min_num; cudaMalloc(&d_min_num, NUM_READS * sizeof(int));

    // will be a list of the index of the minimum edit distance for each read
    int *min_index = (int*) malloc(NUM_READS * sizeof(int));
    int *d_min_index; cudaMalloc(&d_min_index, NUM_READS * sizeof(int));

    resetIndexTable<<<std::pow(4,K)/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(d_index_table);
    cudaDeviceSynchronize();
    computeIndexTable<<<NUM_READS,THREADS_PER_BLOCK>>>(d_reads,d_index_table);

    cudaMemcpy(index_table, d_index_table, std::pow(4,K) * sizeof(DynamicVector), cudaMemcpyDeviceToHost);
    std::cout << std::pow(4,K)/THREADS_PER_BLOCK << std::endl;
    for (int i = 0; i < std::pow(4,K); i++){

        DynamicVectorBlock *current_node = index_table[i].first_node;
        bool printed = false;
        while (current_node != nullptr){
            std::cout << 1 ;
            DynamicVectorBlock temp;
            cudaMemcpy(&temp,current_node,sizeof(DynamicVectorBlock),cudaMemcpyDeviceToHost);
            int n = (temp.next == nullptr)? index_table[i].n : THREADS_PER_BLOCK;
            for (int j = 0; j < n; j++){
                std::cout << temp.x[j] << ", ";
                printed = true;
            }
            current_node = temp.next;
        }
        if (printed){
            std::cout << std::endl;
        }

    }

    computeHistogram<<<NUM_READS/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(d_reads,d_histograms);

//    cudaMemcpy(histograms, d_histograms, NUM_READS * sizeof(Histogram), cudaMemcpyDeviceToHost);

//    for(int i = 0; i < NUM_READS; i++){
//        std::cout << (int)histograms[i].numA << "," << (int)histograms[i].numT << "," <<(int)histograms[i].numG << "," << (int)histograms[i].numC << std::endl;
//    }

    findClosest<<<NUM_READS,THREADS_PER_BLOCK>>>(d_reads,d_min_num, d_min_index, d_histograms);


    cudaMemcpy(min_num, d_min_num, NUM_READS * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(min_index, d_min_index, NUM_READS * sizeof(int), cudaMemcpyDeviceToHost);

//    std::cout << "read index" << ","<< "closest read" << ","<< "edit distance" << std::endl;
//    for(int i = 0; i < NUM_READS; i++){
//        std::cout << i << ","<< min_index[i]<< ","<< min_num[i] << std::endl;
//    }
    return 0;
}