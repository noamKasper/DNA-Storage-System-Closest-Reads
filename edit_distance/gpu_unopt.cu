#include <iostream>

//const int READ_LENGTH = 115;
const int THREADS_PER_BLOCK = 32;

# if !defined DIVIDE_DATA_BY
# define DIVIDE_DATA_BY 1
# endif

# if !defined UCHAR4_OPTIMIZATION
# define UCHAR4_OPTIMIZATION true
# endif

# if !defined FORCE_UCHAR4_OPTIMIZATION
# define FORCE_UCHAR4_OPTIMIZATION false
# endif

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
        arr[j].x = j;
        arr[j].y = j+1;
        arr[j].z = j+2;
        arr[j].w = j+3;
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
                    int new_val = min(diag + (s[i - 1] != t[j - 1]),
                              min(arr[j/4].x + (s[i - 1] != 'P'), arr[j/4 - 1].w + (t[j - 1] != 'P')));
                    diag = arr[j/4].x;
                    arr[j/4].x = new_val;
                    }
                    break;
                case 1:
                    {
                    // y
                    int new_val = min(diag + (s[i - 1] != t[j - 1]),
                              min(arr[j/4].y + (s[i - 1] != 'P'), arr[j/4].x + (t[j - 1] != 'P')));
                    diag = arr[j/4].y;
                    arr[j/4].y = new_val;
                    }
                    break;
                case 2:
                    {
                    // z
                    int new_val = min(diag + (s[i - 1] != t[j - 1]),
                              min(arr[j/4].z + (s[i - 1] != 'P'), arr[j/4].y + (t[j - 1] != 'P')));
                    diag = arr[j/4].z;
                    arr[j/4].z = new_val;
                    }
                    break;
                case 3:
                    {
                    // w
                    int new_val = min(diag + (s[i - 1] != t[j - 1]),
                              min(arr[j/4].w + (s[i - 1] != 'P'), arr[j/4].z + (t[j - 1] != 'P')));
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

__global__ void findClosest(const char *reads, int *min_num, int *min_index, int num_reads){

    __shared__ int min_distance;
    __shared__ int minIdx;
    __shared__ char read[READ_LENGTH];

    if (threadIdx.x == 0) {
        min_distance = READ_LENGTH;
        minIdx = -1;
        for (int i = 0; i < READ_LENGTH; i++) read[i] = reads[(READ_LENGTH * blockIdx.x) + i];
    }
    __syncthreads();

    for (int i = 0; i < num_reads; i += THREADS_PER_BLOCK) {
        int index = threadIdx.x + i;
        if (blockIdx.x != index && index < num_reads) {
            int edit_distance = editDistance(read, (reads + READ_LENGTH * index));
            if(edit_distance < min_distance){
                int previous = atomicMin(&min_distance, edit_distance);
                if(edit_distance < previous && edit_distance == min_distance){
                    minIdx = index;
                }
            }
        }
        __syncthreads();
    }

    if(threadIdx.x == 0){
        min_num[blockIdx.x] = min_distance;
        min_index[blockIdx.x] = minIdx;
    }

}

int main(){
//    cudaSetDevice(4);

    std::string readsStr;
    std::cin >> readsStr;
    int reads_length = readsStr.length();

    int num_reads = reads_length / READ_LENGTH;
    const char *reads = readsStr.c_str();
    char *d_reads; cudaMalloc(&d_reads,reads_length * sizeof(char));
    cudaMemcpy(d_reads, reads, reads_length * sizeof(char), cudaMemcpyHostToDevice);

    // will be a list of the minimum edit distance for each read
    int *min_num = (int*) malloc(num_reads * sizeof(int));
    int *d_min_num; cudaMalloc(&d_min_num, num_reads * sizeof(int));

    // will be a list of the index of the minimum edit distance for each read
    int *min_index = (int*) malloc(num_reads * sizeof(int));
    int *d_min_index; cudaMalloc(&d_min_index, num_reads * sizeof(int));
// divide by 4 blocks to split dataset
    findClosest<<<num_reads/DIVIDE_DATA_BY, THREADS_PER_BLOCK>>>(d_reads,d_min_num, d_min_index, num_reads);

    cudaMemcpy(min_num, d_min_num, num_reads * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(min_index, d_min_index, num_reads * sizeof(int), cudaMemcpyDeviceToHost);

    std::cout<< "settings: " <<-1 << "," <<-1 << "," << -1 << "," << -1 << "," << num_reads << std::endl;
    std::cout << "read index" << ","<< "closest read" << ","<< "edit distance" << std::endl;
    for(int i = 0; i < num_reads; i++){
        std::cout << i << ","<< min_index[i]<< ","<< min_num[i] << std::endl;
    }

    return 0;
}