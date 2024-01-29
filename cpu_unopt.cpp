#include <iostream>

# if !defined DIVIDE_DATA_BY
# define DIVIDE_DATA_BY 1
# endif

# if !defined UCHAR4_OPTIMIZATION
# define UCHAR4_OPTIMIZATION true
# endif

# if !defined FORCE_UCHAR4_OPTIMIZATION
# define FORCE_UCHAR4_OPTIMIZATION false
# endif
//const int READ_LENGTH = 256;

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
                              min(arr[j] + (s[i - 1] != 'P'), arr[j - 1] + (t[i - 1] != 'P')));
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
                                  min(arr[j/4].x + (s[i - 1] != 'P'), arr[j/4 - 1].w + (t[i - 1] != 'P')));
                        diag = arr[j/4].x;
                        arr[j/4].x = new_val;
                        }
                        break;
                    case 1:
                        {
                        // y
                        int new_val = min(diag + (s[i - 1] != t[j - 1]),
                                  min(arr[j/4].y + (s[i - 1] != 'P'), arr[j/4].x + (t[i - 1] != 'P')));
                        diag = arr[j/4].y;
                        arr[j/4].y = new_val;
                        }
                        break;
                    case 2:
                        {
                        // z
                        int new_val = min(diag + (s[i - 1] != t[j - 1]),
                                  min(arr[j/4].z + (s[i - 1] != 'P'), arr[j/4].y + (t[i - 1] != 'P')));
                        diag = arr[j/4].z;
                        arr[j/4].z = new_val;
                        }
                        break;
                    case 3:
                        {
                        // w
                        int new_val = min(diag + (s[i - 1] != t[j - 1]),
                                  min(arr[j/4].w + (s[i - 1] != 'P'), arr[j/4].z + (t[i - 1] != 'P')));
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



int main() {

    std::string readsStr;
    std::cin >> readsStr;
    const char *reads = readsStr.c_str();
    int reads_length = readsStr.length();
    int num_reads = reads_length / READ_LENGTH;

    int minNum = num_reads + 1;
    int minIdx = -1;
    std::cout<< "settings: " <<-1 << "," <<-1 << "," << -1 << "," << -1 << "," << num_reads << std::endl;

    std::cout << "read index" << ","<< "closest read" << ","<< "edit distance" << std::endl;

    for (int readIdx = 0; readIdx < (num_reads / DIVIDE_DATA_BY); readIdx++){
        for (int tempReadIdx = 0; tempReadIdx < num_reads; tempReadIdx++) {
            if (readIdx != tempReadIdx) {
                int edit_distance = editDistance((reads + READ_LENGTH * (readIdx)), (reads + READ_LENGTH * (tempReadIdx)));
                minIdx = (edit_distance < minNum) ? (tempReadIdx) : minIdx;
                minNum = std::min(edit_distance, minNum);

            }
        }
        std::cout << readIdx << "," << minIdx << "," << minNum << std::endl;
        minNum = num_reads + 1;
        minIdx = -1;
    }

    return 0;
}
