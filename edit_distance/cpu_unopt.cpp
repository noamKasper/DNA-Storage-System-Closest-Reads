#include <iostream>

# if !defined DIVIDE_DATA_BY
# define DIVIDE_DATA_BY 1
# endif


int editDistance(const char* s, const char* t){
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

            int new_val = std::min(diag + (s[i - 1] != t[j - 1]),
                              std::min(arr[j] + (s[i - 1] != 'P'), arr[j - 1] + (t[j - 1] != 'P')));
            diag = arr[j];
            arr[j] = new_val;

        }

    }
    return arr[READ_LENGTH];
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
