#include <iostream>
using namespace std;

// A function for listing the absolute paths for each object in a folder
class list_filepaths {

public:
 
void list_filepaths(const char* path) {
    struct dirent* entry;
    DIR* dir = opendir(path);

    if (dir == NULL) {
        return;
    }
    while ((entry = readdir(dir)) != NULL) {
        string ss = entry->d_name;
        if (entry->d_type == DT_REG)
            if (ss[0] != '.') {

                //cout << "C:/dev/MyProg/Test/" << ss << endl;
                ofstream dataset_filepaths;
                dataset_filepaths.open("C:/dev/MyProg/Test/Result/dataset_filepaths.txt", std::ios_base::app);
                dataset_filepaths << "C:/dev/MyProg/Test/" << ss << endl;
            }


    }
    closedir(dir);
}
};