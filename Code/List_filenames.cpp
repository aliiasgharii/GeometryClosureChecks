#include <iostream>
using namespace std;
class list_filenames {
// A function for listing the names of objects in a folder
    void list_filenames(const char* path) {
    struct dirent* entry;
    DIR* dir = opendir(path);
    if (dir == NULL) {
        return;
    }
    while ((entry = readdir(dir)) != NULL) {
        string ss = entry->d_name;
        if (entry->d_type == DT_REG)
            if (ss[0] != '.') {
                // cout << ss << endl;
                ofstream dataset_filenames;
                dataset_filenames.open("C:/dev/MyProg/Test/Result/dataset_filenames.txt", std::ios_base::app);
                dataset_filenames << entry->d_name << endl;
            }
    } closedir(dir);
}
};