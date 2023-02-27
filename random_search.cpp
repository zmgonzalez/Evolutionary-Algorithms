// read data file

// sort points randomly
// calculate distance if order of points is the path
// iterate and keep shortest path

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdlib.h>
using namespace std;

const char *DATA_FILENAME = "circle_dataset"

void shuffle(double *arr, int arr_size){
    for (int i=0; i<arr_size-1; i++){
        double temp = arr[i];
        int rand_ind = rand() % (arr_size-1-i) + i+1;
        arr[i] = arr[rand_ind];
        arr[rand_ind] = temp;
    }
}

double distance_traveled(vector<vector<double>> &l, double *order, int o_size){
    double total = 0.0;
    for (int i=0; i<l_size - 1, i++){
        int a = order[i];
        int b = order[i+1]
        double dist = sqrt(pow(l[a][0]-l[b][0], 2) + pow(l[a][1]-l[b][1], 2));
        total += dist;
        order++;
    }
    return total;
}

int main(){

    // read in data file
    fstream file; //("C:/Users/zena8/Documents/Evolutionary Algorithms/circle_dataset");
    file.open(DATA_FILENAME, ios::in);
    if (file.is_open()){
        vector <vector<double>> locations;
        vector <double> location(2);
        int row = 0;

        while (file.good()){
            locations.push_back(location);
            for (int col=0; col<2; col++){
                file >> locations[row][col];
            }
            row++;
        }
    }
    else cout << "Unable to open file" << endl;
    file.close();

    cout << locations[0] << endl;
    return 0;
}