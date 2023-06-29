#include <iostream>
#include <string>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <random>
using namespace std;
 
int main () {
    vector<int> vec1(5);
    for(int i=0; i<vec1.size(); i++) {
        vec1[i]=i+1;
    }
    // Print old vector
    cout << "Old vector :";
    for(int i=0; i < vec1.size(); i++)
        cout << " " << vec1[i];
    cout << "\n";
	int block_start = 1; //da 1 a size()-1
	int block_finish = 5; //da start +1 a size()
	int n = 3; //tra 1 e finish-start
 
    // rotate function
    rotate(vec1.begin()+block_start, vec1.begin()+n+block_start, vec1.begin()+block_finish); //eventualmente mettere vec1.begin()+1+block_finish
    cout << "New vector after permutation of m continous genes :";
    for (int i=0; i < vec1.size(); i++)
        cout<<" "<<vec1[i];
    cout << "\n\n";
	return 0;


}