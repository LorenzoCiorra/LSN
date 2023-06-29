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
    // Rotate vector left 3 times.
    int rotL=2;
 
    // rotate function
    rotate(vec1.begin()+1, vec1.begin()+1+rotL, vec1.end());
    cout << "New vector after left rotation :";
    for (int i=0; i < vec1.size(); i++)
        cout<<" "<<vec1[i];
    cout << "\n\n";
	return 0;
}