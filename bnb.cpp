#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>
using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <k> <file1>\n";
        return 1;
    }

    int k = stoi(argv[1]);
    string file1 = argv[2];


    string cmd = "./blant -k " + to_string(k)
                + " -s EBE -n 100000 "
                + file1
                + " > " + "data1.txt";

    int ret = system(cmd.c_str());
    if (ret != 0) {
        cerr << "Error running: " << cmd << "\n";
        return 1;
    }

    // --- Read first column from data1.txt ---
    vector<double> a;

    ifstream f1("data1.txt");

    if (!f1) {
        cerr << "Error: could not open data1.txt\n";
        return 1;
    }
	int xtra;
	double man;

	while (f1 >> man >> xtra) {
		a.push_back(man);
	}

    // --- Normalize so each sums to 100 ---
    double suma = 0;
    for (int i = 0; i < a.size(); i++) suma += a[i];

    if (suma == 0) {
        cerr << "Error: data file sums to zero\n";
        return 1;
    }

    for (int i = 0; i < a.size(); i++) {
        a[i] = (a[i] * 100.0) / suma;
    }

    // --- Output a[i] - b[i] for all 6 ---
    for (int i = 0; i < a.size(); i++) cout << a[i] << " ";
	cout << endl;
    return 0;
}
//  + " 2>/dev/null"
