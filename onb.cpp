#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <set>
using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <k> <file1>\n";
        return 1;
    }

    string k = argv[1];
    string file1 = argv[2];

    ifstream fin(file1);
    if (!fin) {
        cerr << "Error opening " << file1 << "\n";
        return 1;
    }

    set<string> nodes;
    string u, v;
    int edges = 0;

    while (fin >> u >> v) {
        nodes.insert(u);
        nodes.insert(v);
        edges++;
    }
    fin.close();

    int n = nodes.size();

    string backup = file1 + ".bak";
    system(("cp " + file1 + " " + backup).c_str());

    // Add #node #edge
    ifstream orig(backup);
    ofstream fout(file1);

    fout << n << " " << edges << "\n";
    string line;
    while (getline(orig, line)) {
        fout << line << "\n";
    }

    orig.close();
    fout.close();

    // Run orca
    string cmd1 = "./orca " + k + " " + file1 + " data_middle.txt" + " > smth.txt";
    if (system(cmd1.c_str()) != 0) {
        cerr << "Error running orca\n";
        return 1;
    }

    system(("mv " + backup + " " + file1).c_str());

    // Run otb2
	string cmd1_1 = "./otb2 " + k + " < data_middle.txt > data1.txt";
    if (system(cmd1_1.c_str()) != 0) {
        cerr << "Error running obt2\n";
        return 1;
    }

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

    double suma = 0;
    for (double x : a) suma += x;

    if (suma == 0) {
        cerr << "Error: data file sums to zero\n";
        return 1;
    }

    for (double &x : a) {
        x = (x * 100.0) / suma;
    }

    for (double x : a) cout << x << " ";
    cout << endl;

    return 0;
}
