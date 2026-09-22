#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <set>
#include <cstdio>
#include <unistd.h>
using namespace std;

string shell_quote(const string& value) {
    string quoted = "'";
    for (char c : value) quoted += c == '\'' ? "'\\''" : string(1, c);
    return quoted + "'";
}

struct HeaderedInput {
    string path;

    HeaderedInput() {
        char name[] = ".onb_input.XXXXXX";
        const int fd = mkstemp(name);
        if (fd >= 0) {
            close(fd);
            path = name;
        }
    }

    ~HeaderedInput() {
        if (!path.empty()) remove(path.c_str());
    }
};

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

    // ORCA needs a header. Keep that copy in the caller's private working
    // directory so neither successful counting nor failures modify the graph.
    HeaderedInput headered;
    if (headered.path.empty()) {
        cerr << "Error creating headered scratch input\n";
        return 1;
    }
    ifstream orig(file1);
    ofstream fout(headered.path);
    if (!orig || !fout) {
        cerr << "Error preparing headered scratch input\n";
        return 1;
    }

    fout << n << " " << edges << "\n";
    string line;
    while (getline(orig, line)) {
        fout << line << "\n";
    }

    fout.close();
    if (orig.bad() || !fout) {
        cerr << "Error writing headered scratch input\n";
        return 1;
    }
    orig.close();

    // Run orca
    string cmd1 = "./orca " + shell_quote(k) + " " + shell_quote(headered.path)
                  + " " + shell_quote("data_middle.txt") + " > " + shell_quote("smth.txt");
    if (system(cmd1.c_str()) != 0) {
        cerr << "Error running orca\n";
        return 1;
    }

    // Run otb2
	string cmd1_1 = "./otb2 " + shell_quote(k) + " < " + shell_quote("data_middle.txt")
                       + " > " + shell_quote("data1.txt");
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
