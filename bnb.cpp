#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>
using namespace std;

string shell_quote(const string& value) {
    string quoted = "'";
    for (char c : value) quoted += c == '\'' ? "'\\''" : string(1, c);
    return quoted + "'";
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <k> <file1>\n";
        return 1;
    }

    int k = stoi(argv[1]);
    string file1 = argv[2];
    const char* output_env = getenv("BNB_OUTPUT_FILE");
    const string output_file = output_env && *output_env ? output_env : "data1.txt";

    // The caller supplies a shell-quoted graph argument. BLANT keeps its
    // resource-directory cwd while its output can live in private scratch.
    string cmd = "./blant -k " + to_string(k)
                + " -s EBE -n 100000 "
                + file1
                + " > " + shell_quote(output_file)  + " 2>/dev/null";

    int ret = system(cmd.c_str());
    if (ret != 0) {
        cerr << "Error running: " << cmd << "\n";
        return 1;
    }

    // --- Read first column from data1.txt ---
    vector<double> a;

    ifstream f1(output_file);

    if (!f1) {
        cerr << "Error: could not open " << output_file << "\n";
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
