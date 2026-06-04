#include <bits/stdc++.h>
using namespace std;

using Vec6 = array<double,6>;

vector<int> wes;

Vec6 add(const Vec6& a, const Vec6& b){
    Vec6 r;
    for(int i=0;i<6;i++) r[i]=a[i]+b[i];
    return r;
}

Vec6 mul(const Vec6& a, double k){
    Vec6 r;
    for(int i=0;i<6;i++) r[i]=a[i]*k;
    return r;
}

Vec6 lerp(const Vec6& P, const Vec6& Q, double u){
    Vec6 R;
    for(int i=0;i<6;i++)
        R[i] = (1-u)*P[i] + u*Q[i];
    return R;
}

//////////////////////////////////////////////////////
// Piecewise displacement
//////////////////////////////////////////////////////

Vec6 displacement(
    double x,
    double ka, double kb,
    const Vec6& A,
    const Vec6& B
){
    Vec6 O{}; // zero vector

    if(x <= ka){
        double u = (ka==0 ? 0 : x/ka);
        return lerp(O, A, u);
    }
    else {
        double u = (x - ka)/(kb - ka);
        return lerp(A, B, u);
    }
}

//////////////////////////////////////////////////////
// Error function
//////////////////////////////////////////////////////

vector<double> build_grid(double kb){
    vector<double> xs;

    for(double x=0; x<=600; x+=24)
        xs.push_back(x);

    for(double x=600; x<=6000; x+=240)
        xs.push_back(x);

    for(double x=6000; x<=60000; x+=2400)
        xs.push_back(x);

    for(double x=60000; x<=kb*10; x+=12000)
        xs.push_back(x);

    return xs;
}

double F6_piece(
    const Vec6& target,
    const Vec6& initial,
    double x,
    double y,
    double ka, double kb,
    const Vec6& A,
    const Vec6& B,
    const Vec6& A2,
    const Vec6& B2
){
    Vec6 v = displacement(x, ka, kb, A, B);
    Vec6 v2 = displacement(y, ka, kb, A2, B2);

    double sum = 0;
    for(int i=0;i<6;i++){
        double val = initial[i] - v[i] - v2[i];
		double d = val - target[i];
		sum += wes[i]*(d*d);
    }

    return sqrt(sum);
}

//////////////////////////////////////////////////////
// Brute search
//////////////////////////////////////////////////////

pair<pair<double,double>, double> minimize_piece(
    const Vec6& target,
    const Vec6& initial,
    double ka, double kb,
    const Vec6& A,
    const Vec6& B,
    const Vec6& A2,
    const Vec6& B2,
	const int e1,
	const int e2
){
    double bestX = 0;
    double bestY = 0;
    double bestF = 200;

	auto xs = build_grid(kb);
	auto ys = build_grid(kb);
	if (e1 == e2 && e1 == 0){
		for(double x : xs){
		    for(double y : ys){
		        double f = F6_piece(target, initial, x, y, ka, kb, A, B, A2, B2);

		        if(f < bestF){
		            bestF = f;
		            bestX = x;
					bestY = y;
		        }
			}
	    }
	}
	else if ((e1 < 0 && e2 > 0) || (e1 > 0 && e2 < 0)){
		for (double x : xs){
			double y = x*(double)e1/(double)e2;
			y = abs(y);
            double f = F6_piece(target, initial, x, y, ka, kb, A, B, A2, B2);

            if(f < bestF){
                bestF = f;
                bestX = x;
                bestY = y;
            }
		}
	}

    return {{bestX, bestY}, bestF};
}

//////////////////////////////////////////////////////
// MAIN
//////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
    int e1 = stoi(argv[1]);
    int e2 = stoi(argv[2]);
    int cm = stoi(argv[3]);

    ifstream file("weights4.txt");
    if (!file) {
        cerr << "Error: Could not open file.\n";
        return 1;
    }

    vector<vector<int>> all_weights(6, vector<int>(6));

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            file >> all_weights[i][j];
        }
    }

    if (cm < 0 || cm >= 6) {
        cerr << "Error: cm out of range.\n";
        return 1;
    }

    wes = all_weights[cm];

    Vec6 target, initial, A, B, C, D, A2, B2, C2, D2;
    double ka = 30000, kb = 300000;

    // target
    for(int i=0;i<6;i++) cin >> initial[i];
    for(int i=0;i<6;i++) cin >> target[i];
	for(int i=0;i<6;i++){
		if(target[i] > 0) target[i] = log(target[i]);
		else target[i] = -15;
		if(initial[i] > 0) initial[i] = log(initial[i]);
		else initial[i] = -15;
	}

    // breakpoints
    //cin >> ka >> kb >> kc;

    // vectors
    for(int i=0;i<6;i++) cin >> A[i];
    for(int i=0;i<6;i++) cin >> B[i];
    for(int i=0;i<6;i++) cin >> A2[i];
    for(int i=0;i<6;i++) cin >> B2[i];

    auto [bestXY, bestF] =
        minimize_piece(target, initial, ka, kb, A, B, A2, B2, e1, e2);
	auto [bestX, bestY] = bestXY;

	cout << fixed << setprecision(0)
	     << bestX << " " << bestY << " ";

	cout << setprecision(10) << bestF << "\n";
}
