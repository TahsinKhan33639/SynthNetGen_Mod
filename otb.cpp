#include <bits/stdc++.h>
using namespace std;

int main(int argc, char* argv[]) {
    int k = stoi(argv[1]);

	vector<vector<int>> mp(2);
	// orca output to blant graphlet id. k=4,5
	mp[0] = {0,0,0,0,6,6,3,3,8,7,7,7,9,9,10};
	mp[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,14,14,14,10,10,10,10,4,4,15,15,15,22,22,22,22,11,11,11,26,16,16,16,16,17,17,17,17,25,25,23,23,23,23,18,18,27,27,27,19,19,24,24,24,28,28,28,29,29,29,30,30,30,31,31,32,32,33};
    long long a;
	int n;
	vector<long long> count(34, 0);
	vector<bool> vis(34, false);
	for(auto i : mp[k-4]){
		vis[i] = true;
	}
	vis[0] = false;

    int row = 0;
	int pp = 73;
	if (k == 4) pp = 15;
    while (true) {
        for (int j = 0; j < pp; j++) {
            if (!(cin >> a)) {
                goto done;
            }
            count[mp[k-4][j]] += a;
        }
        row++;
    }

done:

	for(int i=0; i<34; i++){
		if (vis[i])
			cout << count[i] << " " << i << endl;
	}
}
