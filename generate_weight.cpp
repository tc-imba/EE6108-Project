//
// Created by liu on 29/10/2020.
//

#include <fstream>
#include <sstream>
#include <string>
#include <random>

using namespace std;

int main() {
    ifstream fin("facebook.txt");
    ofstream fout("facebook.weighted.txt");

    mt19937 g;
    uniform_int_distribution<int> dis(1, 10);

    istringstream iss;
    string line;
    while (getline(fin, line)) {
        if (line.empty()) continue;
        iss.clear();
        iss.str(line);
        size_t nodeAId, nodeBId;
        iss >> nodeAId >> nodeBId;
        int weight = dis(g);
        fout << nodeAId << " " << nodeBId << " " << weight << endl;
    }
    return 0;
}