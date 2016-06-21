#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;
class WeightedQuickUnionUF{
public:
	WeightedQuickUnionUF(int N){
		id = new int [N];
		for (int i = 0; i < N; i++)
			id[i] = i;

		sz = new int [N];
		for (int j = 0; j < N; j++)
			sz[j] = 1;
	}


	int find(int p){
		while (p != id[p]){
			p = id [p];
		}
		return p;
	}

	bool connected(int p, int q){
		return find(p) == find(q);
	}

	void unite(int p, int q){
		int i = find(p);
		int j = find(q);
		if (i == j) return;

		if(sz[i] < sz[j]){
			id[i] = j;
			sz[j]  += sz[i];
		}
		else{
			id[j] = i;
			sz[i]  += sz[j];
		}


	}

private:
	int* id;
	int* sz;
};

class Percolation{
public:
	Percolation(int N){
		imaginaryTOP = N*N;
		imaginaryBOT = N*N+1;

		numberOfSingleLine = N;
		opened = new bool* [N];
		for (int row = 0; row < N; row++){
			opened[row] = new bool [N];
			for (int col = 0; col < N; col++)
				opened[row][col] = false;
			}

		uf = new WeightedQuickUnionUF(N*N + 2);
	}


	void open(int i, int j){
		int row = i - 1;
		int col = j - 1;
		if (!opened[row][col])
			opened[row][col] = true;

		if (row == 0){
			uf->unite(imaginaryTOP, col);
			cout << "connect to top!" << endl;
		}

		if (row == numberOfSingleLine-1){
			uf->unite(imaginaryBOT, row*(numberOfSingleLine) + col);
			cout << "connect to bot !" << endl;
		}

		if (row < numberOfSingleLine - 1 && opened[row+1][col] == true){
			uf->unite((row+1)*(numberOfSingleLine) + col, row*(numberOfSingleLine) + col);
			cout << "connect to next row" << endl;
		}

		if (row > 0 && opened[row-1][col] == true){
			uf->unite((row-1)*(numberOfSingleLine) + col, row*(numberOfSingleLine) + col);
			cout << "connect to last row" << endl;
		}

		if (col < numberOfSingleLine - 1 && opened[row][col+1] == true){
			uf->unite((row)*(numberOfSingleLine) + col + 1, row*(numberOfSingleLine) + col);
			cout << "connect to right" << endl;
		}

		if (col > 0 && opened[row][col - 1] == true){
			uf->unite((row)*(numberOfSingleLine) + col - 1, row*(numberOfSingleLine) + col);
			cout << "connect to left" << endl;
		}
	}


	bool isOpen(int i, int j){
		return opened[i-1][j-1];
	}

	bool isFull(int i, int j){
		int row = i - 1;
		int col = j - 1;
		return (uf->connected(imaginaryTOP, row * (numberOfSingleLine) + col));

	}

	bool percolates(){
		return uf->connected(imaginaryTOP, imaginaryBOT);

	}
private:
	WeightedQuickUnionUF* uf;
	bool** opened;
	int imaginaryTOP;
	int imaginaryBOT;
	int numberOfSingleLine;
};

class PercolationStats{
public:
	PercolationStats(int N, int T):experimentNo(T){
		cnt = new double[T];
		srand(time(NULL));
		int colIndex;
		int rowIndex;
		double count;
		double xt;
		for (int i = 0; i < T; i++){
			p = new Percolation(N);
			count = 0.0;
			xt = 0.0;
			while(!p->percolates()){
				colIndex = rand() % N + 1;
				rowIndex = rand() % N + 1;
				cout << rowIndex << " " << colIndex << endl;
				if (!p->isOpen(rowIndex,colIndex)){
					p->open(rowIndex,colIndex);
					count += 1.0;
				}
				xt = count/(N*N);
				cnt[i] = xt;
			}
		}
	}

	double mean(){
		double sum = 0.0;
		for (int i = 0; i < experimentNo; i++){
			sum += cnt[i];
		}
		sum = sum / experimentNo;

		return sum;

	}

	double stddev(){
		double stdsqr = 0.0;
		for (int i = 0; i < experimentNo; i++){
			stdsqr += pow (cnt[i] - this->mean(), 2);
		}
		stdsqr = stdsqr / (experimentNo - 1);
		double stddev = sqrt(stdsqr);

		return stddev;
	}

	double confidenceLo(){
		double CL = 1.96 * this->stddev();
		double sqrtt = sqrt(experimentNo);
		CL = CL/sqrtt;
		CL = this->mean() - CL;

		return CL;


	}

	double confidenceHi(){
		double HL = 1.96 * this->stddev();
		double sqrtt = sqrt(experimentNo);
		HL = HL/sqrtt;
		HL = this->mean() + HL;

		return HL;
	}
private:
	Percolation* p;
	int experimentNo;
	double* cnt;
};

struct CorresspodingPosition{
	string row;
	string col;
};

int main(){
	ifstream file;
	file.open("heart25.txt");

	string totalNumber;
	int N;
	Percolation* p;

	string lineByline;
	int rows;
	int cols;
	CorresspodingPosition cp;
	for (int lineNo = 0; !file.eof(); lineNo++)
	{
		if (lineNo == 0){
			getline(file,totalNumber);
			N = atoi(totalNumber.c_str());
			p = new Percolation(N);
		}
		else{
			getline(file,lineByline);
			istringstream record(lineByline);

			record >> cp.row;
			record >> cp.col;
			rows = atoi(cp.row.c_str());
			cols = atoi(cp.col.c_str());

			p->open(rows,cols);
			if (p->isFull(rows,cols)){
				cout << rows << " " << cols << " " << "FULL" << endl;
			}

		}
	}
	if (p->percolates())
		cout << "percolates" << endl;
	else
		cout << "no percolates" << endl;
	PercolationStats* pStat = new PercolationStats(2, 10000);
	cout << pStat->mean() << endl;
	cout << pStat->confidenceHi()<< ", " << pStat->confidenceLo() << endl;
	cout << pStat->stddev() << endl;

}
