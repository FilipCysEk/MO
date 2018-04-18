#include <iostream>
#include <math.h>
#include <iomanip>

using namespace std;

void printMacierz(double *mac, int size){
	for(int i = 0; i < size; i++){
		//cout<<"";
		for(int j = 0; j < size; j++){
			cout << mac[size*i+j] << "\t";
		}
		cout << endl;
	}
}

void printMacierz(double **mac, int size){
	cout << setprecision(4) << fixed;
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			cout << setw(11) << mac[i][j];
		}
		cout << endl;
	}
}

void printMacierz2(double **mac, int size){
	cout << "{";
	for(int i = 0; i < size; i++){
		cout<<"{";
		cout << mac[i][0];
		for(int j = 1; j < size; j++){
			
			cout << ", " <<mac[i][j];
		}
		cout << "},"<< endl;
	}
	cout << "}" << endl <<endl;
}

void printVector(double *vect, int size){
	cout << setprecision(4) << fixed;
	for(int i = 0; i < size; i++){
		cout<< setw(11) << vect[i];
	}
}

void printVectorPrecission(double *vect, int size, int precision){
	cout << setprecision(precision) << fixed;
	for(int i = 0; i < size; i++){
		cout<< vect[i] << "\t";
	}
}

void allocating(double **mac, int size){
	for(int i = 0; i < size; i++)
		mac[i] = new double[size];
}

void deleting(double **mac, int size){
	for(int i = 0; i < size; i++)
		delete(mac[i]);
		
	delete(mac);
}

void insertDefaultVal(double **mac, double *vect){
	mac[0][0] = 30;
	mac[0][1] = 2/3.;
	mac[1][0] = 3/4.;
	mac[1][1] = 20;
	mac[1][2] = 5/6.;
	mac[2][1] = 7/8.;
	mac[2][2] = 10;
	mac[2][3] = 9/10.;
	mac[3][2] = 11/12.;
	mac[3][3] = 10;
	mac[3][4] = 13/14.;
	mac[4][3] = 15/16.;
	mac[4][4] = 20;
	mac[4][5] = 17/18.;
	mac[5][4] = 19/20.;
	mac[5][5] = 30;
	
	vect[0] = 94/3.;
	vect[1] = 173/4.;
	vect[2] = 581/20.;
	vect[3] = -815/28.;
	vect[4] = -6301/144.;
	vect[5] = -319/10.;
}

void copyMacierz(double **from, double **to, int size){
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			to[i][j] = from[i][j];
		}
	}
}

void zerowanie(double **mac, int size){
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++)
			mac[i][j] = 0;
	}
}
//-----------------------------------------------------

void thomasA(double **A, double *b, int size){
	for(int i = 1; i < size; i++){
		A[i][i] = A[i][i] - A[i][i - 1] / A[i - 1][i - 1] * A[i - 1][i];
		b[i] = b[i] - A[i][i - 1] / A[i - 1][i - 1] * b[i - 1];
		A[i][i - 1] = 0;
	}
}

void thomasX(double **A, double *b, double *x, int size){
	x[size - 1] = b[size - 1] / A[size - 1][size - 1];
	
	for(int i = size -2; i >= 0; i--){
		x[i] = (b[i] - A[i][i + 1] * x[i + 1]) / A[i][i];
	}
}




int main(){
	int n = 6;
	
	double **A = new double *[n];
	allocating(A, n);
	
	double *b = new double [n];
	
	zerowanie(A, n);
	insertDefaultVal(A, b);

	cout << "Macierz A:" << endl;
	printMacierz(A, n);
	cout << endl << endl << "Wektor b:" <<endl;
	printVector(b, n);
	
	thomasA(A, b, n);
	
	
	double *x = new double [n];
	
	thomasX(A, b, x, n);
	cout << endl << endl << "Wektor wyniku:" <<endl;
	printVector(x, n);
	
	deleting(A, n);
	delete(b);
	delete(x);
	
	return 0;
}
