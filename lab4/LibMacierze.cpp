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
	cout << setprecision(precision);
	for(int i = 0; i < size; i++){
		cout<< setw(precision + 5) << vect[i];
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

void copyMacierz(double **from, double **to, int size){
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			to[i][j] = from[i][j];
		}
	}
}

void copyMacierz(double *from, double *to, int size){
	for(int i = 0; i < size; i++){
		to[i] = from[i];
	}
}

void zerowanie(double **mac, int size){
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++)
			mac[i][j] = 0;
	}
}

void checking(double **A, double *x, double *b, double *w, int size){
	for(int i = 0; i < size; i++){
		w[i] = 0;
		for(int j = 0; j < size; j++){
			w[i] += A[i][j] * x[j];
		}
		w[i] = w[i] - b[i];
	}
	
}

//---------------------------------------------------------------------
double det(double **mac, int size){
	if(size == 1){
		return mac[0][0];
	}
	if(size == 2){
		return mac[0][0] * mac[1][1] - mac[0][1] * mac[1][0];
	}
	if(size == 3){
		double det = 0, temp = 1, temp1 = 1;
		for(int i = 0; i < size; i++){
			for(int j = 0; j < size; j++){
				temp *= mac[(i + j) % size][j];
				temp1 *= mac[(j + i) % size][size - 1 - j];
			}
			det += temp;
			det -= temp1; 
			temp = 1;
			temp1 = 1;
		}
		return det;
	}
	cout << "ERROR nie potrafię obliczyć wyznacznika" << endl;
	exit(4);
}

double detTrojkatna(double **mac, int size){
	double det = 1;
	for(int i = 0; i < size; i++){
		det *= mac[i][i];
	}
	return det;
}

void macierzOdwrotna(double **mac, int size, double **odwrotna, double detMac){
	double **tempMac = new double *[size-1];
	allocating(tempMac, size-1);
	
	//iteratory macierzy tymczasowej
	int ti = 0, tj = 0;
	//Eliminacja linii nie liczonych
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			//liczenie linii
			for(int n = 0; n < size; n++){
				if(n != i){
					for(int m = 0; m < size; m++){
						if(m != j){
							tempMac[ti][tj] = mac[n][m];
							tj++;
							if(tj > size -1){ cout<< "tj out of range" <<endl; exit(2);}
						}
					}
					tj = 0;
					ti++;
					if(tj > size -1){ cout<< "ti out of range" <<endl; exit(2);}
				}
			}
			odwrotna[j][i] =  (1/detMac)*det(tempMac, size-1);
			if((i + j) % 2 == 1)
				odwrotna[j][i] *= -1;
			ti = 0;
		}
	}
	deleting(tempMac, size-1);
	
}

void macierzDiagonalnaOdwrotna(double **A, double **w, int size){
	zerowanie(w, size);
	for(int i = 0; i < size; i++){
		w[i][i] = 1/A[i][i];
	}
}

void mnozenieMacierzy(double **A, double **B, double **w, int size){
	for(int iw = 0; iw < size; iw++){
		for(int jw = 0; jw < size; jw++){
			w[iw][jw] = 0;
			for(int i = 0; i < size; i++){
				w[iw][jw] += A[iw][i] * B[i][jw];
			}
		}
	}
}

void mnozenieMacierzy(double **A, double *B, double *w, int size){
	for(int iw = 0; iw < size; iw++){
		w[iw] = 0;
		for(int i = 0; i < size; i++){
			w[iw] += A[iw][i] * B[i];
		}
	}
}

void mnozenieMacierzy(double **A, double b, double **w, int size){
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			w[i][j] = A[i][j] * b;
		}
	}
}

void mnozenieMacierzy(double *A, double b, double *w, int size){
	for(int i = 0; i < size; i++){
		w[i] = A[i] * b;
	}
}

void dodawanieMacierzy(double **A, double **B, double **w, int size){
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++)
			w[i][j] = A[i][j] + B[i][j];
	}
}

void dodawanieMacierzy(double *A, double *B, double *w, int size){
	for(int i = 0; i < size; i++)
		w[i] = A[i] + B[i];
}

void macierzTrojkatnaDolnaZX(double **L, double *x, double*b, int size)
{
	x[0] = b[0]/L[0][0];
	for(int i = 1; i < size; i++){
		x[i] = b[i];
		for(int j = 0; j < i; j++){
			x[i] -= x[j]*L[i][j];
		}
		x[i] /= L[i][i];
	}
}

void macierzTrojkatnaGornaZX(double **L, double *x, double*b, int size)
{
	x[size - 1] = b[size - 1]/L[size - 1][size - 1];
	for(int i = size - 2; i >= 0; i--){
		x[i] = b[i];
		for(int j = size - 1; j > i; j--){
			x[i] -= x[j]*L[i][j];
		}
		x[i] /= L[i][i];
	}
}
