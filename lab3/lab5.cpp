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

void insertDefaultVal(double **mac, double *vect){
	mac[0][0] = 1;
	mac[0][1] = 20;
	mac[0][2] = -30;
	mac[0][3] = -4;
	mac[1][0] = 4;
	mac[1][1] = 20;
	mac[1][2] = -6;
	mac[1][3] = 50;
	mac[2][0] = 9;
	mac[2][1] = -18;
	mac[2][2] = 12;
	mac[2][3] = -11;
	mac[3][0] = 16;
	mac[3][1] = -15;
	mac[3][2] = 14;
	mac[3][3] = 130;
	
	vect[0] = 0;
	vect[1] = 114;
	vect[2] = -5;
	vect[3] = 177;
}

void insertValNa5(double **mac, double *vect, int size, double vare){
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			mac[i][j] = 1;
		}
		mac[i][i] = 1 + vare;
		vect[i] = 6 + (2*vare);
	}
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

void calculatingLU(double **A, int size, double **L, double **U){
	copyMacierz(A, U, size);
	//Zerowana kolumne
	for(int k = 0; k < size; k++){
		//wiersz ; ustawianie 1 na przekątnej L i wypełnianie ierwszego wiersza U
		L[k][k] = 1;
		
		for(int i = k + 1; i < size; i++){
			L[i][k] = U[i][k]/U[k][k];
			
			//Kolumny wiersza
			for(int j = k; j < size; j++){
				U[i][j] = U[i][j] - U[k][j] * L[i][k];
			}
		}
	}
}

void calcX(double **L, double **U, double *b, double *w, int size){
	//Korzystamy z właściwości, że det(A) = det(L)*det(U) oraz, że liczenie wyznacznika macierzy trójkątnej, to tylko liczenie przekątnej.
	double **temp = new double *[size];
	allocating(temp, size);
	double *y = new double [size];
	double det;
	
	//ponieważ na przekątnej są same 1
	det = 1;
	macierzOdwrotna(L, size, temp, det);
	copyMacierz(temp, L, size);
	
	for(int i = 0; i < size; i++){
		y[i] = 0;
		for(int j = 0; j < size; j++){
			y[i] += L[i][j] * b[j];
		}
	}
	
	det = detTrojkatna(U, size);
	macierzOdwrotna(U, size, temp, det);
	
	for(int i = 0; i < size; i++){
		w[i] = 0;
		for(int j = 0; j < size; j++){
			w[i] += temp[i][j] * y[j];
		}
	}
	
	copyMacierz(temp, U, size);
	
	
	deleting(temp, size);
	delete(y);
}




int main(){
	int n = 4;
	
	double **A = new double *[n];
	allocating(A, n);
	
	double **U = new double *[n];
	allocating(U, n);
	
	double **L = new double *[n];
	allocating(L, n);
	
	double *b = new double [n];
	
	insertDefaultVal(A, b);
	
//	zerowanie(U, n);
	zerowanie(L, n);

	cout << "Macierz A:" << endl;
	printMacierz(A, 4);
	cout << endl << endl << "Wektor b:" <<endl;
	printVector(b, 4);
	
	calculatingLU(A, 4, L, U);
	
	cout << endl << endl << "Macierz L:" << endl;
	printMacierz(L, 4);
	cout << endl << endl << "Macierz U:" <<endl;
	printMacierz(U, 4);
	
	double *wynik = new double [n];
	
	calcX(L, U, b, wynik, n);
	cout << endl << endl << "Wektor wyniku:" <<endl;
	printVector(wynik, n);
	
//----------------------------------------------------
	cout << endl << endl << "Zadanie na 5" << endl;
	
	double vare = 10e-5;
	double *x = new double [n];
	
	for(int i = 5; i < 21; i++){
		cout << endl << "10^-" << i << endl;
		
		insertValNa5(A, b, n, vare);
		
		calculatingLU(A, n, L, U);
		calcX(L, U, b, x, n);
		checking(A, x, b, wynik, n);
		
		printVectorPrecission(wynik, n, 20);
		vare /= 10;
	}
	
	
	deleting(A, n);
	deleting(U, n);
	deleting(L, n);
	delete(b);
	delete(wynik);
	delete(x);

	return 0;
}
