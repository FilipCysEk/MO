#include <iostream>
#include <math.h>
#include <iomanip>
#include "LibMacierze.cpp"

using namespace std;

//Zmienne globalne 
double static nOperations = 500;
double static minStepSize = 0.000001;
double static deviation = 0.0001;

void insertDefaultVal(double **mac, double *vect, double *x0){
	mac[0][0] = 100;
	mac[0][1] = 1;
	mac[0][2] = -2;
	mac[0][3] = 3;
	mac[1][0] = 4;
	mac[1][1] = 300;
	mac[1][2] = -5;
	mac[1][3] = 6;
	mac[2][0] = 7;
	mac[2][1] = -8;
	mac[2][2] = 400;
	mac[2][3] = 9;
	mac[3][0] = -10;
	mac[3][1] = 11;
	mac[3][2] = -12;
	mac[3][3] = 200;
	
	vect[0] = 395;
	vect[1] = 603;
	vect[2] = -415;
	vect[3] = -606;
	
	for(int i = 0; i < 4; i++)
		x0[i]=1;
}

void podzialNaLDU(double **A, double **L, double **D, double **U, int size){
	zerowanie(L, size);
	zerowanie(D, size);
	zerowanie(U, size);
	
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			if(j < i){
				L[i][j] = A[i][j];
			}
			else if(i == j){
				D[i][j] = A[i][j];
			}
			else
				U[i][j] = A[i][j];
		}
	}
}

double maxKrokSize(double *A, double *B, int size){
	double max = fabs(A[0] - B[0]);
	
	for(int i = 1; i < size; i++)
		if(max < fabs(A[i] - B[i]))
			max = fabs(A[i] - B[i]);
			
	return max;
}

double maxDeviation(double **A, double *x, double *b, int size){
	double *temp = new double [size];
	mnozenieMacierzy(A, x, temp, size);
	double max = fabs(temp[0] - b[0]);
	
	for(int i = 1; i < size; i++)
		if(temp[i] > max)
			max = fabs(temp[i] - b[i]);
	
	delete(temp);
	return max;
}

void metodaJacobiego(double **A, double **L, double **D, double **U, double *x, double *x0, double *b, int size){
	double *xm1 = new double [size];
	copyMacierz(x0, xm1, size);
	
	double **M = new double *[size];
	allocating(M, size);
	
	double **Dodwr = new double *[size];
	allocating(Dodwr, size);
	
	double **Temp = new double *[size];
	allocating(Temp, size);
	
	double *C = new double [size];
	
	macierzDiagonalnaOdwrotna(D, Dodwr, size);
	
	dodawanieMacierzy(L, U, Temp, size);
	mnozenieMacierzy(Dodwr, Temp, M, size);
	mnozenieMacierzy(M, -1, M, size);
	
	mnozenieMacierzy(Dodwr, b, C, size);
	
	bool endByOperations = 1;
	
	double estymator, residium;
	
	for(int i = 0; i < nOperations; i++){
		mnozenieMacierzy(M, xm1, x, size);
		//mnozenieMacierzy(x, -1, x, size);
		dodawanieMacierzy(x, C, x, size);
		
		//Sprawdzanie, czy wielkość kroku jet odpowiednia
		estymator = maxKrokSize(x, xm1, size);
		residium = maxDeviation(A, x, b, size);
		
		cout<<"Iteracja " << i << endl;
		printVector(x, size);
		cout<< endl << "Estymator: " << estymator << "\t\tResidium: " << residium << endl << endl;
		
		if(estymator < minStepSize){
			endByOperations = 0;
			cout << "Zakonczono, poniewaz krok byl zbyt maly!!" << endl;
			break;
		}
		
		if(residium < deviation){
			endByOperations = 0;
			cout << "Zakonczono, poniewaz znalezlismy sie zbyt blisko wyniku!!" << endl;
			break;
		}
			
		for(int j = 0; j < size; j++)
			xm1[j] = x[j];
	}
	
	if(endByOperations)
		cout << "Zakonczono, poniewaz przekroczono ilosc operacji!!" << endl;
	
	deleting(M, size);
	deleting(Dodwr, size);
	deleting(Temp, size);
	delete(xm1);
}

void metodaGaussaSeidela(double **A, double **L, double **U, double **D, double *x0, double *x, double *b, int size){
	double *xm1 = new double [size];
	copyMacierz(x0, xm1, size);
	
	double **M = new double *[size];
	allocating(M, size);
	
	double **Dodwr = new double *[size];
	allocating(Dodwr, size);
	
	double **Temp = new double *[size];
	allocating(Temp, size);
	
	double *C = new double [size];
	
	dodawanieMacierzy(L, D, Temp, size);
	
	double det_val = detTrojkatna(Temp, size);
	macierzOdwrotna(Temp, size, M, det_val);
	
	mnozenieMacierzy(M, b, C, size);
	
	mnozenieMacierzy(M, U, Temp, size);
	mnozenieMacierzy(Temp, -1, M, size);
	
	bool endByOperations = 1;
	
	double estymator, residium;
	
	for(int i = 0; i < nOperations; i++){
		mnozenieMacierzy(M, xm1, x, size);
		dodawanieMacierzy(x, C, x, size);
		
		//Sprawdzanie, czy wielkość kroku jet odpowiednia
		estymator = maxKrokSize(x, xm1, size);
		residium = maxDeviation(A, x, b, size);
		
		cout<<"Iteracja " << i << endl;
		printVector(x, size);
		cout<< endl << "Estymator: " << estymator << "\t\tResidium: " << residium << endl << endl;
		
		if(estymator < minStepSize){
			endByOperations = 0;
			cout << "Zakonczono, poniewaz krok byl zbyt maly!!" << endl;
			break;
		}
		
		if(residium < deviation){
			endByOperations = 0;
			cout << "Zakonczono, poniewaz znalezlismy sie zbyt blisko wyniku!!" << endl;
			break;
		}
			
		for(int j = 0; j < size; j++)
			xm1[j] = x[j];
	}
	
	if(endByOperations)
		cout << "Zakonczono, poniewaz przekroczono ilosc operacji!!" << endl;
		
	deleting(M, size);
	deleting(Dodwr, size);
	deleting(Temp, size);
	delete(xm1);
}

void metodaSOR(double **A, double **L, double **U, double **D, double *x0, double *x, double *b, double w, int size){
	double *xm1 = new double [size];
	copyMacierz(x0, xm1, size);
	
	double **M = new double *[size];
	allocating(M, size);
	
	double **Dodwr = new double *[size];
	allocating(Dodwr, size);
	
	double **Temp = new double *[size];
	allocating(Temp, size);
	
	double **Temp1 = new double *[size];
	allocating(Temp1, size);
	
	double *C = new double [size];
	
	double det_val, Temp_val;
	
	mnozenieMacierzy(D, 1./w, Temp, size);
	dodawanieMacierzy(L, Temp, Temp, size);
	
	det_val = detTrojkatna(Temp, size);
	macierzOdwrotna(Temp, size, M, det_val);
	
	mnozenieMacierzy(M, b, C, size);
	
	Temp_val = 1 - 1./w;
	
	mnozenieMacierzy(D, Temp_val, Temp, size);
	dodawanieMacierzy(Temp, U, Temp1, size);
	mnozenieMacierzy(M, Temp1, Temp, size);
	mnozenieMacierzy(Temp, -1, M, size);
	
	bool endByOperations = 1;
	
	double estymator, residium;
	
	for(int i = 0; i < nOperations; i++){
		mnozenieMacierzy(M, xm1, x, size);
		dodawanieMacierzy(x, C, x, size);
		
		//Sprawdzanie, czy wielkość kroku jet odpowiednia
		estymator = maxKrokSize(x, xm1, size);
		residium = maxDeviation(A, x, b, size);
		
		cout<<"Iteracja " << i << endl;
		printVector(x, size);
		cout<< endl << "Estymator: " << estymator << "\t\tResidium: " << residium << endl << endl;
		
		if(estymator < minStepSize){
			endByOperations = 0;
			cout << "Zakonczono, poniewaz krok byl zbyt maly!!" << endl;
			break;
		}
		
		if(residium < deviation){
			endByOperations = 0;
			cout << "Zakonczono, poniewaz znalezlismy sie zbyt blisko wyniku!!" << endl;
			break;
		}
			
		for(int j = 0; j < size; j++)
			xm1[j] = x[j];
	}
	
	if(endByOperations)
		cout << "Zakonczono, poniewaz przekroczono ilosc operacji!!" << endl;
		
	deleting(M, size);
	deleting(Dodwr, size);
	deleting(Temp, size);
	deleting(Temp1, size);
	delete(xm1);
	cout<< "???"<<endl;
}



int main(){
	int n = 4;
	
	double **A = new double *[n];
	allocating(A, n);
	
	double **L = new double *[n];
	allocating(L, n);
	
	double **D = new double *[n];
	allocating(D, n);
	
	double **U = new double *[n];
	allocating(U, n);
	
	double **Dodw = new double *[n];
	allocating(Dodw, n);
	//Temp
	
	
	double *b = new double [n];
	double *x0 = new double [n];
	double *x = new double [n];
	
	insertDefaultVal(A, b, x0);
	podzialNaLDU(A, L, D, U, n);
	
	cout<< "A" << endl;
	printMacierz(A, n);
	cout<< "L" << endl;
	printMacierz(L, n);
	cout<< "D" << endl;
	printMacierz(D, n);
	cout<< "U" << endl;
	printMacierz(U, n);
	
	macierzDiagonalnaOdwrotna(D, Dodw, n);
	cout<< "D odwrotna" << endl;
	printMacierz(Dodw, n);
	
	metodaJacobiego(A, L, D, U, x, x0, b, n);
	
	printVector(x, n);
	
	cout<< endl << endl << endl << "Metoda Gaussa -Seidela" << endl;
	
	metodaGaussaSeidela(A, L, U, D, x0, x, b, n);
	
	printVector(x, n);
	
	
	cout<< endl << endl << endl << "Metoda SOR" << endl;
	
	metodaSOR(A, L, U, D, x0, x, b, 1/2., n);
	
	printVector(x, n);
	
	
	deleting(A, n);
	deleting(L, n);
	deleting(D, n);
	deleting(U, n);
	deleting(Dodw, n);
	delete(b);
	delete(x0);
	delete(x);

	return 0;
}
