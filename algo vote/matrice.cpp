#include "matrice.hpp"

using namespace std;

vector<double> multiplication(vector<double> const& X, vector<vector<double>> const& A) {
	int n = X.size();
	vector<double> Y(n, 0);
	for (int i(0); i < n; ++i) {
		double somme = 0;
		for (int j(0); j < n; ++j)
			somme += X[j] * A[j][i];
		Y[i] = somme;
	}

	return Y;
}

vector<vector<double>> multiplication(vector<vector<double>> const& A, vector<vector<double>> const& B) {
	int n = A.size();
	vector<vector<double>> C(n, vector<double>(n, 0.));
	for (int i(0); i < n; ++i)
		for (int j(0); j < n; ++j) {
			double somme = 0;
			for (int k(0); k < n; ++k)
				somme += A[i][k] * B[k][j];
			C[i][j] = somme;
		}
	return C;
}

vector<vector<double>> puissance(vector<vector<double>> matrice, int puissance_n) {
	vector<vector<double>> resultat(matrice.size(), vector<double>(matrice.size(), 0));
	for (int i(0); i < matrice.size(); ++i)
		resultat[i][i] = 1;

	vector<vector<double>> puissance_m = matrice;
	while (puissance_n > 0) {
		if ((puissance_n % 2) == 1)
			resultat = multiplication(resultat, puissance_m);
		puissance_m = multiplication(puissance_m, puissance_m);
		puissance_n = puissance_n / 2;
	}
	return resultat;
}
