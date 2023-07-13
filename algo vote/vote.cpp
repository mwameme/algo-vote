#include <algorithm>
#include <vector>
#include <functional>

#include <cmath>

#include "vote.hpp"

using namespace std;

vector<int> ordre_alea(int n) {
	vector<double> liste(n, 0);
	for (int i(0); i < n; ++i)
		liste[i] = ((double)(rand()) / ((double)(RAND_MAX)));
	return ordre_double(liste);
}

int compter_ordonne(vector<vector<int>> const& pref, vector<int> ordre,int debut) {
	if (debut == 0)
		debut = pref.size();

	int n = pref.size();
	int ordonne = 0;

	for (int i(0); i < debut - 1; ++i)
		for (int j(i + 1); j < debut; ++j)
			if (pref[ordre[i]][ordre[j]] >= pref[ordre[j]][ordre[i]]) //gauche : i préféré à j
				++ordonne;

	return ordonne; //compte le nombre de pair dans le "bon sens"
};

int comparer_listes(vector<vector<int>> const& pref, vector<int> ordre, vector<int> ordre2, int n) { // positif si ordre1 est plus souvent préféré à ordre2 ... 
	int cumul = 0;

	for(int i(0);i<n;++i)
		for (int j(0); j < n; ++j) {
			if (i == j)
				continue;
			int test = pref[ordre[i]][ordre2[j]] - pref[ordre2[j]][ordre[i]];
			if (test > 0)
				++cumul;
			else if(test < 0)
				--cumul;
		}

	return cumul;
}

vector<vector<int>> compresser(vector<vector<int>> const& votes) {
	int N = votes.size();
	int n = votes[0].size();
	vector<vector<int>> preferences(n, vector<int>(n, 0));

	for (int i(0); i < N; ++i) //on parcourt la liste de votes
		for (int j(0); j < n - 1; ++j)
			for (int k(j + 1); k < n; ++k) //on regarde j préféré à k. j < k
				preferences[votes[i][j]][votes[i][k]] += 1;
	return preferences;

}


vector<double> vote(vector<vector<int>> const& pref, double epsilon_) {
	int n = pref.size();

	vector<vector<double>> P(n, vector<double>(n, 0)); // matrice de transition
	int N = 0;
	for (int i(0); i < n; ++i) {
		N += pref[0][i];
		N += pref[i][0];
	}


	double inv_n = 1. / ((double) n);
	double epsilon = epsilon_ /((double) N);
	double N_2 = (double)N * .5;

	for (int i(0); i < n; ++i) // on parcourt la matrice de transitions
		for (int j(0); j < n; ++j)
			P[i][j] = inv_n * (1. + epsilon * ((double)pref[j][i] - N_2));

	for (int i(0); i < n; ++i)
		P[i][i] = inv_n;

	double somme = 0; //calculer la stochasticité
	for (int i(0); i < n; ++i) {
		double somme_temp = 0;
		for (int j(0); j < n; ++j)
			somme_temp += P[i][j];
		if (somme_temp > somme)
			somme = somme_temp;
	}

	double somme_inv = 1 / somme;
	for (int i(0); i < n; ++i)
		for (int j(0); j < n; ++j)
			P[i][j] = P[i][j] * somme_inv;

	//modifier les P_ii, toujours pour la stochasticité
	for (int i(0); i < n; ++i) {
		somme = 0.;
		for (int j(0); j < n; ++j) {
			if (j == i)
				continue;
			somme += P[i][j];
		}
		P[i][i] = 1. - somme;
	}


	vector<double> X(n, inv_n); //calculer la proba stationnaire ...
	X = multiplication(X, puissance(P, 1024));

	somme = 0; //au cas où X n'est pas strictement de masse 1.
	for (int i(0); i < n; ++i)
		somme += X[i];
	somme_inv = 1. / somme;
	for (int i(0); i < n; ++i)
		X[i] *= somme_inv;

	return X;
}

bool est_egal(double temp1, double temp2) { //valeur d'égalité
	if (abs(temp1 - temp2) < 0.0001)
		return true;
	return false;
}

vector<int> get_premiers(vector<double> const& resultat, vector<int> const& ordre) {
	int ordre_temp = ordre[0];
	vector<int> premiers {ordre_temp};
	if (ordre.size() == 0)
		return premiers;
	for (int n(1); n < ordre.size(); ++n) {
		if (est_egal(resultat[ordre[n]], resultat[ordre[n - 1]]))
			premiers.push_back(ordre[n]);
		else
			break;
	}
	sort(premiers.begin(), premiers.end(), [](int a, int b) {return a < b; });
	return premiers;
}

void enlever_liste(vector<vector<int>> & pref, vector<int> const& liste) {//liste croissante
	for (int i(liste.size() - 1); i >= 0; --i)
		pref.erase(pref.begin() + liste[i]);
	
	for (int i(liste.size() - 1); i >= 0; --i)
		for (int j(0); j < pref.size(); ++j)
			pref[j].erase(pref[j].begin() + liste[i]);

//	return pref;
}

void decaler_ordre(vector<int> & ordre, vector<int> const& liste) { //liste croissante
	for (int i(0); i < liste.size(); ++i)
		for (int j(0); j < ordre.size(); ++j)
			ordre[j] = (ordre[j] >= liste[i] ? ordre[j] + 1 : ordre[j]);
	
//	return ordre;
}


vector<int> get_min(vector<vector<int>> pref, double epsilon) {
	int n = pref.size();
	if (n == 1)
		return { 0 };

	vector<double> resultat = vote(pref, epsilon);
	vector<int> ordre = ordre_double(resultat);

	//enlever les premiers ex aequo
	//regarder si il reste un dernier ou non.

	vector<int> i_max = get_premiers(resultat,ordre);

	if (i_max.size() == n)
		return i_max;

	enlever_liste(pref, i_max);

	vector<int> i_min = get_min(pref, epsilon);
	decaler_ordre(i_min, i_max);

	return i_min;
}


std::pair<std::vector<int>, int> get_max_min(std::vector<std::vector<int>>  pref, double epsilon) { //appeler epsilon=2
	//0 : a buggé
	//1 : seul le max (un seul)
	//2 : tout est bon : liste entiere

	int n = pref.size();
	if (n == 1)
		return std::make_pair(vector<int>(1, 0), 2);

	std::vector<int> i_min = get_min(pref, epsilon);

	vector<double> resultat1 = vote(pref, epsilon);
	vector<int> ordre1 = ordre_double(vote(pref, epsilon));

	if ( i_min.size() == n) {
//		if ( est_egal(resultat1[ordre1[0]] , resultat1[ordre1[1]]) )
			return std::make_pair(std::vector<int>(0), 0);

		//la suite de ce if ne sert à rien ... à vérifier !
		// cout << "probleme" << endl;
		/*
		int i_max = ordre1[0];

		//on a sélectionné le premier : c'était quand même possible. On essaye d'obtenir ce qu'il reste des premiers : on appelle la même fonction get_max_min() qui s'occupe seulement des premiers.
		enlever_liste(pref, { i_max });

		vector<int> resultat;
		int flag;
		std::tie(resultat, flag) = get_max_min(pref, epsilon);

		if (flag == 0)
			return std::make_pair(std::vector<int>(1, i_max), 1);;
		if (flag == 1) {
			decaler_ordre(resultat, { i_max });
			resultat.insert(resultat.begin(), i_max);
			if (pref.size() + 1 == resultat.size())
				return make_pair(resultat, 2);
			return make_pair(resultat, 1);
		}
		if (flag == 2) {
			decaler_ordre(resultat, { i_max });
			resultat.insert(resultat.begin(), i_max);
			return make_pair(resultat, 2);
		}
		*/
	}

	vector<vector<int>> pref_save = pref;
	enlever_liste(pref, i_min);

	vector<int> resultat;
	int flag;
	std::tie(resultat, flag) = get_max_min(pref, epsilon);

	if (flag == 0) { //rien n'a fonctionné ... on regarde si avant d'avoir enlevé le i_min on peut trouver un i_max ... et on ré-essaye
		if (est_egal(resultat1[ordre1[0]] , resultat1[ordre1[1]]))
			return std::make_pair(std::vector<int>(0), 0);

		int i_max = ordre1[0];

		//on a sélectionné le premier : c'était quand même possible. On essaye d'obtenir ce qu'il reste des premiers : on appelle la même fonction get_max_min() qui s'occupe seulement des premiers.
		enlever_liste(pref_save, { i_max });

		vector<int> resultat2;
		int flag2;
		std::tie(resultat2, flag2) = get_max_min(pref_save, epsilon);

		if (flag2 == 0)
			return std::make_pair(std::vector<int>(1, i_max), 1);
		if (flag2 == 1) {
			decaler_ordre(resultat2, { i_max });
			resultat2.insert(resultat2.begin(), i_max);

			if (resultat2.size() == pref_save.size() + 1) {
//				resultat2.insert(resultat2.end(), i_min.begin(), i_min.end());
				return make_pair(resultat2, 2);
			}
			return make_pair(resultat2, 1);
		}
		if (flag2 == 2) {
			decaler_ordre(resultat2, { i_max });
			resultat2.insert(resultat2.begin(), i_max);
//			resultat2.insert(resultat2.end(), i_min.begin(), i_min.end());

			return make_pair(resultat2, 2);
		}
	}
	if (flag == 1) {
		decaler_ordre(resultat, i_min);
		if (resultat.size() == pref.size()) {
			resultat.insert(resultat.end(), i_min.begin(), i_min.end());
			return make_pair(resultat, 2);
		}
		return make_pair(resultat, 1);
	}
	// flag == 2
	decaler_ordre(resultat, i_min);
	resultat.insert(resultat.end(), i_min.begin(), i_min.end());

	return make_pair(resultat, 2);
}


std::vector<int> boucle_max_min(std::vector<std::vector<int>>  pref, double epsilon) {
	vector<int> resultat;
	int flag;
	std::tie(resultat, flag) = get_max_min(pref, epsilon);

	if (resultat.size() == pref.size())
		return resultat;
	if (flag == 0)
		return vector<int>(0, 0);
	if (flag == 2)
		return resultat;

	enlever_liste(pref, resultat);
	vector<int> resultat_temp = boucle_max_min(pref, epsilon);
	if (resultat_temp.size() == 0)
		return resultat;
	
	decaler_ordre(resultat_temp, resultat);
	resultat.insert(resultat.end(), resultat_temp.begin(), resultat_temp.end());
	return resultat;

}

