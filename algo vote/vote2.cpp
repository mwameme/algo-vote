#include <algorithm>
#include <vector>
#include <functional>

#include <cmath>

#include "vote2.hpp"

using namespace std;

vector<int> ordre_alea(int n) {
	vector<double> liste(n, 0);
	for (int i(0); i < n; ++i)
		liste[i] = ((double)(rand()) / ((double)(RAND_MAX)));
	return ordre_double(liste);
}

int compter_ordonne(vector<vector<int>> const& pref, vector<int> ordre, int debut) {
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


int comparer_listes(vector<vector<int>> const& pref, vector<int> ordre1, vector<int> ordre2, int n) { // positif si ordre1 est plus souvent préféré à ordre2 ... 
	int cumul = 0;

	for (int i(0); i < n; ++i)
		for (int j(i+1); j < n; ++j) {
			int test = pref[ordre1[i]][ordre2[j]] - pref[ordre2[j]][ordre1[i]];
			if (test > 0)
				++cumul;
			else if (test < 0)
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

vector<vector<int>> compresser_2(vector<vector<int>> const& votes,int n) {//les listes ne sont pas entieres ...
/*	vector<int> total;
	for (int i(0); i < votes.size(); ++i)
		for (int j(0); j < votes[i].size(); ++j)
			total.push_back(votes[i][j]);
	sort(total.begin(), total.end());
	auto last = std::unique(total.begin(), total.end());
	// v now holds {1 2 1 3 4 5 4 x x x}, where 'x' is indeterminate
	total.erase(last, total.end());
	int n = total.size();
	*/
	int N = votes.size();
	vector<vector<int>> preferences(n, vector<int>(n, 0));

	for (int i(0); i < N; ++i) { //on parcourt la liste de votes
		vector<bool> vec_fait(false, n);
		for (int j(0); j < votes[i].size(); ++j) {
			vec_fait[votes[i][j]] = true;
			for (int k(j + 1); k < votes[i].size(); ++k) //on regarde j préféré à k. j < k
				preferences[votes[i][j]][votes[i][k]] += 2;
		}
		for (int j(0); j < votes[i].size(); ++j)
			for (int k(0); k < n; ++k)
				if (!vec_fait[k])
					preferences[votes[i][j]][k] += 2;
		for (int j(0); j < n; ++j) {
			if (vec_fait[j])
				continue;
			for (int k(j + 1); k < n; ++k) {
				if (vec_fait[k])
					continue;
				preferences[j][k] += 1;
				preferences[k][j] += 1;
			}
		}
	}
	return preferences;
}




vector<double> vote(vector<vector<int>> const& pref, double epsilon_) {//chaine de markov ... passe d'un tableau de préférence, à une liste ordonnée
	int n = pref.size();

	vector<vector<double>> P(n, vector<double>(n, 0)); // matrice de transition
	int N = 0;
	for (int i(0); i < n; ++i) {
		N += pref[0][i];
		N += pref[i][0];
	}


	double inv_n = 1. / ((double)n);
	double epsilon = epsilon_ / ((double)N);
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

inline bool est_egal(double temp1, double temp2) { //valeur d'égalité
	if (abs(temp1 - temp2) < 0.00001)
		return true;
	return false;
}

vector<int> get_premiers(vector<double> const& resultat, vector<int> const& ordre) { //"ordre" est ordonné. Renvoyer le tableau des premiers ex-aequo (trié)
	vector<pair<int, double>> liste;
	liste.reserve(resultat.size());
	for (int i(0); i < resultat.size(); ++i)
		liste.push_back(make_pair(ordre[i], resultat[i]));
	return get_premiers(liste);
}

vector<int> get_premiers(vector<pair<int,double>> liste) { //"ordre" est ordonné. Renvoyer le tableau des premiers ex-aequo (trié)
	sort(liste.begin(), liste.end(), [](pair<int, double> a, pair<int, double> b) {return get<1>(a) > get<1>(b); }); //liste triée en décroissant ...

	int ordre_temp = std::get<0>(liste[0]);
	vector<int> premiers{ ordre_temp };
//	if (liste.size() == 0)//bug
//		return premiers;
	for (int n(1); n < liste.size(); ++n) {
		if (est_egal(get<1>(liste[n]),get<1>(liste[n-1])))
			premiers.push_back(get<0>(liste[n]));
		else
			break;
	}
	sort(premiers.begin(), premiers.end(), [](int a, int b) {return a < b; }); //les premiers sont re-triés
	return premiers;
}


void enlever_liste(vector<vector<int>>& pref, vector<int> const& liste) {//liste croissante. Enleve cette liste de votés du tableau de préférences ...
	for (int i(liste.size() - 1); i >= 0; --i)
		pref.erase(pref.begin() + liste[i]);

	for (int i(liste.size() - 1); i >= 0; --i)
		for (int j(0); j < pref.size(); ++j)
			pref[j].erase(pref[j].begin() + liste[i]);
	return;
	//	return pref;
}



void decaler_ordre(vector<int>& ordre, vector<int> const& liste) { //liste croissante. la liste "ordre" est modifiee. On avait enleve "liste" du tableau de preferences.
	for (int i(0); i < liste.size(); ++i)
		for (int j(0); j < ordre.size(); ++j)
			ordre[j] = (ordre[j] >= liste[i] ? ordre[j] + 1 : ordre[j]);
	return;
	//	return ordre;
}


vector<int> get_min(vector<vector<int>> pref, double epsilon) {
	int n = pref.size();
	if (n == 1)
		return { 0 };

	vector<double> resultat = vote(pref, epsilon);
	vector<int> positions(resultat.size(), 0);
	for (int i(0); i < positions.size(); ++i)
		positions[i] = i;

	//enlever les premiers ex aequo
	//regarder si il reste un dernier ou non.

	vector<int> i_max = get_premiers(resultat, positions);

	if (i_max.size() == n)
		return i_max;

	enlever_liste(pref, i_max);

	vector<int> i_min = get_min(pref, epsilon);
	decaler_ordre(i_min, i_max);

	return i_min;

	/*
	vector<pair<int,double>> liste_conjointe;
	liste_conjointe.reserve(resultat.size());
	for(int i(0);i<liste_conjointe.size();++i)
		liste_conjointe.push_back(make_pair(i,-resultat[i]));
	vector<int> derniers = get_premiers(liste_conjointe);
	return derniers;
	*/
}


std::pair<std::vector<int>, int> get_max_min(std::vector<std::vector<int>>  pref, double epsilon) { //appeler epsilon = 2.
	//0 : a buggé. Renvoit la liste 
	//1 : partiel ... renvoit la liste partielle des max
	//2 : tout est bon : liste entiere des max

	int n = pref.size();
	if (n == 1)
		return std::make_pair(vector<int>(1, 0), 2);
	vector<int> save_max; //liste partielle des max.
	vector<vector<int>> pref_save = pref;
	vector<double> resultat_save = vote(pref_save, epsilon);

	while (true) {//boucle pour le cas flag_1 (augmente progressivement la liste save_max)
		if (save_max.size() == n)
			return make_pair(save_max, 2);
		pref = pref_save;
		vector<int> save_max_trie = save_max;
		sort(save_max_trie.begin(), save_max_trie.end());
		enlever_liste(pref, save_max_trie);

		vector<int> i_min = get_min(pref, epsilon);
		vector<double> resultat1 = vote(pref, epsilon);

		if (i_min.size() == (n - save_max.size()) ) { //n'arrive pas à séparer grâce à i_min ...
			decaler_ordre(i_min, save_max_trie);
			if (save_max.size() == 0) {
				return std::make_pair(i_min, 0);
			}
			//si i_min vaut 0 ... simple
			if (i_min.size() == 1) {
				save_max.push_back(i_min[0]);
				return std::make_pair(save_max, 2);
			}
			//arrive-t-on à séparer avec pref_save ?
			std::vector<pair<int, double>> liste_conjointe;
			liste_conjointe.reserve(i_min.size());
			for (int i(0); i < i_min.size(); ++i)
				liste_conjointe.push_back(std::make_pair(i_min[i], resultat_save[i_min[i]]));
			std::vector<int> premiers = get_premiers(liste_conjointe); //dans la liste avec min et sans save_max.

			if (premiers.size() >= 2)
				return make_pair(save_max, 1);
			save_max.push_back(premiers[0]);
			continue;
		}
		
		enlever_liste(pref, i_min);

		vector<int> resultat;
		int flag;
		std::tie(resultat, flag) = get_max_min(pref, epsilon);
		decaler_ordre(resultat, i_min);

		if (flag == 0) { //rien n'a fonctionné ... on regarde si avant d'avoir enlevé le i_min on peut trouver un i_max ... et on ré-essaye
			//trier new_resultat selon pref_save (=resultat1), car pref ne marche pas.
			std::vector<pair<int, double>> liste_conjointe;
			liste_conjointe.reserve(resultat.size());
			for (int i(0); i < resultat.size(); ++i)
				liste_conjointe.push_back(std::make_pair(resultat[i], resultat1[resultat[i]]));
			std::vector<int> premiers = get_premiers(liste_conjointe); //dans la liste avec min et sans save_max.

			if (premiers.size() >= 2) {
				if (save_max.size() == 0) {
//					vector<int> positions(pref_save.size(), 0);
//					for (int i(0); i < positions.size(); ++i)
//						positions[i] = i;
					return make_pair(premiers, 0);
				}
				//n'arrive pas à séparer avant d'enlever i_min .... Donc on essaye de séparer avant d'enlever save_max. Dans le cas save_max !=0
				decaler_ordre(premiers, save_max_trie);
				std::vector<pair<int, double>> liste_conjointe_2;
				liste_conjointe_2.reserve(premiers.size());
				for (int i(0); i < premiers.size(); ++i)
					liste_conjointe_2.push_back(std::make_pair(premiers[i], resultat_save[premiers[i]]));
				vector<int> premiers2 = get_premiers(liste_conjointe_2);

				if (premiers2.size() >= 2) //on n'arrive pas à séparer
					return make_pair(save_max, 1);

				//on arrive à séparer ...
				int i_max = premiers2[0];
				save_max.push_back(i_max);
				continue;
			}

			int i_max = premiers[0];
			for (int i(0); i < save_max_trie.size(); ++i)
				i_max = (i_max >= save_max_trie[i] ? i_max + 1 : i_max);

			save_max.push_back(i_max);
			continue;
		}
		if (flag == 1) {
			decaler_ordre(resultat, save_max_trie);
			save_max.insert(save_max.end(), resultat.begin(), resultat.end());
			if (save_max.size() == n)
				return make_pair(save_max, 2);
			continue;
		}
		if (flag == 2) {
			decaler_ordre(resultat, save_max_trie);
			save_max.insert(save_max.end(), resultat.begin(), resultat.end());
			if (save_max.size() == n)
				return make_pair(save_max, 2);
			continue;
		}
	}
}

vector<pair<int,double>> algo_entier(vector<vector<int>> tableau, double epsilon,int n) { //a partir de la liste des votes
	vector<vector<int>> pref = compresser_2(tableau,n);
	vector<double> resultat = vote(pref, epsilon);
	vector<int> ordre = ordre_double(resultat);
	
//	vector<int> interessant;
//	for (int i(0); (i < ordre.size()) && (i < 10); ++i)
//		interessant.push_back(ordre[i]);
//	sort(interessant.begin(), interessant.end());

	vector<int> ininteressant;
	int fin_premiers = (10< ordre.size()) ? 10 : ordre.size();
	for (int i(fin_premiers); i < ordre.size(); ++i)
		ininteressant.push_back(ordre[i]);

	enlever_liste(pref, ininteressant);
	vector<int> resultat_ordre;
	resultat_ordre = get<0>(get_max_min(pref, epsilon));

	decaler_ordre(resultat_ordre, ininteressant);
	vector<pair<int, double>> liste_conjointe;
	liste_conjointe.reserve(resultat_ordre.size());
	for (int i(0); i < liste_conjointe.size(); ++i)
		liste_conjointe.push_back(make_pair(resultat_ordre[i], resultat[resultat_ordre[i]]));
	return liste_conjointe;
}