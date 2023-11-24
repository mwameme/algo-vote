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

int compter_ordonne(vector<vector<int>> const& pref, const vector<int> & ordre, int debut) {
	if (debut == 0)
		debut = pref.size();

	int n = pref.size();
	int ordonne = 0;

	for (int i(0); i < debut ; ++i)
		for (int j(i + 1); j < debut; ++j)
			if (pref[ordre[i]][ordre[j]] >= pref[ordre[j]][ordre[i]]) //gauche : i préféré à j
				++ordonne;

	return ordonne; //compte le nombre de pair dans le "bon sens"
};


int comparer_listes(vector<vector<int>> const& pref, const vector<int> & ordre1, const  vector<int>&  ordre2, int n) { // positif si ordre1 est plus souvent préféré à ordre2 ... 
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
	int N = votes.size();
	vector<vector<int>> preferences(n, vector<int>(n, 0));

	for (int i(0); i < N; ++i) { //on parcourt la liste des votes
		vector<bool> vec_fait(n, false);
		for (int j(0); j < votes[i].size(); ++j) {
			vec_fait[votes[i][j]] = true;
			for (int k(j + 1); k < votes[i].size(); ++k) //on regarde j préféré à k. j < k
				preferences[votes[i][j]][votes[i][k]] += 2;
		}
		vector<int> vec_non;
		for (int j(0); j < n; ++j)
			if (!vec_fait[j])
				vec_non.push_back(j);

		for (int j(0); j < votes[i].size(); ++j)
			for (int k(0); k < vec_non.size(); ++k)
				preferences[votes[i][j]][vec_non[k]] += 2;
		for (int j(0); j < vec_non.size(); ++j) {
			for (int k(j + 1); k < vec_non.size(); ++k) {
				preferences[vec_non[j]][vec_non[k]] += 1;
				preferences[vec_non[k]][vec_non[j]] += 1;
			}
		}
	}
	return preferences;
}




vector<double> vote(vector<vector<int>> const& pref, double epsilon_,int repete) {//chaine de markov ... passe d'un tableau de préférence, à une liste ordonnée
	int n = pref.size();
	if (n == 0)
		return vector<double>(0);
	if (n == 1)
		return vector<double>(1, 1.);

	vector<vector<double>> P(n, vector<double>(n, 0)); // matrice de transition
	int N = 0;
	N += pref[0][1];
	N += pref[1][0];


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
	X = multiplication(X, puissance(P, repete));

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
	sort(liste.begin(), liste.end(), [](pair<int, double> a, pair<int, double> b) -> bool {return get<1>(a) > get<1>(b); }); //liste triée en décroissant ...

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
	sort(premiers.begin(), premiers.end(), [](int a, int b) -> bool {return a < b; }); //les premiers sont re-triés
	return premiers;
}

/*
void enlever_liste(vector<vector<int>>& pref, vector<int> const& liste) {//liste croissante. Enleve cette liste de votés du tableau de préférences ...
	for (int i(liste.size() - 1); i >= 0; --i)
		pref.erase(pref.begin() + liste[i]);

	for (int i(liste.size() - 1); i >= 0; --i)
		for (int j(0); j < pref.size(); ++j)
			pref[j].erase(pref[j].begin() + liste[i]);
	return;
	//	return pref;
}
*/


void enlever_liste(vector<vector<int>>& pref, vector<int> const& liste) {//liste non forcement croissante. Enleve cette liste de votés du tableau de préférences ...
	int n = pref.size();
	int m = n - liste.size();
	if (m == 0) {
		pref.resize(0);
		return;
	}

	vector<bool> est_garde(n, true);
	for (int i(0); i < liste.size(); ++i)
		est_garde[liste[i]] = false;

	int j = 0;
	for (int k(0); k < m; ++k, ++j) {
		while (est_garde[j] == false ) {
			++j;
			if (j == n)//normalement n'arrive pas !
				throw std::domain_error("probleme enlever_liste");
		}
		swap(pref[k] , pref[j]);
	}
//	after1:
	pref.resize(m);

	for (int i(0); i < pref.size(); ++i) {
		j = 0;
		for (int k(0); k < m; ++k,++j) {
			while (est_garde[j] ==false) {
				++j;
				if (j == n)//normalement n'arrive pas !
					throw std::domain_error("probleme enlever_liste");
			}
			pref[i][k] = pref[i][j];
//			++j;
		}
//		after2:
		pref[i].resize(m);
	}
	return;
	//	return pref;
}

/*
void decaler_ordre(vector<int>& ordre, vector<int> const& liste) { //liste croissante. la liste "ordre" est modifiee. On avait enleve "liste" du tableau de preferences.
	for (int i(0); i < liste.size(); ++i)
		for (int j(0); j < ordre.size(); ++j)
			ordre[j] = (ordre[j] >= liste[i] ? ordre[j] + 1 : ordre[j]);
	return;
}
*/

void decaler_ordre(int& ordre, vector<int> const& liste) { // idem que decaler_ordre, mais pour un seul entier. liste est supposée triée
	for (int i(0); i < liste.size(); ++i)
		if (ordre >= liste[i])
			++ordre;
		else
			break;
	return;
}


void decaler_ordre(vector<int>& ordre, vector<int> const& liste) { //liste croissante. la liste "ordre" est modifiee. On avait enleve "liste" du tableau de preferences.
	if ((liste.size() == 0) || (ordre.size() ==0))
		return;
	vector<int> save_ordre = ordre_int(ordre); //croissant
	sort(ordre.begin(), ordre.end());

	int j = 0;
	for (int i(0); i < ordre.size(); ++i) {
		if (j == liste.size()) {
			ordre[i] += j;
			continue;
		}
		while (ordre[i] + j >= liste[j]) {
			++j;
			if (j == liste.size())
				break;
		}
		ordre[i] += j;
	}

	vector<int> result(save_ordre.size(),0);
	for (int i(0); i < save_ordre.size(); ++i)
		result[i] = ordre[save_ordre[i]]; //vérifier ?
	ordre = result;

	return;
}

/*
void decaler_ordre_back(vector<int>& ordre, vector<int> const& liste) { //on enleve la liste, même action que pour "enlever_liste" ! 
																		//Les deux listes ne doivent pas avoir d'élément communs
																		//ordre doit être trié (croissant), ainsi que liste.
	int pos_liste = 0;
	for (int i(0); i < ordre.size(); ++i) {
		if (pos_liste == liste.size()) {
			ordre[i] -= pos_liste;
			continue;
		}
		while (liste[pos_liste] < ordre[i]) {
			++pos_liste;
			if (pos_liste == liste.size())
				break;
		}
		ordre[i] -= pos_liste;
	}
	return;
}
*/

void decaler_ordre_back(int& ordre, vector<int> const& liste) { // idem que decaler_ordre_back, mais pour un seul entier. Pas d'overlap.
	int i;
	for (i = 0; i < liste.size(); ++i)
		if (liste[i] > ordre)
			break;
	ordre -= i;
	return;
}



vector<int> get_min(vector<vector<int>> pref, double epsilon) {
	int n = pref.size();
	if (n == 0)
		throw std::domain_error("get_min vide");
	if (n == 1)
		return vector<int>{ 0 };

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

	// deuxieme version, non récurrente ...
	/* 
	vector<pair<int,double>> liste_conjointe;
	liste_conjointe.reserve(resultat.size());
	for(int i(0);i<liste_conjointe.size();++i)
		liste_conjointe.push_back(make_pair(i,-resultat[i]));
	vector<int> derniers = get_premiers(liste_conjointe);
	return derniers;
	*/
}

/*
std::tuple<std::vector<int>, int,vector<int>> get_max_min(std::vector<std::vector<int>>  pref, double epsilon) { //appeler epsilon = 2.
	//0 : a buggé. Renvoit la liste 
	//1 : partiel ... renvoit la liste partielle des max
	//2 : tout est bon : liste entiere des max
	vector<int> vide(0);

	int n = pref.size();
	if (n == 1)
		return std::make_tuple(vector<int>(1, 0), 2,vide);
	vector<int> save_max; //liste partielle des max.
	vector<vector<int>> pref_save = pref;
	vector<double> resultat_save = vote(pref_save, epsilon);

	while (true) {//boucle pour le cas flag_1 (augmente progressivement la liste save_max)
		if (save_max.size() == n)
			return make_tuple(save_max, 2,vide);
		pref = pref_save;
		if ((save_max.size() + pref_save.size()) != n)
			throw std::domain_error("probleme dimensions : get_max_min");

		vector<int> save_max_trie = save_max;
		sort(save_max_trie.begin(), save_max_trie.end());
//		enlever_liste(pref, save_max_trie); //on le fait désormais juste avant les "continue" !!! pour accélérer. c'est pref_save qu'on modifie ainsi.

		vector<int> i_min = get_min(pref, epsilon);
		vector<double> resultat1 = vote(pref, epsilon);

		if (i_min.size() == (n - save_max.size()) ) { //n'arrive pas à séparer grâce à i_min ...
			decaler_ordre(i_min, save_max_trie);
			if (save_max.size() == 0) {
				if(i_min.size() >=2 )
					return std::make_tuple(i_min, 0,vide);
				return make_tuple(i_min, 2, vide);
			}
			//si i_min vaut 0 ... simple
			if (i_min.size() == 1) {
				save_max.push_back(i_min[0]);
				return std::make_tuple(save_max, 2,vide);
			}
			//arrive-t-on à séparer avec pref_save ?
			std::vector<pair<int, double>> liste_conjointe;
			liste_conjointe.reserve(i_min.size());
			for (int i(0); i < i_min.size(); ++i)
				liste_conjointe.push_back(std::make_pair(i_min[i], resultat_save[i_min[i]]));
			std::vector<int> premiers = get_premiers(liste_conjointe); //dans la liste avec min et sans save_max.

			if (premiers.size() >= 2)
				return make_tuple(save_max, 1,premiers);

			int i_max = premiers[0];
			save_max.push_back(i_max);
			decaler_ordre_back(i_max, save_max_trie); //sans save_max_trie
			enlever_liste(pref_save, vector<int>{i_max});

			if ((save_max.size() + pref_save.size()) != n)
				throw std::domain_error("probleme dimensions : get_max_min");
			continue;
		}
		
		enlever_liste(pref, i_min);

		vector<int> resultat;
		int flag;
		vector<int> resultat_premiers;
		std::tie(resultat, flag,resultat_premiers) = get_max_min(pref, epsilon);
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
					return make_tuple(premiers, 0,vide);
				}
				//n'arrive pas à séparer avant d'enlever i_min .... Donc on essaye de séparer avant d'enlever save_max. Dans le cas save_max !=0
				decaler_ordre(premiers, save_max_trie);
				std::vector<pair<int, double>> liste_conjointe_2;
				liste_conjointe_2.reserve(premiers.size());
				for (int i(0); i < premiers.size(); ++i)
					liste_conjointe_2.push_back(std::make_pair(premiers[i], resultat_save[premiers[i]]));
				vector<int> premiers2 = get_premiers(liste_conjointe_2);

				if (premiers2.size() >= 2) //on n'arrive pas à séparer
					return make_tuple(save_max, 1,premiers2);

				//on arrive à séparer ...
				int i_max = premiers2[0];
				save_max.push_back(i_max);
				decaler_ordre_back(i_max, save_max_trie);
				enlever_liste(pref_save, vector<int>{i_max});

				if ((save_max.size() + pref_save.size()) != n)
					throw std::domain_error("probleme dimensions : get_max_min");
				continue;
			}

			int i_max = premiers[0];
			enlever_liste(pref_save, vector<int>{i_max});
			decaler_ordre(i_max, save_max_trie);
			save_max.push_back(i_max);

			if ((save_max.size() + pref_save.size()) != n)
				throw std::domain_error("probleme dimensions : get_max_min");
			continue;
		}
		if (flag == 1) {
			enlever_liste(pref_save, resultat);
			decaler_ordre(resultat, save_max_trie);
			save_max.insert(save_max.end(), resultat.begin(), resultat.end());
			if (save_max.size() == n) //liste entiere des save_max
				return make_tuple(save_max, 2,vide);

			//dans ce cas uniquement, on utilise resutlat_premiers //continue sinon
			decaler_ordre(resultat_premiers, i_min);
			decaler_ordre(resultat_premiers, save_max_trie);
			std::vector<pair<int, double>> liste_conjointe;
			liste_conjointe.reserve(resultat_premiers.size());
			for (int i(0); i < resultat_premiers.size(); ++i)
				liste_conjointe.push_back(std::make_pair(resultat_premiers[i], resultat_save[resultat_premiers[i]]));
			std::vector<int> premiers = get_premiers(liste_conjointe); //dans la liste avec min et sans save_max.
			
			if (premiers.size() >= 2)
				return make_tuple(save_max, 1, premiers);
			//on a réussi à séparer
			int i_max = premiers[0];//i_max est la vraie position

			save_max_trie = save_max;
			sort(save_max_trie.begin(), save_max_trie.end());//on a déjà enlevé qqchose de pref_save ... donc on doit faire en fonction du nouveau save_max_trie !
			save_max.push_back(i_max);

			decaler_ordre_back(i_max, save_max_trie);
			enlever_liste(pref_save, vector<int>{i_max});//ici c'est bon.

			if ((save_max.size() + pref_save.size()) != n)
				throw std::domain_error("probleme dimensions : get_max_min");
			continue;
		}
		if (flag == 2) {
			enlever_liste(pref_save, resultat);
			decaler_ordre(resultat, save_max_trie);
			save_max.insert(save_max.end(), resultat.begin(), resultat.end());
			if (save_max.size() == n)
				return make_tuple(save_max, 2,vide);
			//OK

			if ((save_max.size() + pref_save.size()) != n)
				throw std::domain_error("probleme dimensions : get_max_min");
			continue;
		}
	}
}
*/


std::tuple<std::vector<int>, int> get_max_min(std::vector<std::vector<int>>  pref, double epsilon) { //appeler epsilon = 2.
	//0 : a buggé. Renvoit la liste 
	//1 : partiel ... renvoit la liste partielle des max
//	vector<int> vide(0);

	int n = pref.size();
	if (n == 1)
		return std::make_tuple(vector<int>(1, 0), 1);
	vector<int> save_max; //liste partielle des max.
	vector<vector<int>> pref_save = pref;
	vector<double> resultat_save = vote(pref_save, epsilon);

	while (true) {//boucle pour le cas flag_1 (augmente progressivement la liste save_max)
		if (save_max.size() == n)
			return make_tuple(save_max, 1);
		pref = pref_save;
		if ((save_max.size() + pref_save.size()) != n)
			throw std::domain_error("probleme dimensions : get_max_min");

		vector<int> save_max_trie = save_max;
		sort(save_max_trie.begin(), save_max_trie.end());
		//		enlever_liste(pref, save_max_trie); //on le fait désormais juste avant les "continue" !!! pour accélérer. c'est pref_save qu'on modifie ainsi.

		vector<int> i_min = get_min(pref, epsilon);
//		vector<double> resultat1 = vote(pref, epsilon);

		if (i_min.size() == (n - save_max.size())) { //n'arrive pas à séparer grâce à i_min ...
			decaler_ordre(i_min, save_max_trie);
			if (save_max.size() == 0) {
				if (i_min.size() >= 2)
					return std::make_tuple(i_min, 0);
				return make_tuple(i_min, 1);
			}
			//si i_min vaut 0 ... simple
			if (i_min.size() == 1) {
				save_max.push_back(i_min[0]);
				return std::make_tuple(save_max,1 );
			}
			//arrive-t-on à séparer avec pref_save ?
			std::vector<pair<int, double>> liste_conjointe;
			liste_conjointe.reserve(i_min.size());
			for (int i(0); i < i_min.size(); ++i)
				liste_conjointe.push_back(std::make_pair(i_min[i], resultat_save[i_min[i]]));
			std::vector<int> premiers = get_premiers(liste_conjointe); //dans la liste avec min et sans save_max.

			if (premiers.size() >= 2)
				return make_tuple(save_max, 1);

			int i_max = premiers[0];
			save_max.push_back(i_max);
			decaler_ordre_back(i_max, save_max_trie); //sans save_max_trie
			enlever_liste(pref_save, vector<int>{i_max});

			if ((save_max.size() + pref_save.size()) != n)
				throw std::domain_error("probleme dimensions : get_max_min");
			continue;
		}

		enlever_liste(pref, i_min);

		vector<int> resultat;
		int flag;
		std::tie(resultat, flag) = get_max_min(pref, epsilon);
		decaler_ordre(resultat, i_min);

		if (flag == 0) { //rien n'a fonctionné ... on regarde si avant d'avoir enlevé le i_min on peut trouver un i_max ... et on ré-essaye
			//trier new_resultat selon pref_save (=resultat1), car pref ne marche pas.
			vector<double> resultat1 = vote(pref_save, epsilon);

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
					return make_tuple(premiers, 0);
				}
				//n'arrive pas à séparer avant d'enlever i_min .... Donc on essaye de séparer avant d'enlever save_max. Dans le cas save_max !=0
				decaler_ordre(premiers, save_max_trie);
				std::vector<pair<int, double>> liste_conjointe_2;
				liste_conjointe_2.reserve(premiers.size());
				for (int i(0); i < premiers.size(); ++i)
					liste_conjointe_2.push_back(std::make_pair(premiers[i], resultat_save[premiers[i]]));
				vector<int> premiers2 = get_premiers(liste_conjointe_2);

				if (premiers2.size() >= 2) //on n'arrive pas à séparer
					return make_tuple(save_max, 1);

				//on arrive à séparer ...
				int i_max = premiers2[0];
				save_max.push_back(i_max);
				decaler_ordre_back(i_max, save_max_trie);
				enlever_liste(pref_save, vector<int>{i_max});

				if ((save_max.size() + pref_save.size()) != n)
					throw std::domain_error("probleme dimensions : get_max_min");
				continue;
			}

			int i_max = premiers[0];
			enlever_liste(pref_save, vector<int>{i_max});
			decaler_ordre(i_max, save_max_trie);
			save_max.push_back(i_max);

			if ((save_max.size() + pref_save.size()) != n)
				throw std::domain_error("probleme dimensions : get_max_min");
			continue;
		}
		if (flag == 1) { //on a la liste partielle des max (resultat). On l'ajoute à save_max et on recommence
			enlever_liste(pref_save, resultat);
			decaler_ordre(resultat, save_max_trie);
			save_max.insert(save_max.end(), resultat.begin(), resultat.end());
			if (save_max.size() == n)
				return make_tuple(save_max, 1);
			//OK

			if ((save_max.size() + pref_save.size()) != n)
				throw std::domain_error("probleme dimensions : get_max_min");
			continue;
		}
	}
}

pair<vector<pair<int,double>>,int> algo_entier(vector<vector<int>> tableau, double epsilon,int n) { //a partir de la liste des votes
	vector<vector<int>> pref = compresser_2(tableau,n);
	n = pref.size();
	vector<double> resultat = vote(pref, epsilon);
	vector<int> ordre = ordre_double(resultat);
	int max = 7;
	
//	vector<int> interessant;
//	for (int i(0); (i < ordre.size()) && (i < 10); ++i)
//		interessant.push_back(ordre[i]);
//	sort(interessant.begin(), interessant.end());

	vector<int> ininteressant;//on travaille uniquement avec les 7 premiers (selon le calcul de Markov simple). On les choisit ici
	max = (max< ordre.size()) ? max : ordre.size();
	for (int i(max); i < ordre.size(); ++i)
		ininteressant.push_back(ordre[i]);
	sort(ininteressant.begin(), ininteressant.end());

	enlever_liste(pref, ininteressant);//on garde seulement les interessants
	vector<int> resultat_ordre;//les 7 premiers, apres les get_max_min
	resultat_ordre = get<0>(get_max_min(pref, epsilon));
	decaler_ordre(resultat_ordre, ininteressant);//on redécale

	vector<pair<int, double>> liste_conjointe;//on sauve ces premiers dans la liste de retour
	liste_conjointe.reserve(pref.size());
	for (int i(0); i < resultat_ordre.size(); ++i)
		liste_conjointe.push_back(make_pair(resultat_ordre[i], resultat[resultat_ordre[i]]));

	vector<bool> ordre_autre_bool(n,false);//on regarde les éléments qui ne sont pas dans la liste retour
	for (int i(0); i < resultat_ordre.size(); ++i)
		ordre_autre_bool[resultat_ordre[i]] = true;
	vector<int> ordre_autre;
	for (int i(0); i < ordre_autre_bool.size(); ++i)
		if (!ordre_autre_bool[i])
			ordre_autre.push_back(i);
	//on a ceux qui restent, les vrais "ininteressants"

	vector<double> ordre_autre_double;//on regarde les valeurs de ces éléments.
	ordre_autre_double.reserve(ordre_autre.size());
	for (int i(0); i < ordre_autre.size(); ++i)
		ordre_autre_double.push_back(resultat[ordre_autre[i]]);
	vector<int> ordre_autre_trie = ordre_double(ordre_autre_double);//puis on les trie par ordre décroissant

	for (int i(0); i < ordre_autre_trie.size(); ++i)
		liste_conjointe.push_back(make_pair(ordre_autre[ordre_autre_trie[i]], resultat[ordre_autre[ordre_autre_trie[i]]]));//on les rajoute à la liste de retour ...

	return make_pair(liste_conjointe, resultat_ordre.size());
}