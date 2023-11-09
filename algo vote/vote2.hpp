#pragma once

#include <algorithm>
#include <vector>
#include <functional>
#include <stdexcept>

#include "matrice.hpp"
#include <algorithm>
#include <utility>
#include <tuple>

std::vector<int> ordre_alea(int n);

int compter_ordonne(std::vector<std::vector<int>> const& pref, std::vector<int> const& ordre, int debut);//compte le nombre de fois que "ordre" est bien ordonné. Parcourt les paires.

int comparer_listes(std::vector<std::vector<int>> const& pref, std::vector<int> const& ordre, std::vector<int> const& ordre2, int n);// positif si ordre1 est plus souvent préféré à ordre2 ... 

std::vector<std::vector<int>> compresser(std::vector<std::vector<int>> const& votes); //liste de choix de votes. Produit un tableau de préférences

std::vector<std::vector<int>> compresser_2(std::vector<std::vector<int>> const& votes, int n); //idem, mais les listes de votes sont partiels.

std::vector<double> vote(std::vector<std::vector<int>> const& pref, double epsilon); //une seule fois. Passe par la chaine de Markov

bool est_egal(double temp1, double temp2); // quasi-égal : 0.000001

std::vector<int> get_premiers(std::vector<double> const& resultat, std::vector<int> const& ordre); //appelle l'autre version de get_premiers

std::vector<int> get_premiers(std::vector<std::pair<int, double>> liste); //trouve les premiers exeaquo ...

void enlever_liste(std::vector<std::vector<int>>& pref, std::vector<int> const& liste); //liste non forcement croissante. Enleve cette liste de votés du tableau de préférences ...

void decaler_ordre(std::vector<int>& ordre, std::vector<int> const& liste); //liste croissante. la liste "ordre" est modifiee. On avait enleve "liste" du tableau de preferences.

void decaler_ordre(int& ordre, std::vector<int> const& liste);

//void decaler_ordre_back(std::vector<int>& ordre, std::vector<int> const& liste);//liste croissante. on enleve la liste, même action que pour "enlever_liste" ! 
																		//Les deux listes ne doivent pas avoir d'élément communs
																		//ordre doit être trié (croissant).

void decaler_ordre_back(int& ordre, std::vector<int> const& liste); //idem. Pas d'overlap


std::vector<int> get_min(std::vector<std::vector<int>> pref, double epsilon); //enleve le max, itéré, jusqu'à trouver le min. epsilon=2.

//std::tuple<std::vector<int>, int,std::vector<int>> get_max_min(std::vector<std::vector<int>> pref, double epsilon); //enleve le min (de get_min), itéré, jusqu'à obtenir le max.
std::tuple<std::vector<int>, int> get_max_min(std::vector<std::vector<int>> pref, double epsilon);

std::pair<std::vector<std::pair<int, double>>,int> algo_entier(std::vector<std::vector<int>> tableau, double epsilon,int n); //a partir de la liste des votes.
//Prend seulement les 7 premiers pour get_max_min. Les autres sont triés via la chaine de Markov.
//int : taille du get_max_min ...

//ordre 

template<class T> std::vector<int> ordre(std::vector<T> liste_, std::function<bool(T, T)> const& f) {
	int n = liste_.size();
	std::vector<std::pair<int,T>> liste;
	liste.reserve(n);
	for (int i(0); i < n; ++i) {
		liste.push_back( std::make_pair(i, liste_[i]));
	}

	std::sort(liste.begin(), liste.end(), [&](std::pair<int,T> gauche, std::pair<int,T> droite) -> bool {return f(std::get<1>(gauche), std::get<1>(droite)); });
	std::vector<int> resultat;
	resultat.reserve(n);
	for (int i(0); i < n; ++i)
		resultat.push_back(std::get<0>(liste[i]));

	return resultat;
};

inline std::vector<int> ordre_double(std::vector<double> X) { //décroissant
	return ordre<double>(X, [](double gauche, double droite) -> bool {return gauche > droite; });
};

inline std::vector<int> ordre_int(std::vector<int> X) { //croissant
	return ordre<int>(X, [](int gauche, int droite) -> bool {return gauche < droite; });
};

