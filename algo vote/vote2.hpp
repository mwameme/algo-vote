#pragma once

#include <algorithm>
#include <vector>
#include <functional>

#include "matrice.hpp"
#include <algorithm>
#include <utility>
#include <tuple>

std::vector<int> ordre_alea(int n);

int compter_ordonne(std::vector<std::vector<int>> const& pref, std::vector<int> ordre, int debut);

//int comparer_listes(std::vector<std::vector<int>> const& pref, std::vector<int> ordre, std::vector<int> ordre2, int n);// positif si ordre1 est plus souvent préféré à ordre2 ... 

std::vector<std::vector<int>> compresser(std::vector<std::vector<int>> const& votes);

std::vector<double> vote(std::vector<std::vector<int>> const& pref, double epsilon); //une seule fois 

bool est_egal(double temp1, double temp2);

std::vector<int> get_premiers(std::vector<double> const& resultat, std::vector<int> const& ordre);

std::vector<int> get_premiers(std::vector<std::pair<int, double>> liste);

void enlever_liste(std::vector<std::vector<int>>& pref, std::vector<int> const& liste);

void decaler_ordre(std::vector<int>& ordre, std::vector<int> const& liste);

std::vector<int> get_min(std::vector<std::vector<int>> votes, double epsilon); //enleve le max, itéré, jusqu'à trouver le min. epsilon=2.

std::pair<std::vector<int>, int> get_max_min(std::vector<std::vector<int>> pref, double epsilon); //enleve le min (de get_min), itéré, jusqu'à obtenir le max.

std::vector<std::pair<int, double>> algo_entier(std::vector<std::vector<int>> tableau, double epsilon); //a partir de la liste des votes


//ordre 

template<class T> std::vector<int> ordre(std::vector<T> liste_, std::function<bool(T, T)> const& f) {
	int n = liste_.size();
	std::vector<std::pair<int,T>> liste;
	liste.reserve(n);
	for (int i(0); i < n; ++i) {
		liste.push_back( std::make_pair(i, liste_[i]));
	}

	std::sort(liste.begin(), liste.end(), [&](std::pair<int,double> gauche, std::pair<int,double> droite) {return f(std::get<1>(gauche), std::get<1>(droite)); });
	std::vector<int> resultat;
	resultat.reserve(n);
	for (int i(0); i < n; ++i)
		resultat.push_back(std::get<0>(liste[i]));

	return resultat;
};

inline std::vector<int> ordre_double(std::vector<double> X) { //décroissant
	return ordre<double>(X, [](double gauche, double droite) {return gauche >= droite; });
};
