// algo vote.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#include <iostream>
#include "vote2.hpp"
#include <vector>

#include <ctime> // Obligatoire
#include <cstdlib> // Obligatoire


using namespace std;

int main()
{
	srand(static_cast<unsigned>(time(0)));
	int question = 1;

debut:
	int N = 20;
	int n = 15;
	double epsilon = 1.5;
	int debut_test = 4;

	vector<vector<int>> gen;
	for (int i(0); i < N; ++i)
		gen.push_back(ordre_alea(n));
	vector<vector<int>> pref = compresser_2(gen, n);

	if (false) {
		pair<vector<pair<int, double>>, int> resultat_algo = algo_entier(gen, epsilon, n);
		vector<int> resultat_algo_liste;
		resultat_algo_liste.reserve(get<0>(resultat_algo).size());
		for (int i(0); i < get<0>(resultat_algo).size(); ++i)
			resultat_algo_liste.push_back(get<0>(get<0>(resultat_algo)[i]));
		vector<int> resultat_0 = ordre_double(vote(pref, epsilon));

		cout << "compter ordonne : rec    : " << compter_ordonne(pref, resultat_algo_liste, debut_test) << endl;
		cout << "compter ordonne : 0      : " << compter_ordonne(pref, resultat_0, debut_test) << endl;
		cout << "nombre preferences rec>0 : " << comparer_listes(pref, resultat_algo_liste, resultat_0, debut_test) << endl;
		cout << " \t " << endl;
	}
	if (true) {
		vector<int> resultat_0 = ordre_double(vote(pref, 1.5));
		vector<int> resultat_1 = ordre_double(vote(pref, 0.05));
		vector<int> resultat_2 = ordre_double(vote(pref, 0.05,4));
		for (int i(0); i < resultat_0.size(); ++i)
			cout << resultat_0[i] << " ; ";
		cout << endl;
		for (int i(0); i < resultat_1.size(); ++i)
			cout << resultat_1[i] << " ; ";
		cout << endl;
		for (int i(0); i < resultat_2.size(); ++i)
			cout << resultat_2[i] << " ; ";
		cout << endl;
		cout << endl;

	}
	--question;
	if (question == 0) {
		cin >> question;
		if (question < 0)
			question = 0;
		if (question == 0)
			return 0;
	}

	goto debut;
}

// Exécuter le programme : Ctrl+F5 ou menu Déboguer > Exécuter sans débogage
// Déboguer le programme : F5 ou menu Déboguer > Démarrer le débogage

// Astuces pour bien démarrer : 
//   1. Utilisez la fenêtre Explorateur de solutions pour ajouter des fichiers et les gérer.
//   2. Utilisez la fenêtre Team Explorer pour vous connecter au contrôle de code source.
//   3. Utilisez la fenêtre Sortie pour voir la sortie de la génération et d'autres messages.
//   4. Utilisez la fenêtre Liste d'erreurs pour voir les erreurs.
//   5. Accédez à Projet > Ajouter un nouvel élément pour créer des fichiers de code, ou à Projet > Ajouter un élément existant pour ajouter des fichiers de code existants au projet.
//   6. Pour rouvrir ce projet plus tard, accédez à Fichier > Ouvrir > Projet et sélectionnez le fichier .sln.
