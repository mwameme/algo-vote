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

	debut:
	int N = 20;
	int n = 15;
	double epsilon = 1.5;

	vector<vector<int>> gen;
	for (int i(0); i < N; ++i)
		gen.push_back(ordre_alea(n));
	vector<vector<int>> pref = compresser_2(gen,n);

	vector<pair<int,double>> resultat_algo = algo_entier(gen, epsilon,n);
	vector<int> resultat_algo_liste;
	resultat_algo_liste.reserve(resultat_algo.size());
	for (int i(0); i < resultat_algo.size(); ++i)
		resultat_algo_liste.push_back(get<0>(resultat_algo[i]));
	vector<int> resultat_0 = ordre_double(vote(pref, epsilon));

	cout << "compter ordonne : rec    : " << compter_ordonne(pref, resultat_algo_liste,4) << endl;
	cout << "compter ordonne : 0      : " << compter_ordonne(pref, resultat_0,4) << endl;
	cout << "nombre preferences rec>0 : " << comparer_listes(pref, resultat_algo_liste, resultat_0, 4) << endl;
	long question;
	cin >> question;
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
