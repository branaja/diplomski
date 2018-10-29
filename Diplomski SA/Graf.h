#pragma once
using namespace std;
#include<vector>
#include<set>
#include<string>
#include<map>
#include<math.h>
#include<algorithm>
#include <utility> 
#include<boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost\serialization\map.hpp>
#include <boost\archive\text_oarchive.hpp>
#include <boost\archive\text_iarchive.hpp>
#include <fstream>
#include <boost\serialization\vector.hpp>

template < typename UndirectedGraph >
class Graf
{
public:
	Graf();
	virtual ~Graf();
	std::vector<int> simuliranoKaljenje(std::vector<int> zadatak, double temp, int maxIteracija);
	std::vector<int> dijametarAlgoritam(std::vector<int> zadatak);
	std::vector<int> SteinerStabloAlgoritam(UndirectedGraph* graf, std::vector<int> trazeniVrhovi);
	std::vector<int> pokrivacSteinerAlgoritam(std::vector<int> zadatak);
	std::vector<int> poboljsaniSteinerAlgoritam(UndirectedGraph* graf, std::vector<int> zadatak);
	std::vector<int> pcelinjiAlgoritam(std::vector<int> zadatak, int maxIter, int n, int m, int e, int nep, int nsp, int ngh, double alfa);

	UndirectedGraph G;
	UndirectedGraph prosireniGraf;
	bool postojiLiProsireni;
	std::vector<string> vecVjestina;
	std::vector<string> vecImena;
	double parametarHladenja;
	int maxIteracija;
	int pocetnaTemperatura;
	double maxTezina;
	map<string, std::vector<int> > S;
public:
	std::vector<int> nasumicnoRjesenje(std::vector<int> zadatak);
	std::vector<int> nasumicniZadatak();
	std::vector<int> nasumicniZadatakDuljine(int n);
	UndirectedGraph inicijalizirajGraf(int brojOsoba, int brojVjestina);
	std::vector<int> izracunajNovoRjesenje(std::vector<int> trenutnoRjesenje, std::vector<int> zadatak);
	double izracunajEnergiju(std::vector<int> rjesenje);
	bool vjerojatnostPrijelaza(double deltaEnergije, double temp);
	double izracunajTemperaturu(double temp, double paramaterHladenja);
	UndirectedGraph prosiriGraf(UndirectedGraph graf, vector<int>* zadatak);
	double funkcijaDobrote(std::vector<int> rjesenje, double alfa);
	std::pair<int, std::vector<int> > nadjiRjesenjeUSusjedstvu(int brojPcela, int velSusjedstva, std::vector<int> zadatak, std::pair<int, std::vector<int> > parTrenutnoRjesenje, double alfa);
	double zbrojTezina(std::vector<int> rjesenje);
	void napraviGraf();
};

