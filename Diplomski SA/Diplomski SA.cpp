// Diplomski SA.cpp : Defines the entry point for the console application.
//
#include "tinyxml2.h"

#include "stdafx.h"

using namespace std;
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "Graf.cpp"
#include <boost/graph/graphviz.hpp>
#include <chrono>
#include <filesystem>

using namespace boost;
namespace fs = std::experimental::filesystem;


template < typename UndirectedGraph >
void ucitajPodatke(Graf<UndirectedGraph>* g)
{
	std::string path = "C:\\Users\\Branimirko\\Documents\\diplomski_zipovi_SE\\unzip_meta";
	//path = "C:\\Users\\Branimirko\\Documents\\diplomski_zipovi_SE\\test";
	std::map<int, std::map<int, int> > neighbors;
	std::map<int, std::vector<std::string> > S_vjestine;
	int brojac = 0;

	/*
	for (auto & p : fs::directory_iterator(path))
	{
		std::string vjestina = p.path().filename().string();
		if (vjestina.find("desktop") == 0)
			continue;
		vjestina = vjestina.substr(0, vjestina.find("."));

		std::cout << vjestina << std::endl;

		std::string temp = (p / "Posts.xml").string();
		const char * post_path = temp.c_str();
		std::string temp2 = (p / "Users.xml").string();
		const char * users_path = temp2.c_str();

		//std::cout << p / "Posts.xml";

		tinyxml2::XMLDocument posts_xml;
		tinyxml2::XMLDocument users_xml;

		tinyxml2::XMLError eResult = posts_xml.LoadFile(post_path);
		if (eResult != tinyxml2::XML_SUCCESS) break;
		tinyxml2::XMLError usersResult = users_xml.LoadFile(users_path);
		if (usersResult != tinyxml2::XML_SUCCESS) break;

		std::map<int, int> ownerIDs;
		std::map<int, int> accountIDs;
		tinyxml2::XMLNode* root_users = users_xml.FirstChildElement("users");
		if (root_users == nullptr) break;
		tinyxml2::XMLElement* users_element = root_users->FirstChildElement("row");
		if (users_element == nullptr) break;

		for (users_element; users_element != NULL; users_element = users_element->NextSiblingElement())
		{
			int temp = users_element->IntAttribute("AccountId");
			accountIDs[users_element->IntAttribute("Id")] = temp;
			if (temp == -1)
				continue;
			S_vjestine[temp].push_back(vjestina);
		}
		std::ofstream ofs("save_" + vjestina);

		// save data to archive

		boost::archive::text_oarchive oa(ofs);
		// write class instance to archive
		oa << S_vjestine;
		std::cout << "snapshot vjestina spremljen" << std::endl;


		tinyxml2::XMLNode* root = posts_xml.FirstChildElement("posts");
		if (root == nullptr) break;

		tinyxml2::XMLElement* element = root->FirstChildElement("row");
		if (element == nullptr) break;
		do {
			int ownerID = 0, postID = 0, parentID = 0;
			//std::cout << "(ownerId, ParentId (ako ima), Id): ";
			if (element->Attribute("OwnerUserId") != NULL)
			{
				element->QueryIntAttribute("OwnerUserId", &ownerID);
				ownerIDs[element->IntAttribute("Id")] = ownerID;
			}
			if (element->Attribute("ParentId") != NULL)
			{
				element->QueryIntAttribute("ParentId", &parentID);
				parentID = ownerIDs[parentID];
			}
			if (element->Attribute("Id") != NULL)
			{
				element->QueryIntAttribute("Id", &postID);
			}
			ownerID = accountIDs[ownerID];
			parentID = accountIDs[parentID];
			//std::cout << "(" << ownerID;
			//std::cout << ", " << parentID;
			//std::cout << ", " << postID << ")" << std::endl;

			if (parentID != 0)
			{
				neighbors[ownerID][parentID]++;
				neighbors[parentID][ownerID]++;
			}
			element = element->NextSiblingElement();

		} while (element != NULL);

		std::ofstream ofs3("save_" + vjestina + "_snapshot");

		// save data to archive

		boost::archive::text_oarchive oa2(ofs3);
		// write class instance to archive
		oa2 << neighbors;
		std::cout << "snapshot susjeda spremljen" << std::endl;
	}

	std::ofstream ofs2("susjedi_final");

	// save data to archive

	boost::archive::text_oarchive oa2(ofs2);
	// write class instance to archive
	oa2 << neighbors;

	std::cout << "svi susjedi spremljen" << std::endl;

	*/
	/*
	for (auto i = neighbors.begin(); i != neighbors.end(); ++i)
	{
	std::cout << "susjedi od " << i->first << "["<< i->second.size() <<"]: ";
	for (auto j = i->second.begin(); j != i->second.end(); ++j)
	std::cout << j->first << "(" << j->second << ")" <<", ";
	std::cout << std::endl;
	}
	*/

	g->napraviGraf();
	return;
}


int main()
{
	typedef property < edge_weight_t, double >Weight;
	struct Person
	{
		string name;
		vector<string> skills;
	};
	typedef adjacency_list < vecS, vecS, undirectedS,
		Person, Weight > UndirectedGraph;

	//undirected_graph_demo2 < UndirectedGraph >();
	Graf<UndirectedGraph> g;
	//ucitajPodatke<UndirectedGraph>(&g);

	UndirectedGraph grafic = g.inicijalizirajGraf(500,60);

	//g.napraviGraf(&neighbors, &vjestine1, &vjestine2);

	typedef boost::graph_traits<UndirectedGraph>::vertex_iterator v_iter;
	std::pair<v_iter, v_iter> vp;

	typedef property_map<UndirectedGraph, vertex_index_t>::type IndexMap;
	//IndexMap index = get(vertex_index, grafic);

	/*
	for (vp = vertices(grafic); vp.first != vp.second; ++vp.first)
	{
		std::cout << "susjedi od " << grafic[*vp.first].name << " su: ";
		typename boost::graph_traits<UndirectedGraph>::out_edge_iterator e, e_end;

		typename boost::graph_traits<UndirectedGraph>::vertex_descriptor
			s = vertex(index[*vp.first], grafic);
		for (boost::tie(e, e_end) = out_edges(s, grafic); e != e_end; ++e)
			std::cout << grafic[target(*e, grafic)].name << " ";
		std::cout << endl;

		std::cout << "vjestine: ";
		for_each(grafic[s].skills.begin(), grafic[s].skills.end(), [](string name) {std::cout << name << ", "; });
		std::cout << std::endl;
	}
	*/
	/*
	std::vector<int> randomZad = g.nasumicniZadatakDuljine(8);

	//write_graphviz(std::cout, grafic);

	std::cout << "nasumicni zadatak: ";
	for (int i = 0; i < randomZad.size(); ++i)
		std::cout << randomZad[i] << ", ";
	std::cout << endl;

	std::cout << "OVDJE::: " << endl << endl << endl;

	for (int i = 0; i < g.S["C++"].size(); ++i)
		std::cout << g.S["C++"][i] << "" << endl;

	for (int i = 0; i < randomZad.size(); ++i)
	{
		string imeVj = g.vecVjestina[randomZad[i]];
		std::cout << randomZad[i] << "/" << imeVj << "[" << g.S[imeVj].size() <<"]" <<": ";
		for (int j = 0; j < g.S[imeVj].size(); ++j)
		{
			int vel = g.vecImena.size();
			std::cout << g.S[imeVj][j] << "/" << g.vecImena[g.S[imeVj][j] % vel] + to_string((g.S[imeVj][j] / vel)) << ", ";
		}
		std::cout << endl;

	}
	*/
	

	/*
	std::cout << "nasumicno rjesenje: ";

	randomRj = g.nasumicnoRjesenje(randomZad);
	for (int i = 0; i < randomRj.size(); ++i)
		std::cout << randomRj[i] << "/" << g.vecImena[randomRj[i] % g.vecImena.size()] << to_string(randomRj[i] / g.vecImena.size()) << ", ";
	std::cout << endl;

	std::cout << "energija: " << g.izracunajEnergiju(randomRj) << endl;
	
	std::cout << endl;
	*/
	/*
	std::vector<int> randomRj;

	
	std::cout << endl << endl << endl << " DIJAMETAR ALGORITAM :" << endl;

	randomRj = g.dijametarAlgoritam(randomZad);
	for (int i = 0; i < randomRj.size(); ++i)
		std::cout << randomRj[i] << "/" << g.vecImena[randomRj[i] % g.vecImena.size()] << to_string(randomRj[i] / g.vecImena.size()) << ", ";
	std::cout << endl << "energija: " << g.izracunajEnergiju(randomRj) << endl;
	

	std::cout << endl << endl << endl << " POKRIVAC STEINER ALGORITAM :" << endl;


	randomRj = g.pokrivacSteinerAlgoritam(g.G, randomZad);
	for (int i = 0; i < randomRj.size(); ++i)
		std::cout << randomRj[i] << "/" << g.vecImena[randomRj[i] % g.vecImena.size()] << to_string(randomRj[i] / g.vecImena.size()) << ", ";
	std::cout << endl << "energija: " << g.izracunajEnergiju(randomRj) << endl;

	std::cout << endl << endl << endl << " POBOLJSANI STEINER ALGORITAM :" << endl;


	randomRj = g.poboljsaniSteinerAlgoritam(g.G, randomZad);
	for (int i = 0; i < randomRj.size(); ++i)
		std::cout << randomRj[i] << "/" << g.vecImena[randomRj[i] % g.vecImena.size()] << to_string(randomRj[i] / g.vecImena.size()) << ", ";
	std::cout << endl << "energija: " << g.izracunajEnergiju(randomRj) << endl;

	
	std::cout << endl << endl << endl << "SIMULIRANO KALJENJE :" << endl;


	randomRj = g.simuliranoKaljenje(randomZad, 100, 200);
	for (int i = 0; i < randomRj.size(); ++i)
		std::cout << randomRj[i] << "/" << g.vecImena[randomRj[i] % g.vecImena.size()] << to_string(randomRj[i] / g.vecImena.size()) << ", ";
	std::cout << endl << "energija: " << g.izracunajEnergiju(randomRj) << endl;




	std::cout << endl << endl << endl << " PCELINJI ALGORITAM 2:" << endl;


	randomRj = g.pcelinjiAlgoritam(randomZad, 100, 5, 3, 1, 3, 2, 1, 2);
	for (int i = 0; i < randomRj.size(); ++i)
		std::cout << randomRj[i] << "/" << g.vecImena[randomRj[i] % g.vecImena.size()] << to_string(randomRj[i] / g.vecImena.size()) << ", ";
	std::cout << endl << "energija: " << g.izracunajEnergiju(randomRj) << "; funkcija dobrote: " << g.funkcijaDobrote(randomRj, 2) << endl;

	std::cout << endl << endl << endl << " PCELINJI ALGORITAM 3:" << endl;


	randomRj = g.pcelinjiAlgoritam(randomZad, 100, 5, 3, 1, 3, 2, 1, 3);
	for (int i = 0; i < randomRj.size(); ++i)
		std::cout << randomRj[i] << "/" << g.vecImena[randomRj[i] % g.vecImena.size()] << to_string(randomRj[i] / g.vecImena.size()) << ", ";
	std::cout << endl << "energija: " << g.izracunajEnergiju(randomRj) <<"; funkcija dobrote: " << g.funkcijaDobrote(randomRj, 3) << endl;

	std::cout << endl << endl << endl << " PCELINJI ALGORITAM 3.5:" << endl;


	randomRj = g.pcelinjiAlgoritam(randomZad, 100, 5, 3, 1, 3, 2, 1, 3.5);
	for (int i = 0; i < randomRj.size(); ++i)
		std::cout << randomRj[i] << "/" << g.vecImena[randomRj[i] % g.vecImena.size()] << to_string(randomRj[i] / g.vecImena.size()) << ", ";
	std::cout << endl << "energija: " << g.izracunajEnergiju(randomRj) << "; funkcija dobrote: " << g.funkcijaDobrote(randomRj, 3.5) << endl;

	std::cout << endl << endl << endl << " PCELINJI ALGORITAM 4:" << endl;


	randomRj = g.pcelinjiAlgoritam(randomZad, 100, 5, 3, 1, 3, 2, 1, 4);
	for (int i = 0; i < randomRj.size(); ++i)
		std::cout << randomRj[i] << "/" << g.vecImena[randomRj[i] % g.vecImena.size()] << to_string(randomRj[i] / g.vecImena.size()) << ", ";
	std::cout << endl << "energija: " << g.izracunajEnergiju(randomRj) << "; funkcija dobrote: " << g.funkcijaDobrote(randomRj, 4) << endl;

	std::cout << endl << endl << endl << " PCELINJI ALGORITAM 4.5:" << endl;


	randomRj = g.pcelinjiAlgoritam(randomZad, 100, 5, 3, 1, 3, 2, 1, 4);
	for (int i = 0; i < randomRj.size(); ++i)
		std::cout << randomRj[i] << "/" << g.vecImena[randomRj[i] % g.vecImena.size()] << to_string(randomRj[i] / g.vecImena.size()) << ", ";
	std::cout << endl << "energija: " << g.izracunajEnergiju(randomRj) << "; funkcija dobrote: " << g.funkcijaDobrote(randomRj, 4.5) << endl;
	*/

	std::cout << "iteracija, velZadatka, DA e, DA d, DA zbrT, DA time , PSA e, PSA d, PSA zbrT, PSA time,";
	std::cout << "PobSA e, PobSA d, PobSA zbrT, PosB time, SK e, SK d, SK zbrT, SK time,";
	std::cout << "PA e, PA d, PA zbrT, PA time\n";
	for (int i = 0; i < 100; ++i)
	{
		int velicinaZadatka = (i + 10) / 10 + 5;
		std::vector<int> zadatak = g.nasumicniZadatakDuljine(velicinaZadatka);
		int maxIteracija = 1000;
		std::vector<double> energija(5);
		std::vector<double> dobrota(5);
		std::vector<double> zbrojTezina(5);

		std::vector<int> dobivenoRjesenje;
		std::cout << i << "," << velicinaZadatka;
		for (int j = 0; j < 5; ++j)
		{
			auto start = std::chrono::high_resolution_clock::now();
			switch (j) {
			case 0: dobivenoRjesenje = g.dijametarAlgoritam(zadatak); break;
			case 1: dobivenoRjesenje = g.pokrivacSteinerAlgoritam(zadatak); break;
			case 2: dobivenoRjesenje = g.poboljsaniSteinerAlgoritam(&g.G, zadatak); break;
			case 3: dobivenoRjesenje = g.simuliranoKaljenje(zadatak, 1000, maxIteracija); break; 
			case 4: dobivenoRjesenje = g.pcelinjiAlgoritam(zadatak, maxIteracija, 10, 5, 2, 5, 3, 2, 3);
			}
			auto finish = std::chrono::high_resolution_clock::now();

			auto time = (finish - start).count();

			energija[j] = g.izracunajEnergiju(dobivenoRjesenje);
			dobrota[j] = g.funkcijaDobrote(dobivenoRjesenje, 3);
			zbrojTezina[j] = g.zbrojTezina(dobivenoRjesenje);

			std::cout << "," << energija[j] << "," << dobrota[j] << "," << zbrojTezina[j] << "," << time;

		}
		std::cout << std::endl;
		 
	}

    return 0;
}

