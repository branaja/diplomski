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
	//funkcija koja od arhive objava na StackExchange forumu sprema graf osoba i bridova među njima
	std::string path = "C:\\Users\\Branimirko\\Documents\\diplomski_zipovi_SE\\unzip_meta";
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

	
	Graf<UndirectedGraph> g;

	//odkomentirati za slučajan graf
	//UndirectedGraph grafic = g.inicijalizirajGraf(500,60);
	
	//odkomentirati za StackExchange graf
	//g.napraviGraf();

	typedef boost::graph_traits<UndirectedGraph>::vertex_iterator v_iter;
	std::pair<v_iter, v_iter> vp;

	typedef property_map<UndirectedGraph, vertex_index_t>::type IndexMap;

	//testiranje 5 algoritama na danom grafu. 100 iteracija, po 10 za istu veličinu zadatka
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

