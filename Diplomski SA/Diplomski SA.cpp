// Diplomski SA.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

using namespace std;
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "Graf.cpp"
#include <boost/graph/graphviz.hpp>
using namespace boost;

template < typename UndirectedGraph > void
undirected_graph_demo2()
{
	const int V = 4;
	UndirectedGraph undigraph(V);
	typename graph_traits < UndirectedGraph >::vertex_descriptor u, v, x, y;
	typedef typename UndirectedGraph::edge_property_type Weight;
	typename property_map < UndirectedGraph, edge_weight_t >::type
		weight = get(edge_weight, undigraph);
	typename graph_traits < UndirectedGraph >::edge_descriptor e1, e2;
	bool found;

	u = vertex(0, undigraph);
	v = vertex(1, undigraph);
	x = vertex(2, undigraph);
	y = vertex(3, undigraph);
	undigraph[u].name = "Ana";
	undigraph[u].skills.push_back("c++");
	undigraph[u].skills.push_back("c");
	undigraph[u].skills.push_back("python");
	undigraph[v].name = "Bero";
	undigraph[x].name = "Cico";
	undigraph[y].name = "Darko";
	add_edge(u, v, Weight(3.1), undigraph);
	add_edge(x, y, Weight(3.4), undigraph);
	add_edge(u, x, Weight(2.1), undigraph);
	add_edge(x, v, Weight(4.5), undigraph);
	boost::tie(e1, found) = edge(u, v, undigraph);
	boost::tie(e2, found) = edge(v, u, undigraph);
	std::cout << "in an undirected graph is ";
#ifdef __GNUC__
	std::cout << "(u,v) == (v,u) ? " << (e1 == e2) << std::endl;
#else
	std::cout << "(u,v) == (v,u) ? "
		<< std::boolalpha << (e1 == e2) << std::endl;
#endif
	std::cout << "weight[(u,v)] = " << get(weight, e1) << std::endl;
	std::cout << "weight[(v,u)] = " << get(weight, e2) << std::endl;

	std::cout << "the edges incident to 2: ";
	typename boost::graph_traits<UndirectedGraph>::out_edge_iterator e, e_end;
	typename boost::graph_traits<UndirectedGraph>::vertex_descriptor
		s = vertex(2, undigraph);
	for (boost::tie(e, e_end) = out_edges(s, undigraph); e != e_end; ++e)
		std::cout << "(" << source(*e, undigraph)
		<< "," << target(*e, undigraph) << ")" << " ";
	std::cout << endl;

	typedef boost::graph_traits<UndirectedGraph>::vertex_iterator v_iter;
	std::pair<v_iter, v_iter> vp;

	typedef property_map<UndirectedGraph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, undigraph);

	for (vp = vertices(undigraph); vp.first != vp.second; ++vp.first)
	{
		std::cout << "susjedi od " << undigraph[*vp.first].name << " su: ";
		typename boost::graph_traits<UndirectedGraph>::out_edge_iterator e, e_end;

		typename boost::graph_traits<UndirectedGraph>::vertex_descriptor
			s = vertex(index[*vp.first], undigraph);
		for (boost::tie(e, e_end) = out_edges(s, undigraph); e != e_end; ++e)
			std::cout << undigraph[target(*e, undigraph)].name << " ";
		std::cout << endl;

		std::cout << "vjestine: ";
		for_each(undigraph[s].skills.begin(), undigraph[s].skills.end(), [](string name) {std::cout << name << ", "; });
		std::cout << std::endl;
	}
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

	UndirectedGraph grafic = g.inicijalizirajGraf(400,50);

	typedef boost::graph_traits<UndirectedGraph>::vertex_iterator v_iter;
	std::pair<v_iter, v_iter> vp;

	typedef property_map<UndirectedGraph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, grafic);

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
	std::vector<int> randomZad = g.nasumicniZadatak();

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

	

	/*
	std::cout << "nasumicno rjesenje: ";

	randomRj = g.nasumicnoRjesenje(randomZad);
	for (int i = 0; i < randomRj.size(); ++i)
		std::cout << randomRj[i] << "/" << g.vecImena[randomRj[i] % g.vecImena.size()] << to_string(randomRj[i] / g.vecImena.size()) << ", ";
	std::cout << endl;

	std::cout << "energija: " << g.izracunajEnergiju(randomRj) << endl;
	
	std::cout << endl;
	*/
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

    return 0;
}

