#include "Graf.h"


template<typename UndirectedGraph>
Graf<UndirectedGraph>::Graf()
{
	postojiLiProsireni = false;
}

template<typename UndirectedGraph>
Graf<UndirectedGraph>::~Graf()
{
}

template<typename UndirectedGraph>
std::vector<int> Graf<UndirectedGraph>::simuliranoKaljenje(std::vector<int> zadatak, double temp, int maxIteracija)
{
	double parametarHladenja = 0.98;
	std::vector<int> trenutnoRjesenje = this->nasumicnoRjesenje(zadatak);
	int iteracija = 0;
	double energija = izracunajEnergiju(trenutnoRjesenje);
	while (temp > 0 && iteracija < maxIteracija)
	{
		std::vector<int> novoRjesenje = this->izracunajNovoRjesenje(trenutnoRjesenje, zadatak);
		double novaEnergija = izracunajEnergiju(novoRjesenje);
		double deltaEnergije = energija - novaEnergija;
		if (deltaEnergije > 0 || vjerojatnostPrijelaza(deltaEnergije, temp))
		{
			trenutnoRjesenje = novoRjesenje;
			energija = novaEnergija;
		}
		temp = izracunajTemperaturu(temp, parametarHladenja);		
		iteracija++;
	}
	return trenutnoRjesenje;
}

template<typename UndirectedGraph>
std::vector<int> Graf<UndirectedGraph>::dijametarAlgoritam(std::vector<int> zadatak)
{
	typedef graph_traits < UndirectedGraph >::vertex_descriptor vertex_descriptor;
	typedef graph_traits < UndirectedGraph >::edge_descriptor edge_descriptor;
	//S(v) izracunato u inicijalizaciji grafa
	int minVelicina = 0; 
	int indexNajrjedeVjestine = zadatak[0];
	minVelicina = S[vecVjestina[zadatak[0]]].size();
	for (int i = 1; i < zadatak.size(); ++i)
	{
		string trenVj = vecVjestina[zadatak[i]];
		int trenindexzad = zadatak[i];
		int trenutnaVelicina = S[vecVjestina[zadatak[i]]].size();
		if (trenutnaVelicina < minVelicina)
		{
			minVelicina = trenutnaVelicina;
			indexNajrjedeVjestine = zadatak[i];
		}
	}
	string v_rijedak = vecVjestina[indexNajrjedeVjestine];

	int globalnoMinDijametar = 0;
	int globalnoMinOsoba = 0;
	vector<vector<vertex_descriptor> > putovi;
	vector<int> najblizi;
	int gdjejenjabolji = 0;

	for (int i = 0; i < S[v_rijedak].size(); ++i)
	{
		std::vector<vertex_descriptor> p(num_vertices(G));
		std::vector<int> d(num_vertices(G));
		vertex_descriptor s = vertex(S[v_rijedak][i], G);

		property_map<UndirectedGraph, vertex_index_t>::type indexmap = get(vertex_index, G);
		property_map<UndirectedGraph, edge_weight_t>::type weightmap = get(edge_weight, G);
		dijkstra_shortest_paths(G, s, predecessor_map(&p[0])
			.distance_map(&d[0])
			.weight_map(boost::make_constant_property<edge_descriptor>(1ul)));

		int lokalnoMinDijametar = 100;
		int lokalnoMinOsoba = S[v_rijedak][0];

		for (int j = 0; j < zadatak.size(); ++j)
		{
			int lokalnoMaxDijametar = d[S[vecVjestina[zadatak[j]]][0]];
			int lokalnoMaxOsoba = S[vecVjestina[zadatak[j]]][0];
			string vjestina = vecVjestina[zadatak[j]];
			//racuna najmanji dijametar od osobe i. u S[v_rijedak] do neke osobe u S[j]
			for (int k = 0; k < S[vjestina].size(); ++k)
			{
				int osoba = S[vjestina][k];
				
				if (osoba!=S[v_rijedak][i] && d[osoba] < lokalnoMaxDijametar)
				{
					lokalnoMaxDijametar = d[osoba];
					lokalnoMaxOsoba = osoba;
				}
			}
			najblizi.push_back(lokalnoMaxOsoba);
			//naden lokalnoMaxDijametar

			if (j == 0)
			{
				lokalnoMinOsoba = lokalnoMaxOsoba;
				lokalnoMinDijametar = lokalnoMaxDijametar;
			}
			else
			{
				if (lokalnoMaxDijametar < lokalnoMinDijametar)
				{
					lokalnoMinDijametar = lokalnoMaxDijametar;
					lokalnoMinOsoba = lokalnoMaxOsoba;
				}
			}
		}
		//naden loklnoMinDijametar
		if (i == 0)
		{
			globalnoMinOsoba = S[v_rijedak][0];
			globalnoMinDijametar = lokalnoMinDijametar;
		}
		else
		{
			if (lokalnoMinDijametar < globalnoMinDijametar)
			{
				gdjejenjabolji = i;
				globalnoMinDijametar = lokalnoMinDijametar;
				globalnoMinOsoba = S[v_rijedak][i];
			}
		}
	}

	int i_zvjezdica = globalnoMinOsoba;
	std::set<int> skup;
	skup.insert(i_zvjezdica);
	vertex_descriptor zadnji = vertex(i_zvjezdica, G);

	std::vector<vertex_descriptor> p(num_vertices(G));
	std::vector<int> d(num_vertices(G));

	property_map<UndirectedGraph, vertex_index_t>::type indexmap = get(vertex_index, G);
	property_map<UndirectedGraph, edge_weight_t>::type weightmap = get(edge_weight, G);
	dijkstra_shortest_paths(G, zadnji, predecessor_map(&p[0])
		.distance_map(&d[0])
		.weight_map(boost::make_constant_property<edge_descriptor>(1ul)));

	for (int i = 0; i < zadatak.size(); ++i)
	{
		if (i != indexNajrjedeVjestine)
		{
			int najblizaOsoba = najblizi[i+gdjejenjabolji*zadatak.size()];
			vertex_descriptor trenutni = vertex(najblizaOsoba, G);
			while (trenutni != zadnji)
			{
				skup.insert(trenutni);
				trenutni = p[trenutni];
			}
		}
	}
	vector<int> rjesenje;
	rjesenje.assign(skup.begin(), skup.end());
	return rjesenje;
}

template<typename UndirectedGraph>
std::vector<int> Graf<UndirectedGraph>::SteinerStabloAlgoritam(UndirectedGraph* graf, std::vector<int> trazeniVrhovi)
{
	typedef graph_traits < UndirectedGraph >::vertex_descriptor vertex_descriptor;
	typedef graph_traits < UndirectedGraph >::edge_descriptor edge_descriptor;
	std::mt19937 rng;
	rng.seed(std::random_device()());
	std::uniform_int_distribution<std::mt19937::result_type> dist(0, trazeniVrhovi.size() - 1);
	int v = trazeniVrhovi[dist(rng)];
	std::set<int> rjesenje;
	rjesenje.insert(v);
	std::set<int> skupTrazenihVrhova;
	for (unsigned int i = 0; i < trazeniVrhovi.size(); ++i)
		skupTrazenihVrhova.insert(trazeniVrhovi[i]);

	while (!std::includes(rjesenje.begin(), rjesenje.end(), skupTrazenihVrhova.begin(), skupTrazenihVrhova.end()))
	{
		double minUdaljenostGlobalna = 10000;
		int indexMinOsobeGlobalno = 0;
		int v_zvjezdica = 0;
		std::vector<vertex_descriptor> p_najbolji;
		for(auto i =skupTrazenihVrhova.begin(); i != skupTrazenihVrhova.end(); ++i)
			if (std::find(rjesenje.begin(), rjesenje.end(), *i) == rjesenje.end())
			{
				std::vector<vertex_descriptor> p(num_vertices(*graf));
				std::vector<double> d(num_vertices(*graf));
				vertex_descriptor s = vertex(*i, *graf);

				property_map<UndirectedGraph, vertex_index_t>::type indexmap = get(vertex_index, *graf);
				property_map<UndirectedGraph, edge_weight_t>::type weightmap = get(edge_weight, *graf);
				dijkstra_shortest_paths(*graf, s, predecessor_map(&p[0])
					.distance_map(&d[0]));
				double minUdaljenost = 10000;
				int indexMinOsobe = 0;
				for (auto j = rjesenje.begin(); j != rjesenje.end(); ++j)
				{
					if (d[*j] < minUdaljenost)
					{
						minUdaljenost = d[*j];
						indexMinOsobe = *j;
					}
				}
				if (minUdaljenost < minUdaljenostGlobalna)
				{
					minUdaljenostGlobalna = minUdaljenost;
					indexMinOsobeGlobalno = indexMinOsobe;
					p_najbolji = p;
					v_zvjezdica = *i;
				}
			}
		vertex_descriptor trenutni = vertex(indexMinOsobeGlobalno, *graf);
		vertex_descriptor zadnji = vertex(v_zvjezdica, *graf);
		rjesenje.insert(zadnji);
		while (trenutni != zadnji)
		{
			rjesenje.insert(trenutni);
			if (trenutni == p_najbolji[trenutni])
				return std::vector<int>();
			trenutni = p_najbolji[trenutni];
		}
	}
	std::vector<int> konacnoRjesenje;
	konacnoRjesenje.assign(rjesenje.begin(), rjesenje.end());
	return konacnoRjesenje;
}

template<typename UndirectedGraph>
std::vector<int> Graf<UndirectedGraph>::pokrivacSteinerAlgoritam(std::vector<int> zadatak)
{
	std::vector<int> pokrivacZadatka = nasumicnoRjesenje(zadatak);
	std::vector<int> rjesenje = SteinerStabloAlgoritam(&G, pokrivacZadatka);
	return rjesenje;
}

template<typename UndirectedGraph>
std::vector<int> Graf<UndirectedGraph>::poboljsaniSteinerAlgoritam(UndirectedGraph* graf, std::vector<int> zadatak)
{
	prosireniGraf = prosiriGraf(*graf, &zadatak);
	std::vector<int> rjesenje =  SteinerStabloAlgoritam(&prosireniGraf, zadatak);
	rjesenje.resize(rjesenje.size()-zadatak.size());
	return rjesenje;
}

template<typename UndirectedGraph>
std::vector<int> Graf<UndirectedGraph>::pcelinjiAlgoritam(std::vector<int> zadatak, int maxIter, int n, int m, int e, int nep, int nsp, int ngh, double alfa)
{
	typedef std::vector<int> rjesenje;
	int iteracija = 0;
	std::vector<pair<int, rjesenje > > populcijaRjesenja;
	for (int i = 0; i < n; ++i)
	{
		rjesenje temp = nasumicnoRjesenje(zadatak);
		populcijaRjesenja.push_back(std::pair<int, rjesenje >(funkcijaDobrote(temp, alfa), temp));
	}
	std::sort(populcijaRjesenja.begin(), populcijaRjesenja.end(), mySort);

	while (iteracija < maxIter)
	{
		iteracija++;
		for (int i = 0; i < n; ++i)
		{
			if (i < e)
			{
				//obradi elitna rjesenja
				populcijaRjesenja[i] = nadjiRjesenjeUSusjedstvu(nep, ngh, zadatak, populcijaRjesenja[i], alfa);
			}
			else if (i < m)
			{
				//obradi najbolja, neelitna rjesenja
				populcijaRjesenja[i] = nadjiRjesenjeUSusjedstvu(nsp, ngh, zadatak, populcijaRjesenja[i], alfa);
			}
			else
			{
				//obradi ostala rjesenja
				rjesenje temp = nasumicnoRjesenje(zadatak);
				populcijaRjesenja[i] = std::pair<int, rjesenje >(funkcijaDobrote(temp, alfa), temp);
			}
		}
		//sortiraj trenutnu populaciju rjesenja
		std::sort(populcijaRjesenja.begin(), populcijaRjesenja.end(), mySort);
	}
	return populcijaRjesenja[0].second;
}


template<typename UndirectedGraph>
std::vector<int> Graf<UndirectedGraph>::nasumicnoRjesenje(std::vector<int> zadatak)
{
	std::vector<int> rjesenje;
	std::mt19937 rng;
	rng.seed(std::random_device()());
	for (unsigned int j = 0; j < zadatak.size(); j++) 
	{
		string vjestina = vecVjestina[zadatak[j]];
		std::uniform_int_distribution<std::mt19937::result_type> dist(0, S[vjestina].size()-1);
		rjesenje.push_back(S[vjestina][dist(rng)]);
	}
	return rjesenje;
}

template<typename UndirectedGraph>
std::vector<int> Graf<UndirectedGraph>::nasumicniZadatak()
{
	std::vector<int> rjesenje;
	std::mt19937 rng;
	rng.seed(std::random_device()());
	std::uniform_int_distribution<std::mt19937::result_type> dist6(1, 100);
	for (int j = 0; j < vecVjestina.size(); j++) {
		if (S[vecVjestina[j]].size() != 0 && dist6(rng)>50)
			rjesenje.push_back(j);
	}
	return rjesenje;
}

template<typename UndirectedGraph>
std::vector<int> Graf<UndirectedGraph>::nasumicniZadatakDuljine(int n)
{
	std::set<int> zadatakSkup;
	std::mt19937 rng;
	rng.seed(std::random_device()());
	std::uniform_int_distribution<std::mt19937::result_type> dist6(0, vecVjestina.size()-1);
	while (zadatakSkup.size() != n)
	{
		int indexVjestine = dist6(rng);
		string vjestina = vecVjestina[indexVjestine];
		if(S[vjestina].size()!=0)
			zadatakSkup.insert(indexVjestine);
	}
	std::vector<int> zadatak;
	zadatak.assign(zadatakSkup.begin(), zadatakSkup.end());
	return zadatak;
}

template<typename UndirectedGraph>
UndirectedGraph Graf<UndirectedGraph>::inicijalizirajGraf(int brojOsoba, int brojVjestina)
{
	std::vector<string> imena{ "Ana", "Boris", "Cvita", "Dalija", "Ema", "Filip", "Gordana", "Helena", "Ivana", "Jakov", "Konstantin", "Luka",
		"Mia", "Nika", "Petra", "Roko", "Stefan", "Teodora", "Una", "Vanja", "Zlatan" };
	std::vector<string> vjestine{ "C#", "C++", "C", "HTML", "Python", "PHP", "JavaScript", "AI", "AR", "BigData", "Android", "iOS", "Cyber security", "ML",
		"Analytics", "Teamwork", "API", "Design", "Soft skills", "Linux", "Servers", "UX/UI", "Tensorflow" };
	int vel = vjestine.size();
	for (int i = 0; i < brojVjestina; ++i)
		vecVjestina.push_back(vjestine[i % vel] + to_string(i / vel));
	vecImena = imena;

	UndirectedGraph graf(brojOsoba);
	
	int n = brojOsoba;
	int broj_bridova = 0;

	typename boost::graph_traits<UndirectedGraph>::out_edge_iterator e, e_end;
	typedef typename UndirectedGraph::edge_property_type Weight;
	typename property_map < UndirectedGraph, edge_weight_t >::type
		weight = get(edge_weight, graf);

	if (true)
	{
		for (int i = 0; i < n; i++) {
			typename graph_traits < UndirectedGraph >::vertex_descriptor u, v;
			u = vertex(i, graf);
			graf[u].name = vecImena[i%vecImena.size()] + to_string((i/vecImena.size()));
			std::mt19937 rng;
			rng.seed(std::random_device()());
			std::uniform_int_distribution<std::mt19937::result_type> dist6(1, 100);
			for (int j = 0; j < vecVjestina.size(); j++) {
				if (dist6(rng)>95)
				{
					graf[u].skills.push_back(vecVjestina[j]);
					this->S[vecVjestina[j]].push_back(i);
				}
			}
			for (int k = 0; k < n; k++) {
				if (k != i && dist6(rng)>95) {
					v = vertex(k, graf);
					double trenutnaTezina = dist6(rng) / 100.0;
					if(!boost::edge(u,v,graf).second) 
						add_edge(u, v, Weight(trenutnaTezina), graf);
					if (trenutnaTezina > maxTezina)
						maxTezina = trenutnaTezina;
					broj_bridova++;
				}
			}
			
		}

	}
	G = graf;
	return graf;
}

template<typename UndirectedGraph>
std::vector<int> Graf<UndirectedGraph>::izracunajNovoRjesenje(std::vector<int> trenutnoRjesenje, std::vector<int> zadatak)
{
	std::mt19937 rng;
	rng.seed(std::random_device()());
	std::uniform_int_distribution<std::mt19937::result_type> dist6(0, zadatak.size()-1);

	int poredak = dist6(rng);
	int vjestina = zadatak[poredak];
	string imeVjestine = vecVjestina[vjestina];
	std::uniform_int_distribution<std::mt19937::result_type> dist(0, this->S[imeVjestine].size() - 1);
	trenutnoRjesenje[poredak] = S[imeVjestine][dist(rng)];

	return trenutnoRjesenje;
}

template<typename UndirectedGraph>
double Graf<UndirectedGraph>::izracunajEnergiju(std::vector<int> rjesenje)
{
	set<int> s;
	unsigned size = rjesenje.size();
	for (unsigned i = 0; i < size; ++i) s.insert(rjesenje[i]);
	rjesenje.assign(s.begin(), s.end());

	double tezina = 0;
	int brojBridova = 0;
	double beta = 4;
	for (unsigned int i = 0; i < rjesenje.size(); ++i)
	{
		typename graph_traits < UndirectedGraph >::vertex_descriptor u;
		u = vertex(rjesenje[i], G);
		typename boost::graph_traits<UndirectedGraph>::out_edge_iterator e, e_end;
		typedef typename UndirectedGraph::edge_property_type Weight;
		typename property_map < UndirectedGraph, edge_weight_t >::type
			weight = get(edge_weight, G);
		for (boost::tie(e, e_end) = out_edges(u, G); e != e_end; ++e)
		{
			int cvor = target(*e, G);
			if (std::find(rjesenje.begin(), rjesenje.end(), target(*e, G)) != rjesenje.end())
			{
				auto trenutnaTezina = get(weight, *e);
				tezina += trenutnaTezina;
				++brojBridova;
			}
		}
	}
	brojBridova /= 2;
	int brojNepostojecihBridova = rjesenje.size()*(rjesenje.size() - 1) / 2 - brojBridova;
	return tezina/2 + brojNepostojecihBridova * beta * maxTezina;
}

template<typename UndirectedGraph>
bool Graf<UndirectedGraph>::vjerojatnostPrijelaza(double deltaEnergije, double temp)
{
	std::mt19937 rng;
	rng.seed(std::random_device()());
	std::uniform_int_distribution<std::mt19937::result_type> dist(1, 100);
	if (exp(-deltaEnergije / temp) < 1.3)
	{
		return true;
	}
	else 
		return false;
}

template<typename UndirectedGraph>
double Graf<UndirectedGraph>::izracunajTemperaturu(double temp, double paramaterHladenja)
{
	return temp*paramaterHladenja;
}

template<typename UndirectedGraph>
UndirectedGraph Graf<UndirectedGraph>::prosiriGraf(UndirectedGraph graf, vector<int>* zadatak)
{
	typedef graph_traits < UndirectedGraph >::vertex_descriptor vertex_descriptor;
	typedef graph_traits < UndirectedGraph >::edge_descriptor edge_descriptor;
	typedef typename UndirectedGraph::edge_property_type Weight;
	vector<int> indexiZadatka;
	for (auto i = zadatak->begin(); i != zadatak->end(); ++i)
	{
		vertex_descriptor vjestina = add_vertex(graf);
		graf[vjestina].name = "vjestina";
		indexiZadatka.push_back(vjestina);
		for (auto j = S[vecVjestina[*i]].begin(); j != S[vecVjestina[*i]].end(); ++j)
		{
			vertex_descriptor osoba = vertex(*j, graf);
			add_edge(osoba, vjestina, Weight(1000), graf);
		}
	}
	*zadatak = indexiZadatka;
	return graf;
}

template<typename UndirectedGraph>
double Graf<UndirectedGraph>::funkcijaDobrote(std::vector<int> rjesenje, double alfa)
{
	set<int> s;
	unsigned size = rjesenje.size();
	for (unsigned i = 0; i < size; ++i) s.insert(rjesenje[i]);
	rjesenje.assign(s.begin(), s.end());

	double tezina = 0;
	int brojBridova = 0;
	for (unsigned int i = 0; i < rjesenje.size(); ++i)
	{
		typename graph_traits < UndirectedGraph >::vertex_descriptor u;
		u = vertex(rjesenje[i], G);
		typename boost::graph_traits<UndirectedGraph>::out_edge_iterator e, e_end;
		typedef typename UndirectedGraph::edge_property_type Weight;
		typename property_map < UndirectedGraph, edge_weight_t >::type
			weight = get(edge_weight, G);
		for (boost::tie(e, e_end) = out_edges(u, G); e != e_end; ++e)
		{
			int cvor = target(*e, G);
			if (std::find(rjesenje.begin(), rjesenje.end(), target(*e, G)) != rjesenje.end())
			{
				auto trenutnaTezina = 1.0/get(weight, *e);
				tezina += trenutnaTezina;
				++brojBridova;
			}
		}
	}
	brojBridova /= 2;
	int brojNepostojecihBridova = rjesenje.size()*(rjesenje.size() - 1) / 2 - brojBridova;
	return tezina / 2 - brojNepostojecihBridova * alfa * 1.0/maxTezina;
}


static bool mySort(std::pair<int, std::vector<int>> prvi, std::pair<int, std::vector<int>> drugi)
{
	return prvi.first>drugi.first;
}

template<typename UndirectedGraph>
std::pair<int, std::vector<int> > Graf<UndirectedGraph>::nadjiRjesenjeUSusjedstvu(int brojPcela, int velSusjedstva, std::vector<int> zadatak, std::pair<int, std::vector<int> > parTrenutnoRjesenje, double alfa)
{
	typedef std::vector<int> rjesenje;
	std::vector<pair<int, rjesenje > > populcijaRjesenja;
	populcijaRjesenja.push_back(parTrenutnoRjesenje);
	for (int i = 0; i < brojPcela; ++i)
	{
		rjesenje temp;
		for(int j=0; j<velSusjedstva; ++j)
			temp = izracunajNovoRjesenje(parTrenutnoRjesenje.second, zadatak);
		populcijaRjesenja.push_back(std::pair<int, rjesenje >(funkcijaDobrote(temp, alfa), temp));
	}
	std::sort(populcijaRjesenja.begin(), populcijaRjesenja.end(), mySort);

	return populcijaRjesenja[0];
}

template<typename UndirectedGraph>
double Graf<UndirectedGraph>::zbrojTezina(std::vector<int> rjesenje)
{
	set<int> s;
	unsigned size = rjesenje.size();
	for (unsigned i = 0; i < size; ++i) s.insert(rjesenje[i]);
	rjesenje.assign(s.begin(), s.end());

	double tezina = 0;
	int brojBridova = 0;
	for (unsigned int i = 0; i < rjesenje.size(); ++i)
	{
		typename graph_traits < UndirectedGraph >::vertex_descriptor u;
		u = vertex(rjesenje[i], G);
		typename boost::graph_traits<UndirectedGraph>::out_edge_iterator e, e_end;
		typedef typename UndirectedGraph::edge_property_type Weight;
		typename property_map < UndirectedGraph, edge_weight_t >::type
			weight = get(edge_weight, G);
		for (boost::tie(e, e_end) = out_edges(u, G); e != e_end; ++e)
		{
			int cvor = target(*e, G);
			if (std::find(rjesenje.begin(), rjesenje.end(), target(*e, G)) != rjesenje.end())
			{
				auto trenutnaTezina = get(weight, *e);
				tezina += trenutnaTezina;
				++brojBridova;
			}
		}
	}
	brojBridova /= 2;
	return tezina / 2;
}

template<typename UndirectedGraph>
void Graf<UndirectedGraph>::napraviGraf()
{
	std::map<int, std::map<int, int>> susjedi;
	std::map<int, std::vector<std::string>> vjestine1;

	std::ifstream ifs1("save_electronics");
	boost::archive::text_iarchive ia1(ifs1);
	// write class instance to archive
	ia1 >> vjestine1;

	std::ifstream ifs2("save_electronics_snapshot");
	boost::archive::text_iarchive ia2(ifs2);
	// write class instance to archive
	ia2 >> susjedi;

	std::map<int, int> indexi;
	std::set<string> skupVjestina;
	int j = 0;
	for (auto i = susjedi.begin(); i != susjedi.end(); ++i)
	{
		indexi[i->first] = j;
		++j;
		if (i->first == -1 || i->first == 0)
			continue;
		if (i->second.size() < 3)
			for (auto j = i->second.begin(); j != i->second.end(); ++j)
			{
				if(j->first!=i->first)
					susjedi[j->first].erase(i->first);
			}
	}

	int n = susjedi.size();

	UndirectedGraph graf(n);
	this->maxTezina = 0;

	typename boost::graph_traits<UndirectedGraph>::out_edge_iterator e, e_end;
	typedef typename UndirectedGraph::edge_property_type Weight;
	
	int brojBridova = 0;
	for (auto i = susjedi.begin(); i != susjedi.end(); ++i)
	{
		if (i->first == -1 || i->first == 0 || i->second.size()<3)
			continue;

		int tren_index = indexi[i->first];
		
		typename boost::graph_traits < UndirectedGraph >::vertex_descriptor u, v;
		u = vertex(tren_index, graf);
		
		for (auto k = vjestine1.at(i->first).begin(); k != vjestine1.at(i->first).end(); ++k)
		{
			graf[u].skills.push_back(*k);
			S[*k].push_back(tren_index);
			skupVjestina.insert(*k);
		}
		for (auto j = i->second.begin(); j != i->second.end(); ++j)
		{
			int temp_index = indexi[j->first];
			if (temp_index == -1 || temp_index == 0)
				continue;
			v = vertex(temp_index, graf);
			if (tren_index < temp_index)
			{
				double trenTezina = 1.0 / (j->second);
				if (trenTezina > maxTezina)
					maxTezina = trenTezina;
				add_edge(u, v, Weight(1.0 / (j->second)), graf);
				++brojBridova;
			}
		}
	}
	G = graf;
	vecVjestina.assign(skupVjestina.begin(), skupVjestina.end());
}
