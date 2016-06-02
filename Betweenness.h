#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <vector>
#include <stack>
#include <queue>
#include <boost/property_map/property_map.hpp>
//#include <boost/test/minimal.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/lexical_cast.hpp>
#include "Element.h"


//using namespace boost;

#ifndef BETWEENNESS_H
#define BETWEENNESS_H

template<int N>
class CalcBetweenness
{
	Element<N + 1>* KNN;
	XParams* para;

	int m_index[N];

	
	typedef boost::property<boost::edge_weight_t, double, boost::property<boost::edge_index_t, std::size_t> > EdgeProperties;
	typedef boost::adjacency_list<boost::listS, boost::listS, boost::directedS,//undirectedS,
		boost::property<boost::vertex_index_t, int>, EdgeProperties> Graph;

	typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
	typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
	
	string m_filepath, m_filename;

public:


	void test()
	{
		typedef typename boost::graph_traits<Graph>::adjacency_iterator adjacency_iterator;
		Graph g(7);
		std::vector<Vertex> vertices(7);
		std::vector<vertex_iterator> ver_itr(7);
		{
			vertex_iterator v, v_end;
			int index = 0;
			for (boost::tie(v, v_end) = boost::vertices(g); v != v_end; ++v, ++index) 
			{
				put(boost::vertex_index, g, *v, index);
				vertices[index] = *v;
				ver_itr[index] = v;
			}
		}
		boost::add_edge(vertices[0], vertices[1], g);
		boost::add_edge(vertices[0], vertices[2], g);
		boost::add_edge(vertices[0], vertices[3], g);
		boost::add_edge(vertices[4], vertices[5], g);
		boost::add_edge(vertices[3], vertices[4], g);
		boost::add_edge(vertices[2], vertices[4], g);
		
		std::cout << "graph with " << 7 << " vertices and " << num_edges(g) << " edges.\n";


		std::cout << "  Direct translation of Brandes' algorithm...";
		std::vector<double> centrality(7);
		simple_unweighted_betweenness_centrality(g, get(boost::vertex_index, g), make_iterator_property_map(centrality.begin(), get(boost::vertex_index, g), double()));
		std::cout << "DONE.\n";


		vertex_iterator si, si_end;

		cout << "out::\n";
		cout << "finish" << endl;


	}



	void Calc(Element<N + 1>* m_KNN,XParams* m_para ,string filePath,string fileName)
	{
		m_filepath = filePath;
		m_filename = fileName;


		KNN = m_KNN;
		para = m_para;


		// total_num 表示点的总数目
		int total_num = 1;
		for (int i = 0; i < N; ++i)
			total_num *= (KNN->dim[i]-para->patch_w);



		for (int i = 0; i < N; ++i)
		{
			if (i == 0)
				m_index[0] = 1;
			else
				m_index[i] = m_index[i - 1] * (KNN->dim[i - 1]-para->patch_w);
		}


		random_unweighted_test(total_num);
	}


	int GetRowIndex(int *index)
	{
		int sum = index[0];
		for (int i = 1; i < N; ++i)
		{
			sum += index[i] * m_index[i];
		}

		return sum;
	}
	




	// 第一个参数发现没用！！！
	void random_unweighted_test(int n)
	{
		Graph g(n);
		
		int total_num = 1;
		for (int i = 0; i < N; ++i)
			total_num *= (KNN->dim[i]);

		std::vector<Vertex> vertices(total_num);
		// 初始化点的索引
		{
			typename boost::graph_traits<Graph>::vertex_iterator v, v_end;

			boost::tie(v, v_end) = boost::vertices(g);
			int index;

			int m_dim[N + 1];	//记录各个维度的数目
			for (int i = 0; i < N; ++i)
			{
				m_dim[i] = KNN->dim[i] - para->patch_w;
			}

			int m_iter[N + 1];// 用来循环
			if (N == 2)
				m_dim[2] = 1;

			int index_a[N + 1];// 当前点
			int index_b[N + 1];

			for (m_iter[2] = 0; m_iter[2] < m_dim[2]; ++m_iter[2])
			{
				for (m_iter[1] = 0; m_iter[1] < m_dim[1]; ++m_iter[1])
				{
					for (m_iter[0] = 0; m_iter[0] < m_dim[0]; ++m_iter[0])
					{
						for (int i = 0; i < N; ++i)
							index_a[i] = m_iter[i];

						index = GetRowIndex(index_a);

						put(boost::vertex_index, g, *v, index);
						vertices[index] = *v;
						++v;
					}
				}
			}



			std::vector<Edge> edges(n*para->knn);

			int cur_point = 0;
			int map_point = 0;
			for (m_iter[2] = 0; m_iter[2] < m_dim[2]; ++m_iter[2])
			{
				for (m_iter[1] = 0; m_iter[1] < m_dim[1]; ++m_iter[1])
				{
					for (m_iter[0] = 0; m_iter[0] < m_dim[0]; ++m_iter[0])
					{
						for (int i = 0; i < N; ++i)
							index_a[i] = m_iter[i];

						cur_point = GetRowIndex(index_a);

						for (int i = 0; i < para->knn; ++i)
						{
							index_a[N] = i;
							auto nn_ptr = KNN->get(index_a);

							for (int ii = 0; ii < N; ++ii)
								index_b[ii] = IntToIndex<N>(*nn_ptr, ii);
							
							map_point = GetRowIndex(index_b);

							boost::add_edge(vertices[cur_point], vertices[map_point], g);						
						}
					}
				}
			}
		}

		
		std::cout << "Random graph with " << n << " vertices and " << num_edges(g) << " edges.\n";

		std::cout << "  Direct translation of Brandes' algorithm...";
		std::vector<double> centrality(n);
		simple_unweighted_betweenness_centrality(g, get(boost::vertex_index, g), make_iterator_property_map(centrality.begin(), get(boost::vertex_index, g), double()));
		std::cout << "DONE.\n";
	}


	
	template<typename MutableGraph>
	void randomly_add_edges(MutableGraph& g, double edge_probability)
	{

		typedef typename boost::graph_traits<MutableGraph>::directed_category directed_category;
		const bool is_undirected = is_same<directed_category, boost::undirected_tag>::value;

		boost::minstd_rand gen;
		boost::uniform_01<boost::minstd_rand, double> rand_gen(gen);

		typedef typename boost::graph_traits<MutableGraph>::vertex_descriptor vertex;
		typename boost::graph_traits<MutableGraph>::vertex_iterator vi, vi_end;


		for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
		{
			vertex v = *vi;
			typename boost::graph_traits<MutableGraph>::vertex_iterator wi = is_undirected ? vi : vertices(g).first;

			while (wi != vi_end)
			{
				vertex w = *wi++;
				if (v != w) {
					if (rand_gen() < edge_probability)
						add_edge(v, w, g);
				}
			}
		}
	}




	template<typename Graph, typename VertexIndexMap, typename CentralityMap>
	void simple_unweighted_betweenness_centrality(const Graph& g, VertexIndexMap index, CentralityMap centrality)
	{
		
		typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex;
		typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
		typedef typename boost::graph_traits<Graph>::adjacency_iterator adjacency_iterator;
		typedef typename boost::graph_traits<Graph>::vertices_size_type vertices_size_type;
		typedef typename boost::property_traits<CentralityMap>::value_type centrality_type;

		vertex_iterator vi, vi_end;
		for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi)
		{
			boost::put(centrality, *vi, 0);
		}

		vertex_iterator si, si_end;

		for (boost::tie(si, si_end) = vertices(g); si != si_end; ++si) {
			int m_dist = distance(si, si_end);
			cout << "dist: " << m_dist << endl;

			vertex s = *si;

			// S <-- empty stack
			std::stack<vertex> S;

			// P[w] <-- empty list, w \in V
			typedef std::vector<vertex> Predecessors;

			
			std::vector<Predecessors> predecessors(num_vertices(g));
		

			// sigma[t] <-- 0, t \in V
			std::vector<vertices_size_type> sigma(num_vertices(g), 0);

			// sigma[s] <-- 1
			sigma[get(index, s)] = 1;

			// d[t] <-- -1, t \in V
			std::vector<int> d(num_vertices(g), -1);

			// d[s] <-- 0
			d[get(index, s)] = 0;

			// Q <-- empty queue
			std::queue<vertex> Q;

			// enqueue s --> Q
			Q.push(s);


			while (!Q.empty()) {
				// dequeue v <-- Q
				vertex v = Q.front(); Q.pop();

				// push v --> S
				S.push(v);

				adjacency_iterator wi, wi_end;
				for (boost::tie(wi, wi_end) = adjacent_vertices(v, g); wi != wi_end; ++wi) {


					vertex w = *wi;

					// w found for the first time?
					if (d[get(index, w)] < 0) {
						// enqueue w --> Q
						Q.push(w);

						// d[w] <-- d[v] + 1
						d[get(index, w)] = d[get(index, v)] + 1;
					}

					// shortest path to w via v?
					if (d[get(index, w)] == d[get(index, v)] + 1) {
						// sigma[w] = sigma[w] + sigma[v]
						sigma[get(index, w)] += sigma[get(index, v)];

						// append v --> P[w]
						predecessors[get(index, w)].push_back(v);
					}
				}
			}

			// delta[v] <-- 0, v \in V
			std::vector<centrality_type> delta(num_vertices(g), 0);

			// S returns vertices in order of non-increasing distance from s
			while (!S.empty()) {
				// pop w <-- S
				vertex w = S.top(); S.pop();

				const Predecessors& w_preds = predecessors[get(index, w)];
				for (typename Predecessors::const_iterator vi = w_preds.begin();
					vi != w_preds.end(); ++vi) {
					vertex v = *vi;
					// delta[v] <-- delta[v] + (sigma[v]/sigma[w])*(1 + delta[w])
					delta[get(index, v)] +=
						((centrality_type)sigma[get(index, v)] / sigma[get(index, w)])
						* (1 + delta[get(index, w)]);
				}

				if (w != s) {
					// C_B[w] <-- C_B[w] + delta[w]
					centrality[w] += delta[get(index, w)];
				}
			}
		}

		/*typedef typename boost::graph_traits<Graph>::directed_category directed_category;
		const bool is_undirected = is_same<directed_category, boost::undirected_tag>::value;
		if (is_undirected) {
			vertex_iterator v, v_end;
			for (boost::tie(v, v_end) = vertices(g); v != v_end; ++v) {
				put(centrality, *v, get(centrality, *v) / centrality_type(2));
			}
		}*/



		ofstream SaveNetFile(m_filepath + "//" + m_filename+"_Betweeness.data");
		for (boost::tie(si, si_end) = vertices(g); si != si_end; ++si) {
			SaveNetFile << centrality[*si];

			if (distance(si,si_end)!=1)
				SaveNetFile << endl;
		}


		SaveNetFile.close();


	}












};



#endif


