#ifndef GENERALANALYSIS_H
#define GENERALANALYSIS_H

#include "Element.h"
#include "XUnion.h"
#include <map>
#include <string>
#include <fstream>
#include <stack>
#include <queue>
#include <list>
#include "Betweenness.h"
using namespace std;

template <int N>// n for dimenstion
class GeneralAnalysis
{
public:
	void SetData(Element<N + 1>* m_knn, Element<N + 1>* m_d_knn, XParams* m_para, Element<N>* m_e1,Element<N>* m_e2)
	{
		KNN = m_knn;
		D_KNN = m_d_knn;
		para = m_para;
		e1 = m_e1;
		e2 = m_e2;
	}



	// 将一个区域的点集合并！ 这个函数只能对一个element 施行
	void AnalyzeChain(string filePath, string fileName)
	{
		int m_dim[N + 1];	//记录各个维度的数目
		for (int i = 0; i < N ; ++i)
		{
			m_dim[i] = KNN->dim[i] - para->patch_w;
		}

		int total_num = 1;
		for (int i = 0; i < N; ++i)
			total_num *= KNN->dim[i];

		XU = new XUnion(total_num);

		int m_iter[N + 1];// 用来循环
		if (N == 2)
			m_dim[2] = 1;


		int index_a[N + 1];// 当前点
		int index_b[N + 1];

		int XUindex_a, XUindex_b;

		for (m_iter[2] = 0; m_iter[2] < m_dim[2]; ++m_iter[2])
		{
			for (m_iter[1] = 0; m_iter[1] < m_dim[1]; ++m_iter[1])
			{
				for (m_iter[0] = 0; m_iter[0] < m_dim[0]; ++m_iter[0])
				{
					for (int i = 0; i < N; ++i)
						index_a[i] = m_iter[i];				

					XUindex_a = e1->getRowIndex(index_a);
					

					for (int i = 0; i < para->knn; ++i)
					{
						index_a[N] = i;
						auto nn_ptr = KNN->get(index_a);
						
						for (int ii = 0; ii < N; ++ii)
							index_b[ii] = IntToIndex<N>(*nn_ptr, ii);					

						XUindex_b = e1->getRowIndex(index_b);

						XU->xUnion(XUindex_a, XUindex_b);
					}
				}
			}
		}
		OutPutUnion(filePath,fileName);
	}

	void OutPutUnion(string filePath,string fileName)
	{
		int m_dim[N + 1];	//记录各个维度的数目
		for (int i = 0; i < N; ++i)
		{
			m_dim[i] = KNN->dim[i] - para->patch_w;
		}

		Element<N> OutUnion(KNN->dim);


		int m_iter[N + 1];// 用来循环
		if (N == 2)
			m_dim[2] = 1;




		int index_a[N + 1];// 当前点
		map<int, int> MM;
		int num_cout = 0;
		


		for (m_iter[2] = 0; m_iter[2] < m_dim[2]; ++m_iter[2])
		{
			for (m_iter[1] = 0; m_iter[1] < m_dim[1]; ++m_iter[1])
			{
				for (m_iter[0] = 0; m_iter[0] < m_dim[0]; ++m_iter[0])
				{
					for (int i = 0; i < N; ++i)
						index_a[i] = m_iter[i];
					int index_t = e1->getRowIndex(index_a);
					int Parent = XU->find(index_t);

					if(MM.count(Parent)==0)
					{
						*(OutUnion.get(index_a)) = num_cout;
						MM.insert(make_pair(Parent, num_cout++));						
					}
					else
					{

						*(OutUnion.get(index_a)) = MM.count(Parent);
					}
				}
			}
		}

		//string filepath = "G://GeneralizePatchMatch//PatchMatch//KNN01//knn01//GeneralX";
		OutUnion.WriteVolume(filePath, fileName+"_Union", 1);







		// 输出有向图

		ofstream SaveNetFile(filePath+"//"+fileName+"_Net.gml");
		SaveNetFile << "graph\n[\n\tdirected 1\n";
		for (m_iter[2] = 0; m_iter[2] < m_dim[2]; ++m_iter[2])
		{
			for (m_iter[1] = 0; m_iter[1] < m_dim[1]; ++m_iter[1])
			{
				for (m_iter[0] = 0; m_iter[0] < m_dim[0]; ++m_iter[0])
				{
					for (int i = 0; i < N; ++i)
						index_a[i] = m_iter[i];
					int index_t = e1->getRowIndex(index_a);
					SaveNetFile << "\tnode\n\t[\n\t\tid " << index_t << "\n\t]\n";
				}
			}
		}



		int index_b[N + 1];

		int XUindex_a, XUindex_b;


		for (m_iter[2] = 0; m_iter[2] < m_dim[2]; ++m_iter[2])
		{
			for (m_iter[1] = 0; m_iter[1] < m_dim[1]; ++m_iter[1])
			{
				for (m_iter[0] = 0; m_iter[0] < m_dim[0]; ++m_iter[0])
				{
					for (int i = 0; i < N; ++i)
						index_a[i] = m_iter[i];

					XUindex_a = e1->getRowIndex(index_a);


					for (int i = 0; i < para->knn; ++i)
					{
						index_a[N] = i;
						auto nn_ptr = KNN->get(index_a);

						for (int ii = 0; ii < N; ++ii)
							index_b[ii] = IntToIndex<N>(*nn_ptr, ii);
						XUindex_b = e1->getRowIndex(index_b);

						//XU->xUnion(XUindex_a, XUindex_b);
						SaveNetFile << "\tedge\n\t[\n\t\tsource " << XUindex_a << "\n\t\ttarget\t"<<XUindex_b<<"\n\t]\n";
					}
				}
			}
		}

		bool flag = true;



		// 四邻域
		if (flag)
		{
			if (N == 2)
				m_dim[2] = 3;


			for (m_iter[2] = 1; m_iter[2] < m_dim[2]-1; ++m_iter[2])
			{
				for (m_iter[1] = 1; m_iter[1] < m_dim[1]-1; ++m_iter[1])
				{
					for (m_iter[0] = 1; m_iter[0] < m_dim[0]-1; ++m_iter[0])
					{
						for (int i = 0; i < N; ++i)
							index_a[i] = m_iter[i];

						XUindex_a = e1->getRowIndex(index_a);
						

						for (int i = 0; i < N;++i)
						{
							index_a[i] = index_a[i] + 1;

							XUindex_b = e1->getRowIndex(index_a);
							SaveNetFile << "\tedge\n\t[\n\t\tsource " << XUindex_a << "\n\t\ttarget\t" << XUindex_b << "\n\t]\n";

							index_a[i] = index_a[i] - 1;
						}

						for (int i = 0; i < N; ++i)
						{
							index_a[i] = index_a[i] - 1;

							XUindex_b = e1->getRowIndex(index_a);
							SaveNetFile << "\tedge\n\t[\n\t\tsource " << XUindex_a << "\n\t\ttarget\t" << XUindex_b << "\n\t]\n";

							index_a[i] = index_a[i] + 1;
						}

					}
				}
			}
		}


		SaveNetFile << "]";
		SaveNetFile.close();


		// 分析入度！！
		

		Element<N> InDegree(KNN->dim);
		InDegree.clear2val(0);
		index_b[N + 1];

		for (m_iter[2] = 0; m_iter[2] < m_dim[2]; ++m_iter[2])
		{
			for (m_iter[1] = 0; m_iter[1] < m_dim[1]; ++m_iter[1])
			{
				for (m_iter[0] = 0; m_iter[0] < m_dim[0]; ++m_iter[0])
				{
					for (int i = 0; i < N; ++i)
						index_a[i] = m_iter[i];


					for (int i = 0; i < para->knn; ++i)
					{
						index_a[N] = i;
						auto nn_ptr = KNN->get(index_a);

						for (int ii = 0; ii < N; ++ii)
							index_b[ii] = IntToIndex<N>(*nn_ptr, ii);

						*InDegree.get(index_b) = *InDegree.get(index_b)+1;
					}
				}
			}
		}

		InDegree.WriteVolume(filePath, fileName+"_InDegree", 1);
	}


	void BoostBetweenness(string filePath,string fileName)
	{
		CalcBetweenness<N> CBB;
		CBB.Calc(KNN, para,filePath,fileName);		
	}






private:

	Element<N>* e1;
	Element<N>* e2;
	Element<N + 1>* KNN;
	Element<N + 1>* D_KNN;
	XParams* para;
	XUnion* XU;
};


#endif

