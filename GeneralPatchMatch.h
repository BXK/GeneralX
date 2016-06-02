#include "Element.h"
#include <set>
#include <vector>
#include "Patch.h"
#include <algorithm>
#include <ctime>
using namespace std;


#ifndef GENERALPATCHMATCH_H
#define GENERALPATCHMATCH_H


template <int N>// n for dimenstion
class GeneralPatch
{
private:
	XParams* para;
	Element<N>* e1;
	Element<N>* e2;
	Element<N + 1>* KNN;
	Element<N + 1>* D_KNN;

public:
	void SetData(Element<N>* v1, Element<N>* v2, XParams* mpara)
	{
		e1 = v1;
		e2 = v2;
		para = mpara;
	}

	void DoGeneralPatch()
	{
		Init_Patch_Match();
		Init_P_M_Dist();
		Do_Patch_Match();
	}

	void SaveData(string filepath,string fileName)
	{
		KNN->WriteVolume(filepath, fileName+"_knn", 1);
		D_KNN->WriteVolume(filepath, fileName + "_d_knn", 1);
	}




	void Replace(string filepath,string fileName)
	{
		//// 替换
		int index_a[N + 1];
		int index_b[N + 1];

		int iter_index[N];


		int* dim_1 = new int[e1->num + 1];

		for (int i = 0; i < e1->num; ++i)
		{
			dim_1[i] = e1->dim[i] - para->patch_w;
		}


		for (int iter = 0; iter < para->knn; ++iter)
		{
			e1->clear2val(0);

			if (N == 2)
			{
				dim_1[2] = 1;
			}

			for (iter_index[2] = 0; iter_index[2] < dim_1[2]; ++iter_index[2])
			{				
				for (iter_index[1] = 0; iter_index[1] < dim_1[1]; ++iter_index[1])
				{
					for (iter_index[0] = 0; iter_index[0] < dim_1[0]; ++iter_index[0])
					{
						for (int ii = 0; ii < N; ++ii)
						{
							index_a[ii] = iter_index[ii];
						}
						index_a[N] = iter;

						auto knn_ptr = KNN->get(index_a);
						for (int i = 0; i < N; i++)
							index_b[i] = IntToIndex<N>(*knn_ptr, i);

						*e1->get(index_a) = *e2->get(index_b);
					}
				}
			}
			string num = to_string(iter);

			string ff = GetParams(para);
			e1->WriteVolume(filepath, fileName + "_replace" + num, 0);
		}
	}


private:

	// init the knn get the ann
	void Init_Patch_Match()
	{
		int* e_dim = new int[e1->num + 1];
		for (int i = 0; i < e1->num; ++i)
			e_dim[i] = e1->dim[i];
		e_dim[e1->num] = para->knn;


		KNN = new Element<N + 1>(e_dim);

		int* dim_1 = new int[e1->num + 1];
		int* dim_2 = new int[e1->num];

		for (int i = 0; i < e1->num; ++i)
		{
			dim_1[i] = e1->dim[i] - para->patch_w;
			dim_2[i] = e2->dim[i] - para->patch_w;
		}


		int m_index[N + 1];
		int new_index[N];

		int* iter_index = new int[N];


		if (N == 2)
		{
			dim_1[2] = 1;
		}
		for (iter_index[2] = 0; iter_index[2] < dim_1[2]; ++iter_index[2])
		{
			for (iter_index[1] = 0; iter_index[1] < dim_1[1]; ++iter_index[1])
			{
				for (iter_index[0] = 0; iter_index[0] < dim_1[0]; ++iter_index[0])
				{
					for (int jj = 0; jj < N; ++jj)
						m_index[jj] = iter_index[jj];

					set<int> SS;
					for (auto i = 0; i < para->knn; ++i)
					{
						for (;;)
						{
							// 随机坐标
							for (int j = 0; j < N; ++j)
								new_index[j] = rand() % dim_2[j];

							if (para->ToSelf)
							{
								if (!DistPoint<N>(m_index, new_index, para))
									continue;
							}

							// 这个没有判断 点之间的距离
							if (SS.count(IndexToInt<N>(new_index)) == 0)
							{
								m_index[N] = i;
								SS.insert(IndexToInt<N>(new_index));
								*KNN->get(m_index) = IndexToInt<N>(new_index);
								break;
							}
						}
					}
				}
			}
		}
	}

	void Init_P_M_Dist()
	{
		cout << "Calc Init Dist" << endl;

		if (para->patch_w > 32)
		{
			printf("patch width unsupported!!!");
			return;
		}

		int* e_dim = new int[e1->num + 1];
		for (int i = 0; i < e1->num; ++i)
			e_dim[i] = e1->dim[i];
		e_dim[e1->num] = para->knn;

		D_KNN = new Element<N + 1>(e_dim);
		D_KNN->clear2val(INT_MAX);

		int* dim_1 = new int[e1->num + 1];
		int* dim_2 = new int[e1->num];

		for (int i = 0; i < e1->num; ++i)
		{
			dim_1[i] = e1->dim[i] - para->patch_w;
			dim_2[i] = e2->dim[i] - para->patch_w;
		}

		int index_a[N + 1];
		int index_b[N];

		int* iter_index = new int[N + 1];


		vector<knn_pair<int>> v;
		v.reserve(para->knn);


		if (N == 2)
		{
			dim_1[2] = 1;
		}


		for (iter_index[2] = 0; iter_index[2] < dim_1[2]; ++iter_index[2])
		{
			for (iter_index[1] = 0; iter_index[1] < dim_1[1]; ++iter_index[1])
			{
				for (iter_index[0] = 0; iter_index[0] < dim_1[0]; ++iter_index[0])
				{
					for (int jj = 0; jj < N; ++jj)
						index_a[jj] = iter_index[jj];


					v.clear();

					for (auto i = 0; i < para->knn; i++)
					{
						index_a[N] = i;
						auto ann_ptr = KNN->get(index_a);

						for (int ii = 0; ii < N; ++ii)
							index_b[ii] = IntToIndex<N>(*ann_ptr, ii);

						auto d = patch_dist_ab(e1, e2, index_a, index_b, INT_MAX / 2, para);

						v.push_back(knn_pair<int>(d, *ann_ptr));
					}

					make_heap(v.begin(), v.end());


					for (auto i = 0; i < para->knn; i++)
					{
						index_a[N] = i;
						auto dnn_ptr = D_KNN->get(index_a);
						auto ann_ptr = KNN->get(index_a);

						auto current = v[i];
						*dnn_ptr = current.err;
						*ann_ptr = current.index;
					}
				}
			}
		}
	}



	// index1 是源点 index2 是index1 匹配的点
	void Attemp_P_M(vector<knn_pair<int>>& q, int* index1, int* index2)
	{
		if (q.size() != para->knn)
		{
			printf("q size is wrong (%d, %d)\n", static_cast<int>(q.size()), para->knn);
			exit(1);
		}

		/*int num = 0;
		int index_tp[N];
		for (auto iter = m_hash.begin(); iter != m_hash.end(); ++iter)
		{
			int index_Int = *iter;
			for (int i = 0; i < N; ++i)
				index_tp[i] = IntToIndex<N>(index_Int, i);

			if (!DistPoint<N>(index_tp, index2, para))
			{
				num++;
			}
			if (num > 1)
				return;
		}*/

		int num = 0;
		int index_tp[N];
		for (int i = 0; i < para->knn;++i)
		{
			int index_Int = q[i].index;
			for (int i = 0; i < N; ++i)
				index_tp[i] = IntToIndex<N>(index_Int, i);

			if (!DistPoint<N>(index_tp, index2, para))
			{
				num++;
			}
			if (num > 1)
				return;
		}





		if (para->ToSelf)
		{
			if (!DistPoint<N>(index1, index2, para))
				return;
		}



		auto pos = IndexToInt<N>(index2);


		// 这个处理可能有问题
		auto err = q[0].err;
		auto current = patch_dist_ab(e1, e2, index1, index2, err, para);

		if (current < err)
		{

			// 大根堆
			q.push_back(knn_pair<int>(current, pos));
			push_heap(q.begin(), q.end());
			pop_heap(q.begin(), q.end());
			q.pop_back();
		}
	}

	void Do_Patch_Match()
	{
		auto nn_iter = 0;
		for (; nn_iter < para->nn_iters; ++nn_iter)
		{
			double dur;
			clock_t start, end;
			start = clock();


			int Index_Start[N + 1];
			int Index_End[N + 1];
			int Index_Chan[N + 1];

			if (nn_iter % 2 == 0)
			{
				for (int i = 0; i < N; ++i)
				{
					Index_Start[i] = 0;
					Index_End[i] = e1->dim[i] - para->patch_w;
					Index_Chan[i] = 1;
				}
			}
			else
			{
				for (int i = 0; i < N; ++i)
				{
					Index_Start[i] = e1->dim[i] - para->patch_w;
					Index_End[i] = 0;
					Index_Chan[i] = -1;
				}
			}


			vector<knn_pair<int>> v;
			v.reserve(para->knn+1);
			for (auto i = 0; i < para->knn; i++)
				v.push_back(knn_pair<int>(0, 0));


			int* iter_index = new int[N + 1];


			int index_a[N + 1];


			if (N == 2)
			{
				Index_Start[2] = 0;
				Index_End[2] = 1;
				Index_Chan[2] = 1;
			}

			for (iter_index[2] = Index_Start[2]; iter_index[2] != Index_End[2]; iter_index[2] += Index_Chan[2])
			{
				for (iter_index[1] = Index_Start[1]; iter_index[1] != Index_End[1]; iter_index[1] += Index_Chan[1])
				{
					for (iter_index[0] = Index_Start[0]; iter_index[0] != Index_End[0]; iter_index[0] += Index_Chan[0])
					{
						for (int i = 0; i < N; ++i)
							index_a[i] = iter_index[i];


						//set<int> m_Hash;

						for (int ii = 0; ii < para->knn; ++ii)
						{
							index_a[N] = ii;
							int* knn_ptr = KNN->get(index_a);
							int* dknn_ptr = D_KNN->get(index_a);

							v[ii] = knn_pair<int>(*dknn_ptr, *knn_ptr);
							//m_Hash.insert(*knn_ptr);
						}


						/*********************************************************************************
						*
						*						propagate
						*
						*
						**********************************************************************************/
						int index_a_tp[N + 1];
						int index_b_tp[N + 1];


						for (auto tt = 0; tt < N; ++tt)
						{
							if (static_cast<unsigned>(iter_index[tt] + Index_Chan[tt]) < static_cast<unsigned>(e1->dim[tt] - para->patch_w))
							{
								for (auto ii = 0; ii < para->knn; ++ii)
								{
									for (int i = 0; i < N; ++i)
										index_a_tp[i] = iter_index[i];

									index_a_tp[tt] += Index_Chan[tt];
									index_a_tp[N] = ii;

									auto tp_ptr = KNN->get(index_a_tp);

									for (int i = 0; i < N; ++i)
										index_b_tp[i] = IntToIndex<N>(*tp_ptr, i);


									index_b_tp[tt] -= Index_Chan[tt];


									if (static_cast<unsigned>(index_b_tp[tt]) < static_cast<unsigned>(e1->dim[tt] - para->patch_w))
										Attemp_P_M(v, index_a, index_b_tp);
								}
							}
						}


						// random search
						int rs_start = INT_MAX;
						int max_rang = -1;
						for (int i = 0; i < N; ++i)
							if (e2->dim[i] > max_rang)
								max_rang = e2->dim[i];
						if (rs_start > max_rang)
							rs_start = max_rang;


						int BestCord[N + 1];
						int CordMin[N];
						int CordMax[N];
						int Cord_b[N];


						for (auto mag = rs_start; mag >= 1; mag /= 2)
						{
							for (auto ii = 0; ii < para->knn; ++ii)
							{
								for (int t = 0; t < N; ++t)
								{
									BestCord[t] = IntToIndex<N>(v[ii].index, t);
									CordMin[t] = BestCord[t] - mag > 0 ? BestCord[t] - mag : 0;
									CordMax[t] = BestCord[t] + mag > e2->dim[t] - para->patch_w ? e2->dim[t] - para->patch_w : BestCord[t] + mag;
									Cord_b[t] = CordMin[t] + rand() % (CordMax[t] - CordMin[t]);
								}
								Attemp_P_M(v, index_a, Cord_b);
							}
						}

						sort(v.begin(), v.end());
						for (auto i = 0; i < para->knn; i++)
						{
							index_a[N] = i;
							int* knn_ptr = KNN->get(index_a);
							int* dknn_ptr = D_KNN->get(index_a);
							*knn_ptr = v[i].index;
							*dknn_ptr = v[i].err;
						}
					}
				}
			}


			end = clock();
			dur = static_cast<double>(end - start);
			printf("knn iter: %d, cost: %f\n", nn_iter, (dur / CLOCKS_PER_SEC));
		}
	}
};


#endif

