#include "Element.h"
#ifndef PATCH_H
#define PATCH_H

template <int N>
inline int patch_dist_ab(Element<N>* e1, Element<N>* e2, int* index1, int* index2, int maxVal, XParams* para)
{
	auto ans = 0;

	int m_index1[N];
	int m_index2[N];
	int iter_index[N];

	if (N == 3)
	{
		for (iter_index[2] = 0; iter_index[2] < para->patch_w; ++iter_index[2])
		{
			for (iter_index[1] = 0; iter_index[1] < para->patch_w; ++iter_index[1])
			{
				for (iter_index[0] = 0; iter_index[0] < para->patch_w; ++iter_index[0])
				{
					for (int ii = 0; ii < N; ++ii)
					{
						m_index1[ii] = iter_index[ii] + index1[ii];
						m_index2[ii] = iter_index[ii] + index2[ii];
					}
					auto c1 = *(e1->get(m_index1));
					auto c2 = *(e2->get(m_index2));

					auto dr = (c1 & 255) - (c2 & 255);
					ans += dr * dr;
					if (ans >= maxVal)
					{
						return ans;
					}
				}
			}
		}
	}
	else if (N == 2)
	{
		for (iter_index[1] = 0; iter_index[1] < para->patch_w; ++iter_index[1])
		{
			for (iter_index[0] = 0; iter_index[0] < para->patch_w; ++iter_index[0])
			{
				for (int ii = 0; ii < N; ++ii)
				{
					m_index1[ii] = iter_index[ii] + index1[ii];
					m_index2[ii] = iter_index[ii] + index2[ii];
				}
				auto c1 = *e1->get(m_index1);
				auto c2 = *e2->get(m_index2);

				auto dr = (c1 & 255) - (c2 & 255);
				ans += dr * dr;
				if (ans >= maxVal)
				{
					return ans;
				}
			}
		}
	}
	else
	{
		printf("error dimension error!\n");
		system("pause");
	}


	return ans;
}


#endif

