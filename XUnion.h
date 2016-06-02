#ifndef DUNION_H
#define DUNION_H
#include <utility>
#include <iostream>

using namespace std;
class XUnion
{
public:
	int count;
	int num_group;
	int *M;


	XUnion(int N) :count(N), num_group(N)
	{
		M = new int[N];
		if (M == NULL)
		{
			cout << "error!!" << endl;
			return;
		}
		for (int i = 0; i < N; i++)
			M[i] = i;
	}

	~XUnion()
	{
		delete[] M;
	}


	bool connectd(int p, int q)
	{
		return find(p) == find(q);
	}



	int find(int p)
	{
		while (p != M[p])
		{
			M[p] = M[M[p]];
			p = M[p];
		}
		return p;
	}



	void xUnion(int p, int q)
	{
		int pRoot = find(p);
		int qRoot = find(q);

		if (pRoot == qRoot)
			return;
		M[pRoot] = qRoot;
		num_group--;
	}

};



#endif
