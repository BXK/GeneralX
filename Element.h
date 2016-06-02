#include <iostream>


#include <string>
using namespace std;


#ifndef ELEMENT_H
#define ELEMENT_H


class XParams
{
public:
	int patch_w;
	int knn;
	int nn_iters;
	bool ToSelf;
	float point_dist;
};




// distance nearest neighbor
template <class T>
class knn_pair
{
public:
	T err, index;

	knn_pair(const T dist_, const T pt_) : err(dist_), index(pt_)
	{
	}

	int operator <(const knn_pair& b) const
	{
		return err < b.err;
	}

	int operator ==(const knn_pair& b) const
	{
		return err == b.err;
	}
};











template <int N>
class Element
{
public:
	int dim[N];
	int m_index[N];
	int num;
	void* data;

	explicit Element(int* m_dim)
	{
		num = N;
		Create_Volume(m_dim);
	}

	~Element()
	{
		delete[] data;
	}


	void CopyData(Element& II)
	{
		if (II.num != num)
		{
			cout << "copy error!" << endl;;
			system("pause");
			return;
		}
		delete[] data;

		Create_Volume(II.dim);
	}


	int getRowIndex(int *index)
	{
		try
		{
			for (int i = 0; i < N; ++i)
			{
				if (index[i]>dim[i])
				{
					throw 1;
				}
			}
		}
		catch (...)
		{
			cout << "error! index out of range!" << endl;
			system("pause");
		}



		

		int sum = index[0];
		for (int i = 1; i < N; ++i)
		{
			sum += index[i] * m_index[i];
		}

		return sum;
	}




	int* get(int* index)
	{
		/*try
		{
			for (int i = 0; i < N; ++i)
			{
				if (index[i]>dim[i])
				{
					throw 1;
				}
			}
		}
		catch (...)
		{
			cout << "error! index out of range!" << endl;
			system("pause");
		}
		


		int* tp = (int*)data;

		int sum = index[0];
		for (int i = 1; i < N; ++i)
		{
			sum += index[i] * m_index[i];
		}*/

		int* tp = (int*)data;
		int sum = getRowIndex(index);
		return &tp[sum];
	}


	void clear2val(int val)
	{
		int sum = 1;
		for (int i = 0; i < N; ++i)
		{
			sum *= dim[i];
		}


		int* tp = (int*)data;
		for (int i = 0; i < sum; ++i)
		{
			tp[i] = val;
		}
	}


	inline void LoadVolume(std::string filepath, int flag)
	{
		auto fp = fopen(filepath.c_str(), "rb");


		int num_total = 1;
		for (int i = 0; i < N; ++i)
		{
			num_total *= dim[i];
		}

		auto vol_data = reinterpret_cast<int*>(data);
		if (flag == 0)// char
		{
			unsigned char* m_data = new unsigned char[num_total];
			fread(m_data, sizeof(char), num_total, fp);

			for (auto i = 0; i < num_total; i++)
				vol_data[i] = static_cast<int>(m_data[i]);
		}
		else if (flag == 1)// int
		{
			fread((int*)data, sizeof(int), num_total, fp);
		}
		fclose(fp);
	}


	//@param flag 0 unsigned char, 1 int
	// filepath: C://data
	// filename aa
	inline void WriteVolume(std::string filePath, string filename, int flag)
	{
		auto raw_path = filePath + "//" + filename + ".raw";
		auto fp_raw = fopen(raw_path.c_str(), "wb");
		int num_total = 1;
		for (int i = 0; i < N; ++i)
		{
			num_total *= dim[i];
		}


		if (flag == 0)// char
		{
			unsigned char ch = 0;
			auto vol_data = reinterpret_cast<int*>(data);
			for (auto i = 0; i < num_total; i++)
			{
				auto t_data = vol_data[i];
				ch = static_cast<unsigned char>(t_data);
				fputc(ch, fp_raw);
			}
			fclose(fp_raw);
		}
		else if (flag == 1)// int
		{
			auto vol_data = (unsigned char*)(data);
			//cout << "me: " << vol_data[0] << "," << vol_data[1] << endl;
			for (auto i = 0; i < 4*num_total; i++)
			{
				fputc(vol_data[i], fp_raw);
			}
			fclose(fp_raw);
		}
	}


private:
	void Create_Volume(int* m_dim)
	{
		

		int sum = 1;
		for (int i = 0; i < N; ++i)
		{
			sum *= m_dim[i];
			dim[i] = m_dim[i];
			if (i == 0)
				m_index[0] = 1;
			else
				m_index[i] = m_index[i - 1] * dim[i - 1];
		}

		data = (void *)new int[sum];
		memset(data, 0, sizeof(int) * sum);
	}
};

















template<int N>
bool DistPoint(int *cord1, int *cord2, XParams *para)
{
	long long sum = 0;
	for (int i = 0; i < N;++i)
	{
		sum += (cord1[i] - cord2[i])*(cord1[i] - cord2[i]);
	}

	if (sum > para->point_dist)
		return true;
	return false;
}




template<int N>
int IndexToInt(int *cord)
{
	if (N > 3||N<2)
	{
		cout << "error" << endl;
		system("pause");
	}

	if (N==2)
	{
		return ((cord[1] << 10) | cord[0]);
	}

	if (N==3)
	{
		return ((cord[2] << 20) | (cord[1] << 10) | cord[0]);
	}

}


template<int N>
int IntToIndex(int val, int i)
{
	if (i < 0||i>N-1)
	{
		cout << "error" << endl;
		system("pause");
	}

	switch (i)
	{
	case 0:
		return ((val)&((1 << 10) - 1));
	case 1:
		return (((val) >> 10)&((1 << 10) - 1));
	case 2:
		return ((val) >> 20);
	}
}






string GetParams(XParams *para)
{
	string m_knn = to_string(para->knn);
	string m_pw = to_string(para->patch_w);
	string m_iter = to_string(para->nn_iters);
	string m_pd = to_string(int(para->point_dist));
	return "_k=" + m_knn + "_p" + m_pw + "_i" + m_iter + "_d" + m_pd;
}













#endif
