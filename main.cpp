#include <iostream>

#include "Element.h"
#include "GeneralPatchMatch.h"
#include <string>
#include "GeneralAnalysis.h"


string filename = "sugar_499¡Á499_2D_unit8";
string filepath1 = "G://GeneralizePatchMatch//PatchMatch//KNN01//Matlab2DpicGenerate//RawData//2D//" + filename + ".raw";
string out_PatchMatch_Path = "G://GeneralizePatchMatch//PatchMatch//GeneralX//GeneralX//GeneralX//Result//PatchMatch";
string out_Replaced_Path = "G://GeneralizePatchMatch//PatchMatch//GeneralX//GeneralX//GeneralX//Result//Replaced";
string out_Net_Path = "G://GeneralizePatchMatch//PatchMatch//GeneralX//GeneralX//GeneralX//Result//Net";





int const DIM = 2;
XParams *paras = new XParams;
Element<DIM> *TT1;
Element<DIM> *TT2;
int m2d[4] = { 499, 499, 499 };



Element<DIM + 1>* KNN;
Element<DIM + 1>* D_KNN;



void Init()
{
	TT1 = new Element<DIM>(m2d);
	TT1->clear2val(10);
	TT1->LoadVolume(filepath1, 0);


	TT2 = new Element<DIM>(m2d);
	TT2->clear2val(10);
	TT2->LoadVolume(filepath1, 0);
}


void Do_Patch_Match()
{
	GeneralPatch<DIM> GG;
	GG.SetData(TT1, TT2, paras);
	GG.DoGeneralPatch();

	string ff = GetParams(paras);
	string m_file = filename + ff;

	GG.Replace(out_Replaced_Path, m_file);
	GG.SaveData(out_PatchMatch_Path, m_file);
}


void LoadData()
{
	m2d[DIM] = paras->knn;
	KNN = new Element<DIM+1>(m2d);
	D_KNN = new Element<DIM + 1>(m2d);

	string ff = GetParams(paras);
	string m_file = filename + ff;

	KNN->LoadVolume(out_PatchMatch_Path + "//" + m_file + "_knn.raw", 1);
	D_KNN->LoadVolume(out_PatchMatch_Path + "//" + m_file + "_d_knn.raw", 1);


	GeneralAnalysis<DIM> GA;
	GA.SetData(KNN, D_KNN, paras,TT1,TT2);

}

void Do_Analyze()
{
	LoadData();
	GeneralAnalysis<DIM>  AN;

	string ff = GetParams(paras);
	string m_file = filename + ff;

	AN.SetData(KNN, D_KNN, paras, TT1, TT2);
	AN.AnalyzeChain(out_Net_Path, m_file);
}


void DoBetweenness()
{
	LoadData();
	GeneralAnalysis<DIM>  AN;
	AN.SetData(KNN, D_KNN, paras, TT1, TT2);

	string ff = GetParams(paras);
	string m_file = filename + ff;

	AN.BoostBetweenness(out_Net_Path, m_file);
}





int main()
{
	paras->knn = 5;
	paras->patch_w = 1;
	paras->nn_iters = 60;
	paras->ToSelf = true;
	paras->point_dist = 400;



	Init();

	cout << "IF you want do KNN please enter y:" << endl;
	char ch;
	cin >> ch;

	if (ch == 'y')
		Do_Patch_Match();
	

	cout << "IF you want do AN please enter y:" << endl;
	cin >> ch;
	if (ch == 'y')
		Do_Analyze();


	cout << "IF you want do betweenness y:" << endl;
	cin >> ch;
	if (ch == 'y')
		DoBetweenness();



	system("pause");
}
