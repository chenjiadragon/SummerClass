// Lesson1-RW.cpp : 定义控制台应用程序的入口点。
#include "../include/MyBaseModel.hpp"
#include "../include/MyHalfEdgeModel.hpp"
#include "../include/Dijkstra.hpp"
#include "../include/FastMarching.hpp"
#include "../include/AddBase.hpp"

int main()
{
	//CGAL
	MyBaseModel model;
	model.ReadObjFile("../modeldata/sphere.obj");
	MyBaseModel final_model = AddBase(model, 0.1);
	final_model.WriteObjFile("../modeldata/final_model.obj");
	return 0;
}


int main3()
{
	MyHalfEdgeModel model;
	model.ReadObjFile("../modeldata/sphere.obj");


	//dij和fasmarching并纹理贴图
	set<int> sources;
	sources.insert(387);
	sources.insert(202);


	Dijkstra dij(model, sources);
	FastMarching fas(model, sources);
	dij.Run();
	fas.Run();
	auto disF = dij.GetDistanceField();
	auto disF_fas = fas.GetDistanceField();

	auto maxDis = *max_element(disF.begin(), disF.end());
	for (int i = 0; i < disF.size(); ++i)
		disF[i] /= maxDis * 1.01;
	auto maxFas = *max_element(disF_fas.begin(), disF_fas.end());
	for (int i = 0; i < disF_fas.size(); ++i)
		disF_fas[i] /= maxFas * 1.01;
	model.WriteTexturedObjFile("../modeldata/sphere_distanceField.obj", disF);
	model.WriteTexturedObjFile("../modeldata/sphere_distanceField_fas.obj", disF_fas);


	return 0;
}

int main2()
{
	MyHalfEdgeModel model;
	model.ReadObjFile("../modeldata/sphere.obj");

	//等值线
	vector<double> scalarField;
	for (auto v : model.GetVertices())
	{
		scalarField.push_back(v.z());
	}
	auto isoline = model.ExtractIsoline(scalarField, 0.1);
	model.SaveIsoline("../modeldata/isoline.obj", isoline);
	pair<MyBaseModel, MyBaseModel> ans = model.SplitModelByIsoline(scalarField, 0.1);
	ans.first.WriteObjFile("../modeldata/sphere_less.obj");
	ans.second.WriteObjFile("../modeldata/sphere_larger.obj");

	return 0;
}

int main1()
{
	MyHalfEdgeModel model;
	model.ReadObjFile("../modeldata/sphere.obj");
	model.ReadObjFile("../mdoeldata/bunny.obj");

	return 0;
}


