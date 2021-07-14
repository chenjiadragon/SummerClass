#pragma once
#include "MyBaseModel.hpp"
#include <cmath>
#include <queue>

class MyHalfEdgeModel :	public MyBaseModel
{
public:
	struct MyMeshEdge
	{
		int leftVert;
		int rightVert;
		int indexOfFrontFace;
		int indexOfOppositeVertex;
		int indexOfPreviousEdge;   //pre
		int indexOfNextEdge;  //next
		int indexOfReverseEdge;   //twin
		double leftAngle;        //根据向量求夹角
		double rightAngle;
		double length;
	};
protected:
	vector<MyMeshEdge> m_edges;
	map<pair<int, int>, int> fromEdge2ID;  //根据两个点找半边
	vector<int> fromVert2Edge;

public:
	void ReadObjFile(const char* filename);
	vector<int> GetNeighboringVertices(int v) const;
	void CreateMeshEdges();
	MyMeshEdge index_of_edge(int ID) //根据边的下标找边
	{
		return m_edges[ID];
	}
	int index_of_ver(int IDOfVer)   //根据点的下标找边
	{
		return fromVert2Edge[IDOfVer];
	}
	int index_of_face(int IDOfFace) //根据面的下标找边(其中一条)
	{
		return IDOfFace * 3;
	}
	int GetNumberOfComponents();
	vector<MyBaseModel> Decompose();

	//V-E+F=2-2g
	int GetNumberOfGenii() const;

};

inline void MyHalfEdgeModel::ReadObjFile(const char* filename)
{
	MyBaseModel::ReadObjFile(filename);
	CreateMeshEdges();
}

//得到顶点的相邻顶点
vector<int> MyHalfEdgeModel::GetNeighboringVertices(int v) const
{
	vector<int> neighbors;
	int firstEdge = fromVert2Edge[v];
	int nxtEdge = firstEdge;
	do
	{
		neighbors.push_back(m_edges[nxtEdge].rightVert);
		nxtEdge = m_edges[m_edges[nxtEdge].indexOfPreviousEdge].indexOfReverseEdge;
	} while (nxtEdge != firstEdge);
	return neighbors;
}

//V-E+F=2-2g
int MyHalfEdgeModel::GetNumberOfGenii() const
{
	return 1 - (m_verts.size() + m_faces.size() - m_edges.size() / 2) / 2;
}


void MyHalfEdgeModel::CreateMeshEdges()
{
	m_edges.resize(3 * m_faces.size());
	fromVert2Edge.resize(m_verts.size());
	for (int i = 0; i < m_faces.size(); ++i)
	{
		//MyMeshEdge temp;
		for (int j = 0; j < 3; j++)
		{
			//3*i+j
			m_edges[3 * i + j].indexOfFrontFace = i;
			m_edges[3 * i + j].leftVert = m_faces[i][j];
			fromVert2Edge[m_edges[3 * i + j].leftVert] = 3 * i + j;
			m_edges[3 * i + j].rightVert = m_faces[i][(j + 1) % 3];
			m_edges[3 * i + j].indexOfOppositeVertex = m_faces[i][(j + 2) % 3];
			m_edges[3 * i + j].indexOfPreviousEdge = 3 * i + (j + 3 - 1) % 3;
			m_edges[3 * i + j].indexOfNextEdge = 3 * i + (j + 3 + 1) % 3;
			fromEdge2ID[make_pair(m_edges[3 * i + j].leftVert, m_edges[3 * i + j].rightVert)] = 3 * i + j;
			//求边长
			m_edges[3 * i + j].length = (m_verts[m_edges[3 * i + j].leftVert] - m_verts[m_edges[3 * i + j].rightVert]).norm();


		}
	}
	for (int i = 0; i < m_edges.size(); i++)
	{
		m_edges[i].indexOfReverseEdge = fromEdge2ID[make_pair(m_edges[i].rightVert, m_edges[i].leftVert)];
	}	
}


int MyHalfEdgeModel::GetNumberOfComponents()
{
	set<int> unVisistedVertices;


	for (int i = 0; i < m_verts.size(); ++i)
	{
		unVisistedVertices.insert(i);
		VertsId2Component.push_back(0);
	}
		
	int numberOfComponents = 0;
	while (!unVisistedVertices.empty())
	{
		++numberOfComponents;
		auto firstVertex = *unVisistedVertices.begin();
		queue<int> tobedeleted;
		tobedeleted.push(firstVertex);

		VertsId2Component[firstVertex] = numberOfComponents;

		unVisistedVertices.erase(firstVertex);
		while (!tobedeleted.empty())
		{
			auto top = tobedeleted.front();
			tobedeleted.pop();
			auto firstEdge = fromVert2Edge[top];
			auto nxtEdge = firstEdge;
			do
			{
				if (unVisistedVertices.find(m_edges[nxtEdge].rightVert) != unVisistedVertices.end()) //若由v2e找到的边的右端点未被访问过
				{
					tobedeleted.push(m_edges[nxtEdge].rightVert);
					VertsId2Component[m_edges[nxtEdge].rightVert] = numberOfComponents;

					unVisistedVertices.erase(m_edges[nxtEdge].rightVert);
					nxtEdge = m_edges[m_edges[nxtEdge].indexOfPreviousEdge].indexOfReverseEdge; //找到这条边的前边的反边
				}
				else
				{
					nxtEdge = m_edges[m_edges[nxtEdge].indexOfPreviousEdge].indexOfReverseEdge; //找到这条边的前边的反边
				}
			} while (nxtEdge != firstEdge);
		}
	}
	return numberOfComponents;
}


vector<MyBaseModel> MyHalfEdgeModel::Decompose()
{
	int numberOfComponents = GetNumberOfComponents();  //连通分支数

	vector<vector<Eigen::Vector3d>> verts_multi_component;
	verts_multi_component.clear();
	verts_multi_component.resize(numberOfComponents + 1);
	vector<vector<Eigen::Vector3i>> faces_multi_component;
	faces_multi_component.clear();
	faces_multi_component.resize(numberOfComponents + 1);

	vector<MyBaseModel> resModels;

	vector<map<int,int>> vertmap(numberOfComponents);    //根据旧点找新点

	//vector<int> StartCounter;
	//StartCounter.push_back(0);
	for (int i = 0; i < m_verts.size(); i++)
	{
		verts_multi_component[VertsId2Component[i] - 1].push_back(m_verts[i]);
		//if (i > 0 && VertsId2Component[i] > VertsId2Component[i - 1])
		//{
		//	StartCounter.push_back(i);
		//}
		//都是从0开始的
		vertmap[VertsId2Component[i] - 1][i] = verts_multi_component[VertsId2Component[i] - 1].size() - 1;

	}

	for (int i = 0; i < m_faces.size(); i++)
	{
		int index = VertsId2Component[m_faces[i].x()] - 1;//面在哪个连通分支上
		Eigen::Vector3i new_face;
		new_face.x() = vertmap[index][m_faces[i].x()];
		new_face.x() = vertmap[index][m_faces[i].x()];
		new_face.x() = vertmap[index][m_faces[i].x()];
		//new_face.x() = m_faces[i].x();
		//new_face.y() = m_faces[i].y() - StartCounter[index];
		//new_face.z() = m_faces[i].z() - StartCounter[index];
		faces_multi_component[index].push_back(new_face);
	}

	for (int i = 0; i < numberOfComponents; i++)
	{
		MyBaseModel MBM(verts_multi_component[i], faces_multi_component[i]);
		resModels.push_back(MBM);
	}
	return resModels;
}
