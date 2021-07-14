#pragma once
#include "ComputeDistanceField.hpp"
#include <queue> //priority_queue
#include <limits>
#include <cfloat>
using namespace std;

class Dijkstra :
	public ComputeDistanceField
{
protected:
	vector<int> m_parents;
public:
	Dijkstra(const MyHalfEdgeModel& model, const set<int>&sources);
	void Run();
	vector<Eigen::Vector3d> BacktracePath(int destination) const;
};

Dijkstra::Dijkstra(const MyHalfEdgeModel& model, const set<int>&sources)
	: ComputeDistanceField(model, sources)
{
}

void Dijkstra::Run()
{
	m_distanceField.clear();
	m_distanceField.resize(m_model.GetVertices().size(), DBL_MAX);
	m_ancestors.clear();
	m_ancestors.resize(m_model.GetVertices().size(), -1);
	m_parents.clear();
	m_parents.resize(m_model.GetVertices().size(), -1);

	struct Event
	{
		int v;
		double dis;
		int parent;
		int ancestor;
		bool operator>(const Event& other) const
		{
			return dis > other.dis;
		}
	};
	priority_queue<Event, vector<Event>, greater<Event> > pending;
	for (auto source : m_sources)
	{
		Event evt;
		evt.v = source;
		evt.dis = 0;
		evt.parent = -1;
		evt.ancestor = source;
		pending.push(evt);
	}

	while (!pending.empty())
	{
		auto top = pending.top();
		pending.pop();
		if (m_distanceField[top.v] < FLT_MAX) // updated already
		{
			continue;
		}
		m_distanceField[top.v] = top.dis;
		m_ancestors[top.v] = top.ancestor;
		m_parents[top.v] = top.parent;

		vector<int> neighbors = m_model.GetNeighboringVertices(top.v);
		for (auto neigh : neighbors)
		{
			auto newDis = top.dis + (m_model.GetVertices()[top.v] - m_model.GetVertices()[neigh]).norm();
			if (newDis < m_distanceField[neigh])
			{
				Event evt;
				evt.v = neigh;
				evt.parent = top.v;
				evt.ancestor = top.ancestor;
				evt.dis = newDis;
				pending.push(evt);
			}
		}
	}
}

vector<Eigen::Vector3d> Dijkstra::BacktracePath(int destination) const
{
	vector<Eigen::Vector3d> path;
	int nxt = destination;
	while (nxt != -1)
	{
		path.push_back(m_model.GetVertices()[nxt]);
		nxt = m_parents[nxt];
	}
	reverse(path.begin(), path.end());
	return path;
}