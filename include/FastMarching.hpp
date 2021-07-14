#pragma once
#include "ComputeDistanceField.hpp"
#include <queue> //priority_queue
#include <limits>
#include <cfloat>
#include <Eigen\dense>
#include <math.h>
//#include "MyMathSupport.hpp"

using namespace std;

class FastMarching :
	public ComputeDistanceField
{
protected:
	vector<int> m_parents;
public:
	FastMarching(const MyHalfEdgeModel& model, const set<int>&sources);
	void Run();
	double GetFastMarchingDistance(int p1, int p2, int p3, double p1_dist, double p2_dist);
	void FastMarchingUnfoldTriangleto2D(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3,
		double p1_dist, double p2_dist, Eigen::Vector2d &q1, Eigen::Vector2d &q2, Eigen::Vector2d &q3);
	vector<Eigen::Vector3d> BacktracePath(int destination) const;
};

FastMarching::FastMarching(const MyHalfEdgeModel& model, const set<int>&sources)
	: ComputeDistanceField(model, sources)
{
}

void FastMarching::Run()
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
		for (int i = 0; i < neighbors.size(); i++)
		{
			auto neigh = neighbors[i];
			if (m_distanceField[neigh] > FLT_MAX) {
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
			else {
				auto preFMPoint = neighbors[(i - 1 + neighbors.size()) % neighbors.size()];
				auto nxtFMPoint = neighbors[(i + 1) % neighbors.size()];
				double newPreDis = GetFastMarchingDistance(neigh, top.v, preFMPoint, m_distanceField[neigh], top.dis);
				if (newPreDis < m_distanceField[preFMPoint])
				{
					Event evt;
					evt.v = preFMPoint;
					evt.parent = top.v;
					evt.ancestor = top.ancestor;
					evt.dis = newPreDis;
					pending.push(evt);
				}
				double newNxtDis = GetFastMarchingDistance(neigh, top.v, nxtFMPoint, m_distanceField[neigh], top.dis);
				if (newNxtDis < m_distanceField[nxtFMPoint])
				{
					Event evt;
					evt.v = nxtFMPoint;
					evt.parent = top.v;
					evt.ancestor = top.ancestor;
					evt.dis = newNxtDis;
					pending.push(evt);
				}
			}
		}
	}
}

double FastMarching::GetFastMarchingDistance(int p1, int p2, int p3, double p1_dist, double p2_dist) {
	Eigen::Vector2d q1, q2, q3;
	FastMarchingUnfoldTriangleto2D(m_model.GetVertices()[p1], m_model.GetVertices()[p2], m_model.GetVertices()[p3],
		p1_dist, p2_dist, q1, q2, q3);
	Eigen::Vector2d leftPoint, rightPoint;
	if (q1.x() < q2.x()) {
		leftPoint = q1;
		rightPoint = q2;
	}
	else if (q1.x() == q2.x()) {
		leftPoint = q1;
		rightPoint = q1;
	}
	else {
		leftPoint = q2;
		rightPoint = q1;
	}
	double dist;
	if (q3.x() > rightPoint.x()) {
		dist = rightPoint.y() + (rightPoint - q3).norm();
	}
	else if (q3.x() < leftPoint.x()) {
		dist = leftPoint.y() + (leftPoint - q3).norm();
	}
	else {
		dist = q3.y();
	}
	return max(max(q1.y(), q2.y()), dist);
}

void FastMarching::FastMarchingUnfoldTriangleto2D(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3,
	double p1_dist, double p2_dist, Eigen::Vector2d &q1, Eigen::Vector2d &q2, Eigen::Vector2d &q3)
{
	double a = (p3 - p2).norm();
	double b = (p3 - p1).norm();
	double c = (p1 - p2).norm();
	q1 = Eigen::Vector2d(0, p1_dist);
	q2 = Eigen::Vector2d(sqrt(c * c - (p2_dist - p1_dist) * (p2_dist - p1_dist)), p2_dist);
	double angle = acos((c * c + a * a - b * b) / 2 * c * a) - asin((p2_dist - p1_dist) / c);
	q3 = Eigen::Vector2d(q2.x() - a * cos(angle), q2.y() + a * sin(angle));
}

vector<Eigen::Vector3d> FastMarching::BacktracePath(int destination) const
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