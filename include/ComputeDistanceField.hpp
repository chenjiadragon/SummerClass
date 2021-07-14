#pragma once
#include <vector>
#include <set>
#include "MyHalfEdgeModel.hpp"
using namespace std;
class ComputeDistanceField
{
protected:
	vector<double> m_distanceField;
	set<int> m_sources;
	vector<int> m_ancestors;
	const MyHalfEdgeModel& m_model;
public:
	ComputeDistanceField(const MyHalfEdgeModel& model, const set<int>&sources);
	virtual vector<Eigen::Vector3d> BacktracePath(int destination) const = 0;
	int GetAncestor(int destination) const;
	virtual void Run() = 0;
	vector<double> GetDistanceField() const { return m_distanceField; }
};

ComputeDistanceField::ComputeDistanceField(const MyHalfEdgeModel& model, const set<int>&sources)
	: m_model(model), m_sources(sources)
{
}



int ComputeDistanceField::GetAncestor(int destination) const
{
	return m_ancestors[destination];
}