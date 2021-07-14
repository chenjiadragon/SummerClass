#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
struct FaceInfo2
{
	FaceInfo2() {}
	int nesting_level;
	bool in_domain() {
		return nesting_level % 2 == 1;
	}
};
typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>               TDS;
typedef CGAL::Exact_predicates_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
typedef CDT::Point                                                Point;
typedef CGAL::Polygon_2<K>                                        Polygon_2;
typedef CDT::Vertex_handle                                          Vertex_handle;
typedef CDT::Face_handle                                          Face_handle;
void
mark_domains(CDT& ct,
	Face_handle start,
	int index,
	std::list<CDT::Edge>& border)
{
	if (start->info().nesting_level != -1) {
		return;
	}
	std::list<Face_handle> queue;
	queue.push_back(start);
	while (!queue.empty()) {
		Face_handle fh = queue.front();
		queue.pop_front();
		if (fh->info().nesting_level == -1) {
			fh->info().nesting_level = index;
			for (int i = 0; i < 3; i++) {
				CDT::Edge e(fh, i);
				Face_handle n = fh->neighbor(i);
				if (n->info().nesting_level == -1) {
					if (ct.is_constrained(e)) border.push_back(e);
					else queue.push_back(n);
				}
			}
		}
	}
}
//explore set of facets connected with non constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
void
mark_domains(CDT& cdt)
{
	for (CDT::Face_handle f : cdt.all_face_handles()) {
		f->info().nesting_level = -1;
	}
	std::list<CDT::Edge> border;
	mark_domains(cdt, cdt.infinite_face(), 0, border);
	while (!border.empty()) {
		CDT::Edge e = border.front();
		border.pop_front();
		Face_handle n = e.first->neighbor(e.second);
		if (n->info().nesting_level == -1) {
			mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
		}
	}
}

#include "MyBaseModel.hpp"

MyBaseModel AddBase(const MyBaseModel& model, double thickness)
{
	double zMin = DBL_MAX;
	double xMin = DBL_MAX, xMax = -DBL_MAX;
	double yMin = DBL_MAX, yMax = -DBL_MAX;
	vector<double> scalarField;
	for (int i = 0; i < model.GetVertices().size(); ++i)
	{
		scalarField.push_back(model.GetVertices()[i].z());
		if (model.GetVertices()[i].z() < zMin)
		{
			zMin = model.GetVertices()[i].z();
		}

		if (model.GetVertices()[i].x() < xMin)
		{
			xMin = model.GetVertices()[i].x();
		}
		if (model.GetVertices()[i].x() > xMax)
		{
			xMax = model.GetVertices()[i].x();
		}
		if (model.GetVertices()[i].y() < yMin)
		{
			yMin = model.GetVertices()[i].y();
		}
		if (model.GetVertices()[i].y() > yMax)
		{
			yMax = model.GetVertices()[i].y();
		}
	}

	auto isoline = model.ExtractIsoline(scalarField, zMin + 0.5 * thickness);   //得到那一圈点
	MyBaseModel upperPart = model.SplitModelByIsoline(scalarField, zMin + 0.5 * thickness).second;

	xMin -= 2 * thickness;
	xMax += 2 * thickness;
	yMin -= 2 * thickness;
	yMax += 2 * thickness;

	Polygon_2 polygon_outer;
	polygon_outer.push_back(Point(xMin, yMin));
	polygon_outer.push_back(Point(xMax, yMin));
	polygon_outer.push_back(Point(xMax, yMax));
	polygon_outer.push_back(Point(xMin, yMax));

	CDT cdt;
	cdt.insert_constraint(polygon_outer.vertices_begin(), polygon_outer.vertices_end(), true);

	//2021/7/8 RA
	for (auto loop : isoline)
	{
		Polygon_2 polygon_in;
		for (auto p : loop) {
			polygon_in.push_back(Point(p.x(), p.y()));
		}
		cdt.insert_constraint(polygon_in.vertices_begin(), polygon_in.vertices_end(), true);
	}
	mark_domains(cdt);
	vector<Eigen::Vector3d> final_verts = upperPart.GetVertices();
	vector<Eigen::Vector3i> final_faces = upperPart.GetFaces();
	cout << final_verts.size() << endl;
	//base thicked
	Eigen::Vector3d v1, v2, v3, v4, v5, v6, v7, v8;
	v1 << xMin, yMin, zMin + 0.5 * thickness;
	v2 << xMax, yMin, zMin + 0.5 * thickness;
	v3 << xMax, yMax, zMin + 0.5 * thickness;
	v4 << xMin, yMax, zMin + 0.5 * thickness;
	v5 << xMin, yMin, zMin - 0.5 * thickness;
	v6 << xMax, yMin, zMin - 0.5 * thickness;
	v7 << xMax, yMax, zMin - 0.5 * thickness;
	v8 << xMin, yMax, zMin - 0.5 * thickness;
	final_verts.push_back(v1);
	final_verts.push_back(v2);
	final_verts.push_back(v3);
	final_verts.push_back(v4);
	final_verts.push_back(v5);
	final_verts.push_back(v6);
	final_verts.push_back(v7);
	final_verts.push_back(v8);
	int Pnum = final_verts.size() - 1;
	cout << Pnum << endl;
	Eigen::Vector3i f1, f2, f3, f4, f5, f6, f7, f8, f9, f10;
	f1 << Pnum - 3, Pnum - 2, Pnum - 6;
	f2 << Pnum - 3, Pnum - 6, Pnum - 7;
	f3 << Pnum - 3, Pnum, Pnum - 1;
	f4 << Pnum - 3, Pnum - 1, Pnum - 2;
	f5 << Pnum - 5, Pnum - 6, Pnum - 1;
	f6 << Pnum - 6, Pnum - 2, Pnum - 1;
	f7 << Pnum - 7, Pnum - 4, Pnum;
	f8 << Pnum - 7, Pnum, Pnum - 3;
	f9 << Pnum - 4, Pnum - 5, Pnum;
	f10 << Pnum, Pnum - 5, Pnum - 1;
	cout << f10 << endl;
	final_faces.push_back(f1);
	final_faces.push_back(f2);
	final_faces.push_back(f3);
	final_faces.push_back(f4);
	final_faces.push_back(f5);
	final_faces.push_back(f6);
	final_faces.push_back(f7);
	final_faces.push_back(f8);
	final_faces.push_back(f9);
	final_faces.push_back(f10);
	cout << Pnum << endl;
	MyBaseModel* final_model;


	//points connected 
	for (Face_handle f : cdt.finite_face_handles())
	{
		if (f->info().in_domain())
		{
			Eigen::Vector3i f_verts;
			for (int i = 0; i < 3; i++)
			{
				Vertex_handle p = f->vertex(i);
				auto P_inbase = p->point();
				Eigen::Vector3d pt(P_inbase.x(), P_inbase.y(), zMin + 0.5 * thickness);
				double minDis = DBL_MAX;
				int Pid = -1;
				for (int j = 0; j < final_verts.size(); j++)
				{
					auto v = final_verts[j];
					double dis = sqrt((v.x() - pt.x()) * (v.x() - pt.x()) + (v.y() - pt.y()) * (v.y() - pt.y()) + (v.z() - pt.z()) * (v.z() - pt.z()));
					if (dis < minDis)
					{
						minDis = dis;
						Pid = j;
					}
				}
				cout << Pid << endl;
				if (i == 0)
					f_verts.x() = Pid;
				if (i == 1)
					f_verts.y() = Pid;
				if (i == 2)
					f_verts.z() = Pid;

			}
			final_faces.push_back(f_verts);
		}
	}
	final_model = new MyBaseModel(final_verts, final_faces);
	return *final_model;
	//mark_domains(cdt);

	//return upperPart;
}