#pragma once
#include <Eigen\dense>
namespace MyMathSupport
{
	Eigen::Vector3d GetNormalVector(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3);
	void UnfoldTriangleFrom3Dto2D(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3,
		Eigen::Vector2d &q1, Eigen::Vector2d &q2, Eigen::Vector2d &q3);
	Eigen::Vector3d GetBarycentricCoord(Eigen::Vector2d q1, Eigen::Vector2d q2, Eigen::Vector2d q3, Eigen::Vector2d query);
	Eigen::Vector3d ComposePoint3DFromBarycentricCoord(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3, Eigen::Vector3d coord);
};

Eigen::Vector3d MyMathSupport::GetNormalVector(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3)
{
	auto vec1 = p2 - p1;
	auto vec2 = p3 - p2;
	return vec1.cross(vec2);
}

void MyMathSupport::UnfoldTriangleFrom3Dto2D(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3,
	Eigen::Vector2d &q1, Eigen::Vector2d &q2, Eigen::Vector2d &q3)
{
	auto dir1 = p2 - p1;
	dir1.normalize();
	auto dir2 = GetNormalVector(p1, p2, p3);
	dir2.normalize();
	auto dir3 = dir1.cross(dir2);
	dir3.normalize();

	//q1 = Eigen::Vector2d(dir3.dot(p1), dir1.dot(p1));      //Ì¯Æ½Çó×ø±ê
	//q2 = Eigen::Vector2d(dir3.dot(p2), dir1.dot(p2));
	//q3 = Eigen::Vector2d(dir3.dot(p3), dir1.dot(p3));

	q1 = Eigen::Vector2d(0, 0);
	q2 = Eigen::Vector2d(dir3.dot(p2- p1), dir1.dot(p2 - p1));
	q3 = Eigen::Vector2d(dir3.dot(p3 - p1), dir1.dot(p3 - p1));
}

Eigen::Vector3d MyMathSupport::GetBarycentricCoord
(Eigen::Vector2d q1, Eigen::Vector2d q2, Eigen::Vector2d q3, Eigen::Vector2d query)
{
	Eigen::Matrix3d m;
	m << q1.x(), q2.x(), q3.x(),
		q1.y(), q2.y(), q3.y(),
		1, 1, 1;
	Eigen::Vector3d b;
	b << query.x(), query.y(), 1;

	return m.inverse() * b;
}

Eigen::Vector3d MyMathSupport::ComposePoint3DFromBarycentricCoord(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3, Eigen::Vector3d coord)
{
	return coord.x() * p1 + coord.y() * p2 + coord.z() * p3;
}