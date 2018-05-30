#ifndef WrapDoubleCylinder_H
#define WrapDoubleCylinder_H

/*
* WrapDoubleCylinder.hpp
*
* Obstacle Set Algorithm Simulation for Double Cylinder Obstacles
*
*/

#include "WrapObst.h"

class WrapDoubleCylinder : public WrapObst
{
private:
	Eigen::Vector3d
		vec_z_U,      // U Cylinder Positive z axis
		point_U,      // U Cylinder Origin
		vec_z_V,      // V Cylinder Positive z axis
		point_V,      // V Cylinder Origin
		point_g,
		point_h;

	Eigen::MatrixXd
		M_U,          // Obstacle Coord Transformation Matrix for U
		M_V;          // Obstacle Coord Transformation Matrix for V

	double
		radius_U,     // U Cylinder Radius
		radius_V;     // V Cylinder Radius

public:
	// default constructor
	WrapDoubleCylinder()
	{
		vec_z_U = point_U = vec_z_V = point_V = point_g = point_h =
			Eigen::Vector3d(0.0, 0.0, 0.0);
		type = double_cylinder;
	}

	void setCylinderConfig(const Eigen::Vector3d &U,
		const Eigen::Vector3d &Z_U,
		const Eigen::Vector3d &V,
		const Eigen::Vector3d &Z_V)
	{
		this->point_U = U;
		this->vec_z_U = Z_U;
		this->point_V = V;
		this->vec_z_V = Z_V;
	}

	// constructor
	WrapDoubleCylinder(const Eigen::Vector3d &P,
		const Eigen::Vector3d &S,
		const Eigen::Vector3d &U,
		const Eigen::Vector3d &Z_U,
		const double R_U,
		const Eigen::Vector3d &V,
		const Eigen::Vector3d &Z_V,
		const double R_V)
		: WrapObst(P, S, Eigen::Vector3d(0.0, 0.0, 0.0), 0.0),
		point_U(U), vec_z_U(Z_U), radius_U(R_U),
		point_V(V), vec_z_V(Z_V), radius_V(R_V)
	{
		type = double_cylinder;
	}

	using WrapObst::compute;
	void compute();

	using WrapObst::getPoints;
	Eigen::MatrixXd getPoints(int num_points);
};

#endif