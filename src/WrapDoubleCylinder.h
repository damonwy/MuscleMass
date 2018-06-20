#pragma once
#ifndef MUSCLEMASS_SRC_WRAPDOUBLECYLINDER_H_
#define MUSCLEMASS_SRC_WRAPDOUBLECYLINDER_H_

/*
* WrapDoubleCylinder.hpp
*
* Obstacle Set Algorithm Simulation for Double Cylinder Obstacles
*
*/
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

#include "WrapObst.h"
#include "Particle.h"
#include "Vector.h"

class Particle;
class Vector;
class Shape;
class Program;
class MatrixStack;
class Rigid;

class WrapDoubleCylinder : public WrapObst
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
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

	Status
		status_U,     // U Wrapping Status
		status_V;     // V Wrapping Status

	std::shared_ptr<Particle> U;	// U Cylinder Origin
	Eigen::Matrix4d E_W_U;			// Where current transform is wrt world(updated)
	Eigen::Matrix4d E_P_U;			// Where the local frame is wrt parent(fixed)

	std::shared_ptr<Rigid> parent_U;
	std::shared_ptr<Vector> z_U;		// Z axis of Cylinder U(direction matters)

	std::shared_ptr<Particle> V;		// V Cylinder Origin
	Eigen::Matrix4d E_W_V;		
	Eigen::Matrix4d E_P_V;		
 
	std::shared_ptr<Rigid> parent_V;
	std::shared_ptr<Vector> z_V;		// Z axis of Cylinder V(direction matters)

	std::shared_ptr<Particle> P;
	std::shared_ptr<Particle> S;

	Eigen::MatrixXd arc_points;	// each col stores the position of a point in that arc
	const std::shared_ptr<Shape> cylinder_shape;
	int num_points;

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
	WrapDoubleCylinder(const std::shared_ptr<Shape> s, 	
		const double R_U,
		const double R_V,
		const int _num_points)
		: WrapObst(),
		radius_U(R_U), radius_V(R_V), cylinder_shape(s), num_points(_num_points)
	{
		type = double_cylinder;	
		this->arc_points.resize(3, 3 * this->num_points + 1);
	}

	using WrapObst::compute;
	void compute();

	using WrapObst::getPoints;
	Eigen::MatrixXd getPoints(int num_points);

	Status get_status_u() const { return status_U; }	
	Status get_status_v() const { return status_V; }

	void step();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;

	// get
	Eigen::Vector3d getp_U() const { return this->E_W_U.block<3, 1>(0, 3); }
	Eigen::Matrix3d getR_U() const { return this->E_W_U.block<3, 3>(0, 0); }
	Eigen::Matrix4d getE_U() const { return this->E_W_U; }
	Eigen::Matrix4d getE_P_U() const { return this->E_P_U; }

	Eigen::Vector3d getp_V() const { return this->E_W_V.block<3, 1>(0, 3); }
	Eigen::Matrix3d getR_V() const { return this->E_W_V.block<3, 3>(0, 0); }

	Eigen::Matrix4d getE_V() const { return this->E_W_V; }
	Eigen::Matrix4d getE_P_V() const { return this->E_P_V; }

	std::shared_ptr<Rigid> getParent_U() const { return this->parent_U; }
	std::shared_ptr<Rigid> getParent_V() const { return this->parent_V; }

	// set
	void setE_U(Eigen::Matrix4d E) { this->E_W_U = E; }
	void setE_V(Eigen::Matrix4d E) { this->E_W_V = E; }

	void setP(std::shared_ptr<Particle> _P) { this->P = _P; this->point_P = P->x; }
	void setS(std::shared_ptr<Particle> _S) { this->S = _S; this->point_S = S->x; }

	void setU(std::shared_ptr<Particle> _U) { this->U = _U; this->point_U = U->x; }
	void setV(std::shared_ptr<Particle> _V) { this->V = _V; this->point_V = V->x; }

	void setZ_U(std::shared_ptr<Vector> _z_U) { this->z_U = _z_U; this->vec_z_U = this->z_U->dir; }
	void setZ_V(std::shared_ptr<Vector> _z_V) { this->z_V = _z_V; this->vec_z_V = this->z_V->dir; }
	void setNumPoints(int _num_points) { this->num_points = _num_points; }

	void setParent_U(std::shared_ptr<Rigid> _parent_U) { this->parent_U = _parent_U; }
	void setParent_V(std::shared_ptr<Rigid> _parent_V) { this->parent_V = _parent_V; }

	void setE_P_U(Eigen::Matrix4d E) { this->E_P_U = E; }
	void setE_P_V(Eigen::Matrix4d E) { this->E_P_V = E; }	
};

#endif // MUSCLEMASS_SRC_WRAPDOUBLECYLINDER_H_