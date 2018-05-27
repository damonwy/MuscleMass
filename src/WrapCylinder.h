#ifndef WrapCylinder_H
#define WrapCylinder_H

/*
* WrapCylinder.hpp
*
* Obstacle Set Algorithm Simulation for Cylinder Obstacles
*
*/
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Shape;
class Program;
class MatrixStack;
class Rigid;

#include "WrapObst.h"

class WrapCylinder : public WrapObst
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
	Eigen::Vector3d vec_z;      // Cylinder Positive z axis
	Eigen::Matrix4d E_W_0;		// Where current transform is wrt world
	std::shared_ptr<Rigid> parent;
	const std::shared_ptr<Shape> cylinder_shape;

public:
	// default constructor
	WrapCylinder()
	{
		vec_z = Eigen::Vector3d(0.0, 0.0, 0.0);
		type = cylinder;
	}

	void setCylinderConfig(const Eigen::Vector3d &O,
		const Eigen::Vector3d &Z)
	{
		this->point_O = O;
		this->vec_z = Z;
	}

	// constructor
	WrapCylinder(const std::shared_ptr<Shape> s, const Eigen::Vector3d &P,
		const Eigen::Vector3d &S,
		const Eigen::Vector3d &O,
		const Eigen::Vector3d &Z,
		const double R)
		: WrapObst(P, S, O, R), vec_z(Z), cylinder_shape(s)
	{
		type = cylinder;
	}

	using WrapObst::compute;
	void compute();

	using WrapObst::getPoints;
	Eigen::MatrixXd getPoints(int num_points, double &theta_s, double &theta_e, Eigen::Matrix3d &_M);

	void reset();
	void step(double h);
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog) const;

	// get
	Eigen::Vector3d getP() const;
	Eigen::Matrix3d getR() const;

	double r; // radius

};

#endif