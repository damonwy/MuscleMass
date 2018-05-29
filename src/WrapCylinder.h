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
class Particle;
class Vector;
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
	Eigen::Matrix4d E_P_0;		// Where the local frame is wrt parent
	std::shared_ptr<Rigid> parent;
	std::shared_ptr<Particle> P;
	std::shared_ptr<Particle> S;
	std::shared_ptr<Particle> O;
	std::shared_ptr<Vector> Z;
	const std::shared_ptr<Shape> cylinder_shape;
	int num_points;

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
	WrapCylinder(const std::shared_ptr<Shape> s, 
		const Eigen::Vector3d &p,
		const Eigen::Matrix3d &R,
		const double _r,
		const int _num_points)
		: WrapObst(), cylinder_shape(s)
	{		
		//WrapObst(_P->x, _S->x, _O->x, _r),	

		this->num_points = _num_points;
		E_P_0.setIdentity();
		this->E_P_0.block<3, 3>(0, 0) = R;
		this->E_P_0.block<3, 1>(0, 3) = p;
		this->E_W_0 = E_P_0;
		
		this->type = cylinder;
		this->r = _r;
	}

	using WrapObst::compute;
	void compute();

	using WrapObst::getPoints;
	Eigen::MatrixXd getPoints(int num_points, double &theta_s, double &theta_e, Eigen::Matrix3d &_M) const;

	void reset();
	void step(double h);
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;

	// get
	Eigen::Vector3d getP() const;
	Eigen::Matrix3d getR() const;
	Eigen::Matrix4d getE() const;
	Eigen::Matrix4d getE_P_0() const;

	// set
	void setp(Eigen::Vector3d p);
	void setR(Eigen::Matrix3d R);
	void setE(Eigen::Matrix4d E);
	void setP(std::shared_ptr<Particle> _P);
	void setS(std::shared_ptr<Particle> _S);
	void setO(std::shared_ptr<Particle> _O);
	void setZ(std::shared_ptr<Vector> _Z);
	void setNumPoints(int _num_points);

	double r; // radius

};

#endif