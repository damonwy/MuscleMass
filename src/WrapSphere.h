#pragma once
#ifndef MUSCLEMASS_SRC_WRAPSPHERE_H_
#define MUSCLEMASS_SRC_WRAPSPHERE_H_

/*
* WrapSphere.hpp
*
* Obstacle Set Algorithm Simulation for Sphere Obstacles
*
*/

#include "WrapObst.h"

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


class WrapSphere : public WrapObst
{
private:
	Eigen::MatrixXd arc_points; // each col stores the position of a point in that arc

	std::shared_ptr<Rigid> parent;
	std::shared_ptr<Particle> P;
	std::shared_ptr<Particle> S;
	std::shared_ptr<Particle> O;	// the origin
	
	int num_points;

public:
	// default constructor
	WrapSphere()
	{
		type = sphere;
	}

	void setSphereConfig(const Eigen::Vector3d &O)
	{
		this->point_O = O;
	}

	// constructor
	WrapSphere(const double R, const int _num_points)
		: WrapObst(), num_points(_num_points)
	{
		this->radius = R;
		type = sphere;
		arc_points.resize(3, this->num_points + 1);
	}

	using WrapObst::compute;
	void compute();

	using WrapObst::getPoints;
	Eigen::MatrixXd getPoints(int num_points);

	void reset();
	void step();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;

	void setP(std::shared_ptr<Particle> _P);
	void setS(std::shared_ptr<Particle> _S);
	void setO(std::shared_ptr<Particle> _O);
	void setParent(std::shared_ptr<Rigid> _parent);

};

#endif //MUSCLEMASS_SRC_WRAPSPHERE_H_