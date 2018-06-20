#pragma once

#ifndef MUSCLEMASS_SRC_VECTOR_H_
#define MUSCLEMASS_SRC_VECTOR_H_

#include <vector>
#include <memory>
#include "Particle.h"

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Particle;
class Program;
class MatrixStack;

class Vector
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	Vector();
	~Vector();
	void reset();
	void setP(std::shared_ptr<Particle> _p) { this->p = _p; }
	void draw(std::shared_ptr<MatrixStack> MV, std::shared_ptr<MatrixStack> P, const std::shared_ptr<Program> p) const;
	void update(Eigen::Matrix4d E);

	Eigen::Vector3d dir0; // initial direction
	Eigen::Vector3d dir;  // current direction

	bool fixed;

private:
	std::shared_ptr<Particle> p;	// starting point
};

#endif // MUSCLEMASS_SRC_VECTOR_H_