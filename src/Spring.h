#pragma once
#ifndef __Spring__
#define __Spring__

#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Program;
class MatrixStack;
class Particle;

class Spring
{
public:
	Spring(std::shared_ptr<Particle> p0, std::shared_ptr<Particle> p1, double _mass);
	virtual ~Spring();
	double computeLength();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;
	//Eigen::MatrixXd computeJacobianMatrix();
	void setPosBeforePert();

	std::shared_ptr<Particle> p0;
	std::shared_ptr<Particle> p1;
	double E;
	double L;
	double mass;

	Eigen::Vector3d p0_b;
	Eigen::Vector3d p1_b;
	
};

#endif
