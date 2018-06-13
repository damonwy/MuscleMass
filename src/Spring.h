#pragma once
#ifndef MUSCLEMASS_SRC_SPRING_H_
#define MUSCLEMASS_SRC_SPRING_H_

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
	void setPosBeforePert(); // Save the postions of p0 and p1 before perturbation
	Eigen::Vector3d getP0BeforePert() const { return this->p0_b; }
	Eigen::Vector3d getP1BeforePert() const { return this->p1_b; }

	std::shared_ptr<Particle> p0;
	std::shared_ptr<Particle> p1;
	double E;
	double L;
	double mass;

	Eigen::Vector3d p0_b; // the position of p0 before perturbation
	Eigen::Vector3d p1_b; // the position of p1 before perturbation
	
};

#endif // MUSCLEMASS_SRC_SPRING_H_
