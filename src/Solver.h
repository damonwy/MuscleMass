#pragma once

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Rigid;

class Solver
{
public:
	Solver(std::vector< std::shared_ptr<Rigid> > _boxes);
	virtual ~Solver();
	void step(double h);

private:
	std::vector< std::shared_ptr<Rigid> > boxes;
	Eigen::MatrixXd A;
	Eigen::VectorXd x;
	Eigen::VectorXd b;
};

