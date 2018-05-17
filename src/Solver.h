#pragma once

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Rigid;

class Solver
{
public:
	Solver(std::vector< std::shared_ptr<Rigid> > _boxes, bool _isReduced);
	virtual ~Solver();
	void step(double h);

private:
	std::vector< std::shared_ptr<Rigid> > boxes;
	Eigen::MatrixXd A;
	Eigen::MatrixXd M;
	Eigen::MatrixXd J;
	Eigen::VectorXd x;
	Eigen::VectorXd b;
	Eigen::VectorXd f;
	bool isReduced;
};

