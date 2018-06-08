#pragma once

#ifndef MUSCLEMASS_SRC_RKF45INTEGRATOR_H_
#define MUSCLEMASS_SRC_RKF45INTEGRATOR_H_
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Rigid;
class Spring;
class Particle;

class RKF45Integrator {
public:
	RKF45Integrator(std::vector< std::shared_ptr<Rigid> > _boxes, std::vector< std::shared_ptr<Spring> > _springs, bool _isReduced);

	double m;
	double n;
	const int num_joints;

private:
	std::vector< std::shared_ptr<Rigid> > boxes;
	std::vector< std::shared_ptr<Spring> > springs;

	Eigen::MatrixXd A;
	Eigen::MatrixXd M;
	Eigen::MatrixXd J;
	Eigen::VectorXd x;
	Eigen::VectorXd b;
	Eigen::VectorXd f;
	double epsilon;
	bool isReduced;
	Eigen::Vector3d grav;

};

#endif // MUSLEMASS_SRC_RKF45INTEGRATOR
