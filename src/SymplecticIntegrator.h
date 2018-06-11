#pragma once

#ifndef MUSCLEMASS_SRC_SYMPLECTICINTEGRATOR_H_
#define MUSCLEMASS_SRC_SYMPLECTICINTEGRATOR_H_
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Rigid;
class Spring;
class Particle;
class Joint;

class SymplecticIntegrator {
public:
	SymplecticIntegrator(std::vector< std::shared_ptr<Rigid> > _boxes, std::vector< std::shared_ptr<Joint>> _joints, std::vector< std::shared_ptr<Spring> > _springs, bool _isReduced);
	void step(double h);
	Eigen::MatrixXd getJ_twist_thetadot();
	virtual ~SymplecticIntegrator();

	double m;
	double n;
	const int num_joints;

private:
	std::vector< std::shared_ptr<Rigid> > boxes;
	std::vector< std::shared_ptr<Spring> > springs;
	std::vector< std::shared_ptr<Joint> > joints;
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

#endif // MUSLEMASS_SRC_SYMPLECTICINTEGRATOR

