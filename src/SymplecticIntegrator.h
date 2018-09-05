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
class Program;
class MatrixStack;

class SymplecticIntegrator {
public:
	SymplecticIntegrator(std::vector< std::shared_ptr<Rigid> > _boxes, std::vector< std::shared_ptr<Joint>> _joints, std::vector< std::shared_ptr<Spring> > _springs, bool _isReduced, bool _isMuscle, int _num_samples, Eigen::Vector3d _grav, double _epsilon);
	void step(double h);
	void step_matlab(double h);
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P) const;
	Eigen::MatrixXd getJ();
	Eigen::MatrixXd getJdot(Eigen::VectorXd thetadotlist);
	Eigen::MatrixXd getGlobalJacobian(Eigen::VectorXd thetalist);
	virtual ~SymplecticIntegrator();

	double m;
	double n;
	const int num_joints;

private:
	std::vector< std::shared_ptr<Rigid> > boxes;
	std::vector < std::shared_ptr<Rigid> > moving_boxes;
	std::vector< std::shared_ptr<Spring> > springs;
	std::vector< std::shared_ptr<Joint> > joints;
	Eigen::MatrixXd A;
	Eigen::MatrixXd M;
	Eigen::MatrixXd J;
	Eigen::MatrixXd Jdot;
	Eigen::VectorXd x;
	Eigen::VectorXd b;
	Eigen::VectorXd f;


	double epsilon;
	bool isReduced;
	bool isMuscle;
	Eigen::Vector3d grav;
	std::vector< std::shared_ptr<Particle> > debug_points;
	int num_samples;	// The number of samples along the muscle lines
};

#endif // MUSLEMASS_SRC_SYMPLECTICINTEGRATOR_H_

