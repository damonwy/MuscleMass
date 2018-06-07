#pragma once

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Rigid;
class Spring;
class Particle;
enum Integrator { RKF45, SYMPLECTIC};

class Solver
{
public:
	Solver(std::vector< std::shared_ptr<Rigid> > _boxes, std::vector< std::shared_ptr<Spring> > _springs, bool _isReduced, Integrator _time_integrator);
	virtual ~Solver();
	void step(double h);
	void dynamics(double t, double y[], double yp[]);
	Eigen::MatrixXd getJ_twist_thetadot();
	std::vector<double> solve(double y[], double yp[], const int neqn, double t_s, double t_end);
	Eigen::VectorXd solve_once(double y[], double yp[], const int neqn, double t_s, double t_end, int n_step);
	Integrator time_integrator;

	float r4_abs(float x);
	float r4_epsilon();
	void r4_fehl(void f(float t, float y[], float yp[]), int neqn,
		float y[], float t, float h, float yp[], float f1[], float f2[], float f3[],
		float f4[], float f5[], float s[]);
	float r4_max(float x, float y);
	float r4_min(float x, float y);
	int r4_rkf45(void f(float t, float y[], float yp[]), int neqn,
		float y[], float yp[], float *t, float tout, float *relerr, float abserr,
		int flag);
	float r4_sign(float x);
	double r8_abs(double x);
	double r8_epsilon();
	void r8_fehl(void (Solver::*f)(double t, double y[], double yp[]), int neqn,
		double y[], double t, double h, double yp[], double f1[], double f2[], double f3[],
		double f4[], double f5[], double s[]);
	double r8_max(double x, double y);
	double r8_min(double x, double y);
	int r8_rkf45(void (Solver::*f)(double t, double y[], double yp[]), int neqn,
		double y[], double yp[], double *t, double tout, double *relerr, double abserr,
		int flag);
	double r8_sign(double x);
	void timestamp();
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
	
	
};

