#pragma once
#ifndef MUSCLEMASS_SRC_PARTICLE_H_
#define MUSCLEMASS_SRC_PARTICLE_H_
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Shape;
class Program;
class MatrixStack;
class Rigid;

class Particle
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Particle();
	Particle(const std::shared_ptr<Shape> shape);
	virtual ~Particle();
	void tare();
	void reset();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void update(Eigen::Matrix4d E);
	void updateTemp(Eigen::Matrix4d E);

	std::shared_ptr<Rigid> getParent() const { return this->parent; }
	Eigen::Vector3d getTempPos() const { return this->x_temp; }
	Eigen::MatrixXd getJacobianMatrix() const { return this->J; }
	double getPotentialEnergy() const { return this->V; }
	double getKineticEnergy() const { return this->K; }

	void setPotentialEnergy(double _V) { this->V = _V; }
	void setKineticEnergy(double _K) { this->K = _K; }
	void setJacobianMatrix(Eigen::MatrixXd _J) { this->J = _J; }
	void setJacobianMatrixCol(Eigen::Vector3d p, int icol) { this->J.col(icol) = p; }
	void setParent(std::shared_ptr<Rigid> _parent) { this->parent = _parent; }

	void clearJacobianMatrix() { this->J.setZero(); }
	
	double computePotentialEnergy(Eigen::Vector3d grav);
	double computeKineticEnergy(Eigen::VectorXd phi);

	bool fixed;
	double r;					// radius
	double m;					// mass
	int i;						// starting index
	Eigen::Vector3d x0;			// initial position
	Eigen::Vector3d v0;			// initial velocity
	Eigen::Vector3d x;			// position
	Eigen::Vector3d v;			// velocity
	Eigen::Vector3d x_temp;		// temporary position, for computation 
	Eigen::Vector3d normal;	
	double s;					// non-dimensional material coordinate [0,1]
								// 0 at muscle origin and 1 at insertion; remain fixed.
	
private:
	const std::shared_ptr<Shape> sphere;
	double V;					// potential energy
	double K;					// kinetic energy
	Eigen::MatrixXd J;			// Jacobian Matrix
	std::shared_ptr<Rigid> parent;
};

#endif // MUSCLEMASS_SRC_PARTICLE_H_
