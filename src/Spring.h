#pragma once
#ifndef MUSCLEMASS_SRC_SPRING_H_
#define MUSCLEMASS_SRC_SPRING_H_
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "TypeDef.h"

class Program;
class MatrixStack;
class Particle;
class Rigid;

class Spring
{
public:
	Spring(std::shared_ptr<Particle> p0, std::shared_ptr<Particle> p1, double _mass, int num_samples, Eigen::Vector3d _grav, double _epsilon);
	virtual ~Spring();
	
	void step();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;
	void setPosBeforePert(); // Save the postions of p0 and p1 before perturbation
	Eigen::Vector3d getP0BeforePert() const { return this->p0_b; }
	Eigen::Vector3d getP1BeforePert() const { return this->p1_b; }
	
	void setSamples(std::vector < std::shared_ptr<Particle> > _samples) { this->samples = _samples; }
	void updateSamples();
	
	Vector12d getBoxTwists() const { return this->phi_box; }
	Eigen::Vector2d getBoxID() const { return this->box_id; };
	double getPotentialEnergy() const { return this->V; }
	double getKineticEnergy() const { return this->K; }
	std::vector<std::shared_ptr<Particle> > getSamples() const { return this->samples; }

	double computeLength();
	void computeEnergy();
	static Eigen::MatrixXd computeMassMatrix(std::vector<std::shared_ptr<Spring> > springs, int num_boxes);
	static Eigen::VectorXd computeGravity(std::vector<std::shared_ptr<Spring> > springs, int num_boxes);

	Eigen::Vector2d box_id;

	std::shared_ptr<Particle> p0;
	std::shared_ptr<Particle> p1;
	double E;	// stiffness
	double L;	// initial length
	double l;	// current length
	double mass;
	double epsilon;

	Eigen::Vector3d p0_b; // the position of p0 before perturbation
	Eigen::Vector3d p1_b; // the position of p1 before perturbation

	std::vector<std::shared_ptr<Particle> > samples;	// sample points along the spring
	Eigen::Vector3d grav;

	double V;
	double K;
	Vector12d phi_box;	
};

#endif // MUSCLEMASS_SRC_SPRING_H_
