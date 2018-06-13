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
	void setParent(std::shared_ptr<Rigid> _parent) { this->parent = _parent; }
	std::shared_ptr<Rigid> getParent() const { return this->parent; }
	Eigen::Vector3d getTempPos() const { return this->x_temp; }
	
	double r; // radius
	double m; // mass
	int i;  // starting index
	Eigen::Vector3d x0; // initial position
	Eigen::Vector3d v0; // initial velocity
	Eigen::Vector3d x;  // position
	Eigen::Vector3d v;  // velocity
	Eigen::Vector3d x_temp; // temporary position, for computation 
	bool fixed;
	Eigen::Vector3d normal;
	std::shared_ptr<Rigid> parent;
	
private:
	const std::shared_ptr<Shape> sphere;
};

#endif // MUSCLEMASS_SRC_PARTICLE_H_
