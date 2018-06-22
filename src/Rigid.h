#pragma once

#ifndef MUSCLEMASS_SRC_RIGID_H_
#define MUSCLEMASS_SRC_RIGID_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "TypeDef.h"
#include "Joint.h"

class Shape;
class Program;
class MatrixStack;
class Joint;
class WrapSphere;
class WrapCylinder;
class Particle;
class WrapDoubleCylinder;

class Rigid
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Rigid();
	Rigid(const std::shared_ptr<Shape> shape, Eigen::Matrix3d _R, Eigen::Vector3d _p, Eigen::Vector3d _dimension, double _r, double _m, bool _isReduced, Eigen::Vector3d _grav);
	virtual ~Rigid();
	void tare();
	void reset();
	void step(double h);
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;
	void computeForces();	  // Use current E
	void computeTempForces(); // Use temporary E

	// set
	void setIndex(int _i){ this->i = _i; }
	void setP(Eigen::Vector3d p) { this->E_W_0.block<3, 1>(0, 3) = p; }
	void setR(Eigen::Matrix3d R) { this->E_W_0.block<3, 3>(0, 0) = R; }

	void setTwist(Vector6d _twist) { this->twist = _twist; }
	void setForce(Vector6d _force) { this->force = _force; }
	void setParent(std::shared_ptr<Rigid> _parent) { this->parent = _parent; }

	void setCylinderStatus(bool _isCylinder) { this->isCylinder = _isCylinder; }
	void setDoubleCylinderStatus(bool _isDoubleCylinder) { this->isDoubleCylinder = _isDoubleCylinder; }
	void setSphereStatus(bool _isSphere) { this->isSphere = _isSphere; }

	void setEtemp(Eigen::Matrix4d E) { this->E_W_0_temp = E; }
	void setE(Eigen::Matrix4d E) { this->E_W_0 = E; }
	void setJoint(std::shared_ptr<Joint> _joint) { this->joint = _joint; }
	void setRotationAngle(double _dtheta) { this->joint->setDTheta(_dtheta); }
	void setThetadot(double _thetadot) { this->joint->setThetadot(_thetadot); }

	void addPoint(std::shared_ptr<Particle> _point) { this->points.push_back(_point); }
	void addChild(std::shared_ptr<Rigid> _child) { this->children.push_back(_child); }
	void addSphere(std::shared_ptr<WrapSphere> _sphere) { this->spheres.push_back(_sphere); }
	void addCylinder(std::shared_ptr<WrapCylinder> _cylinder) { this->cylinders.push_back(_cylinder); }
	void addDoubleCylinder(std::shared_ptr<WrapDoubleCylinder> _double_cylinders) { this->double_cylinders.push_back(_double_cylinders); } 
	
	void setJointAngle(double _theta, bool isDrawing);
	void setSingleJointAngle(double _theta);
	
	void updateCylinders();
	void updateDoubleCylinders();
	void updateSpheres();
	void updatePoints();
	void updateTempPoints();

	// get
	Eigen::Matrix4d getE() const { return this->E_W_0; }
	Eigen::Vector3d getP() const { return this->E_W_0.block<3, 1>(0, 3); }
	Eigen::Matrix3d getR() const { return this->E_W_0.block<3, 3>(0, 0); }
	Eigen::Matrix4d getEtemp() const { return this->E_W_0_temp; }
	Eigen::Vector3d getDimension() const { return this->dimension; }
	double getAngle() const { return this->joint->getTheta(); }
	double getThetadot() const { return this->joint->getThetadot(); }
	double getPotentialEnergy() const { return this->V; }
	double getKineticEnergy() const { return this->K; }

	Vector6d getTwist() const { return this->twist; }
	Vector6d getForce() const { return this->force; }
	Matrix6d getMassMatrix() const { return this->mass_mat; }
	
	std::shared_ptr<Rigid> getParent() const { return this->parent; }
	std::shared_ptr<Joint> getJoint() const { return this->joint; }
	std::vector<std::shared_ptr<Particle>> getPoints() const { return this->points; }
	int getIndex() const{ return this->i; }

	static Eigen::Matrix4d inverse(const Eigen::Matrix4d &E);
	static Matrix3x6d gamma(const Eigen::Vector3d &r);
	static Matrix6d adjoint(const Eigen::Matrix4d &E);
	static Eigen::Matrix3d bracket3(const Eigen::Vector3d &a);
	static Eigen::Matrix4d bracket6(const Vector6d &a);
	static Eigen::Vector3d unbracket3(const Eigen::Matrix3d &A);
	static Vector6d unbracket6(const Eigen::Matrix4d &A);
	static Eigen::Matrix4d integrate(const Eigen::Matrix4d &E0, const Eigen::VectorXd &phi, double h);

	double r; // radius
	double m; // mass
	bool isReduced;		// Use reduced coord or maximal coord
	bool isCylinder;	// drawing?
	bool isDoubleCylinder;
	bool isSphere;
	
	Eigen::Vector3d grav;	

private:
	const std::shared_ptr<Shape> box;

	std::shared_ptr<Rigid> parent;
	std::shared_ptr<Joint> joint;

	std::vector< std::shared_ptr<Particle> > points;
	std::vector< std::shared_ptr<Rigid> > children;
	std::vector< std::shared_ptr<WrapSphere> > spheres;
	std::vector< std::shared_ptr<WrapCylinder> > cylinders;
	std::vector< std::shared_ptr<WrapDoubleCylinder> > double_cylinders;

	Eigen::Matrix4d E_W_0;		// Where current transform is wrt world
	Eigen::Matrix4d E_W_0_0;	// Where current transform is wrt world at the start
	Eigen::Matrix4d E_W_0_temp; // for computational purpose, not for drawing
	Eigen::Vector3d dimension;
	Matrix6d mass_mat;
	Vector6d twist;
	Vector6d force;

	double V;
	double K;

	int i;  // index
};

#endif // MUSCLEMASS_SRC_RIGID_H_