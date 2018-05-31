#pragma once

#ifndef MUSCLEMASS_SRC_RIGID_H_
#define MUSCLEMASS_SRC_RIGID_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Shape;
class Program;
class MatrixStack;
class Joint;
class WrapSphere;
class WrapCylinder;
class Particle;
class WrapDoubleCylinder;

typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 3, 6> Matrix3x6d;
typedef Eigen::Matrix<double, 6, 3> Matrix6x3d;
typedef Eigen::Matrix<double, 5, 6> Matrix5x6d;

class Rigid
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Rigid();
	Rigid(const std::shared_ptr<Shape> shape, Eigen::Matrix3d _R, Eigen::Vector3d _p, Eigen::Vector3d _dimension, double _r, double _m, bool _isReduced);
	virtual ~Rigid();
	void tare();
	void reset();
	void step(double h);
	
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;
	void computeForces();
	void computeTempForces();

	// set

	void setIndex(int _i);
	void setP(Eigen::Vector3d p);
	void setR(Eigen::Matrix3d R);
	void setTwist(Vector6d _twist);
	void setForce(Vector6d _force);
	void setParent(std::shared_ptr<Rigid> _parent);
	void setCylinderStatus(bool _isCylinder);
	void setDoubleCylinderStatus(bool _isDoubleCylinder);
	void setSphereStatus(bool _isSphere);

	void addChild(std::shared_ptr<Rigid> _child);
	void setJointAngle(double _theta, bool isDrawing);
	void setRotationAngle(double _theta);
	void addSphere(std::shared_ptr<WrapSphere> _sphere);
	void addCylinder(std::shared_ptr<WrapCylinder> _cylinder);
	void addDoubleCylinder(std::shared_ptr<WrapDoubleCylinder> _double_cylinders);
	void addPoint(std::shared_ptr<Particle> _point);

	// get
	Eigen::Vector3d getP() const;
	Eigen::Matrix3d getR() const;
	Matrix6d getMassMatrix() const;
	Vector6d getTwist() const;
	Vector6d getForce() const;
	std::shared_ptr<Joint> getJoint() const;
	Eigen::Matrix4d getE() const;
	Eigen::Matrix4d getEtemp() const;
	std::shared_ptr<Rigid> getParent() const;
	int getIndex() const;
	Eigen::Vector3d getDimension() const;

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
	bool isReduced; 
	bool isCylinder;
	bool isDoubleCylinder;
	bool isSphere;
	
	
	Eigen::Vector3d grav;	

private:
	const std::shared_ptr<Shape> box;
	std::shared_ptr<Rigid> parent;
	std::vector< std::shared_ptr<Rigid> > children;
	std::vector< std::shared_ptr<WrapSphere> > spheres;
	std::vector< std::shared_ptr<WrapCylinder> > cylinders;
	std::vector< std::shared_ptr<WrapDoubleCylinder> > double_cylinders;

	std::vector< std::shared_ptr<Particle> > points;
	std::shared_ptr<Joint> joint;
	Eigen::Matrix4d E_W_0; // Where current transform is wrt world
	Eigen::Matrix4d E_W_0_0; // Where current transform is wrt world at the start
	Eigen::Matrix4d E_W_0_temp; // for computational purpose, not for drawing
	Eigen::Vector3d dimension;
	Matrix6d mass_mat;
	Vector6d twist;
	Vector6d force;
	int i;  // index
};
#endif // MUSCLEMASS_SRC_RIGID_H_