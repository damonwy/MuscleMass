#pragma once
#ifndef __Rigid__
#define __Rigid__

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Shape;
class Program;
class MatrixStack;

typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 3, 6> Matrix3x6d;
typedef Eigen::Matrix<double, 6, 3> Matrix6x3d;

class Rigid
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	Rigid();
	Rigid(const std::shared_ptr<Shape> shape);
	virtual ~Rigid();
	void tare();
	void reset();
	void step(double h);
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;

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
	int i;  // index
	Vector6d twist;

	bool fixed;	

private:
	const std::shared_ptr<Shape> box;
};

#endif
