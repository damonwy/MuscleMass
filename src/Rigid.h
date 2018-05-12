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

	double r; // radius
	double m; // mass
	int i;  // starting index
	Eigen::Vector3d x0; // initial position
	Eigen::Vector3d v0; // initial velocity
	Eigen::Vector3d x;  // position
	Eigen::Vector3d v;  // velocity
	bool fixed;
	Eigen::Vector3d normal;

private:
	const std::shared_ptr<Shape> box;
};

#endif
