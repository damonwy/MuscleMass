#pragma once
#ifndef MUSCLEMASS_SRC_WRAPOBST_H_
#define MUSCLEMASS_SRC_WRAPOBST_H_

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <complex>

enum Status { wrap, inside_radius, no_wrap, empty };
enum Type { none, sphere, cylinder, double_cylinder };
#define PI 3.141593

class WrapObst
{
protected:
	Eigen::Vector3d
		point_P,    // Bounding-Fixed Via Point 1
		point_S,    // Bounding-Fixed Via Point 2
		point_O,    // Obstacle Center Point
		point_q,    // Obstacle Via Point 1 in Obstacle Frame
		point_t;    // Obstacle Via Point 2 in Obstacle Frame

	Eigen::MatrixXd M;  // Obstacle Coord Transformation Matrix
	Status status;      // Wrapping Status
	Type type;          // Obstacle Type
	double path_length,  // Wrapping Path Length
		radius;       // obstacle sphere radius

public:
	// set muscle origin point
	void setOrigin(const Eigen::Vector3d &P)
	{
		this->point_P = P;
	}

	// set muscle insertion point
	void setInsertion(const Eigen::Vector3d &S)
	{
		this->point_S = S;
	}

	// default constructor
	WrapObst()
	{
		point_P = point_S = point_O = point_q = point_t
			= Eigen::Vector3d(0.0, 0.0, 0.0);
		M = Eigen::MatrixXd(3, 3);
		status = empty;
		path_length = 0.0;
		radius = 0.0;
		type = none;
	}

	// constructor
	WrapObst(const Eigen::Vector3d &P,
		const Eigen::Vector3d &S,
		const Eigen::Vector3d &O,
		const double R) :
		point_P(P), point_S(S), point_O(O), radius(R)
	{
		point_q = point_t = Eigen::Vector3d(0.0, 0.0, 0.0);
		M = Eigen::MatrixXd(3, 3);
		status = empty;
		path_length = 0.0;
		type = none;
	}

	// wrap calculation
	void compute() {}

	double getLength()
	{
		return this->path_length;
	}

	Status getStatus() const
	{
		return this->status;
	}

	double getRadius()
	{
		return this->radius;
	}

	Eigen::MatrixXd getPoints() {}

};

#endif // MUSCLEMASS_SRC_WRAPOBST_H_