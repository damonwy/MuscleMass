#pragma once
#ifndef MUSCLEMASS_SRC_MLPARAMETRICSHAPE_H_
#define MUSCLEMASS_SRC_MLPARAMETRICSHAPE_H_

#include <vector>
#include <Eigen/Dense>

class MLShapeInfo;
class MLError;
class MLDerivInfo;

// This is an abstract class that defines a parametric shape
// A parametric shape is abstracticly defined by a feasible set and a function that maps point
// in this feasible set to shapes. In our code, the feasible set is always a hypercube defined by a set
// of ranges. The functions bellow return information about the feasible set (validate parameters, help sample
// the feasible set, etc) and also return the shape (and the derivative of the shape) for a given
// parameters configuration.

class MLParametricShape {
public:
	MLParametricShape();

	virtual int getNParams() const = 0;

	virtual double getMinRange(int iParam) const = 0;

	virtual double getMaxRange(int iParam) const = 0;

	virtual MLError evalShapeInfo(const std::vector<double> &params, MLShapeInfo **result) = 0;

	virtual MLError evalDeriv(const std::vector<double> &params, int direction, MLDerivInfo **result) = 0;

	// returns the parameters at the center of the parameter space
	MLError getMidPoint(std::vector<double> *result);
	// checks if parameters are valid
	MLError validateParameters(const std::vector<double> &parameters) const;
	// maps a value between 0 and 1 to the ranges of param i 
	MLError getRelativeParamVal(int i, double alpha, double *result);
	// maps a point in the hypercube with ranges [-1 1] to the parameter space (given by the ranges)
	std::vector<double> MLParametricShape::mapToStandardHypercube(const std::vector<double> &params);
	// maps a point in parameter space to the hypercube with ranges [-1 1]
	std::vector<double> MLParametricShape::mapFromStandardHypercube(const std::vector<double> &unitParams);
	// returns the vector of ranges 
	void getRanges(std::vector<double> * ranges);

};


#endif // MUSCLEMASS_SRC_MLPARAMETRICSHAPE_H_