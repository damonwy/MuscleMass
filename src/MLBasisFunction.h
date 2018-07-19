#pragma once
#ifndef MUSCLEMASS_SRC_MLBASICFUNCTION_H_
#define MUSCLEMASS_SRC_MLBASICFUNCTION_H_

#include <map>
#include <vector>

#include "MLCommon.h"

class MLError;
class MLParametricShape;
struct MLBasisFunctionKey;
class MLAdaptiveGridCell;
class MLPrecomputedSample;

class MLBasisFunction {
public:
	enum MLBasisFunctionType { LINEAR_BSPLINE, CUBIC_BSPLINE };
	MLBasisFunction();

	virtual MLBasisFunctionType getType() = 0;
	virtual MLError eval(const Eigen::VectorXd &pos, double *result) = 0;
	virtual MLError evalDeriv(const Eigen::VectorXd &pos, int direction, double *result) = 0;
	virtual MLError refine(int dir, std::vector<MLBasisFunction*> *basisFunctions) = 0;
	virtual void log() = 0;

	Eigen::VectorXd& getSupport();
	Eigen::VectorXd& getCenter();

	double getWeight();
	void setWeightToZero();
	void setWeight(double weight);
	void addWeight(double val);
	bool isEqual(MLBasisFunction* other);
	bool checkRegionOverlap(const std::vector<double> &uncoveredRegion);
	MLError refineAllDirections(std::vector<MLBasisFunction*> *basisFunctions);
	MLError checkIsCoveredByCells(const std::vector<MLAdaptiveGridCell*> &adjacentCells, const MLParametricShape* paramShape);

	MLError computeReplacingWeights(const std::vector<std::vector<MLPrecomputedSample*>> &replacingSamples,
		int splitDirection, const Eigen::VectorXd &cellSize, int *replacingSamplesIndex, std::vector<double> *replacingWeights);

	static MLError newFromCopyWithWeights(MLBasisFunction* basisFunction, double weight, MLBasisFunction** result);
	static MLError newFromCopyKey(const MLBasisFunctionKey & basisFunctionKey, double weight, MLBasisFunction** result);

protected:
	Eigen::VectorXd m_center;
	Eigen::VectorXd m_support;
	double m_weight;
};

class MLLinearBSpline : public MLBasisFunction
{
public:
	MLLinearBSpline(const Eigen::VectorXd &center, const Eigen::VectorXd &support, double weigth);
	MLBasisFunctionType getType() { return LINEAR_BSPLINE; }
	MLError eval(const Eigen::VectorXd &pos, double *result);
	MLError evalDeriv(const Eigen::VectorXd &pos, int direction, double *result);
	MLError refine(int dir, std::vector<MLBasisFunction*> *basisFunctions);
	void log();

};

class MLCubicBSpline : public MLBasisFunction
{
public:
	MLCubicBSpline(const Eigen::VectorXd &center, const Eigen::VectorXd &support, double weigth);
	MLBasisFunctionType getType() { return CUBIC_BSPLINE; }
	MLError eval(const Eigen::VectorXd &pos, double *result);
	MLError evalDeriv(const Eigen::VectorXd &pos, int direction, double *result);
	MLError refine(int dir, std::vector<MLBasisFunction*> *basisFunctions);
	void log();

};


#endif MUSCLEMASS_SRC_MLBASICFUNCTION_H_
