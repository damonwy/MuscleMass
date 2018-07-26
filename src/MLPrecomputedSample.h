#pragma once

#ifndef MUSCLEMASS_SRC_MLPRECOMPUTEDSAMPLE_H_
#define MUSCLEMASS_SRC_MLPRECOMPUTEDSAMPLE_H_

#include <Eigen\Dense>
#include <vector>
#include <unordered_map>

#include "MLBasisFunction.h"

class MLShapeInfo;
class MLError;
class MLBasisFunction;
class MLParametricShape;
// class MLAdaptiveGridCell;

class GroupFunctionInfo {
public:
	GroupFunctionInfo() {
		m_initialized = false;
		m_weightSum = 0;
	}

	void add(std::shared_ptr<MLPrecomputedSample> sample, double weight);
	Eigen::VectorXd getCenter() {
		return m_centerSum / m_weightSum;
	}

	double getWeight() {
		return m_weightSum;
	}

private:
	Eigen::VectorXd m_centerSum;
	double m_weightSum;
	bool m_initialized;
};

class MLPrecomputedSample {
public:
	MLPrecomputedSample(const Eigen::VectorXd center, std::shared_ptr<MLBasisFunction> basisFunction, std::shared_ptr<MLShapeInfo> shapeInfo);
	MLPrecomputedSample(const Eigen::VectorXd center, std::shared_ptr<MLShapeInfo> shapeInfo);
	MLPrecomputedSample(const Eigen::VectorXd center, std::vector<std::shared_ptr<MLBasisFunction>> &basisFunctions, std::shared_ptr<MLShapeInfo> shapeInfo);
	void removeBasisFunction(std::shared_ptr<MLBasisFunction> basisFunctions);

	MLError getInterpolationWeight(const Eigen::VectorXd pos, double * result);
	MLError getDerivInterpolationWeight(const Eigen::VectorXd pos, int direction, double * result);
	Eigen::VectorXd& getCenter();
	std::shared_ptr<MLShapeInfo> getShapeInfo();

	MLError getBasisType(MLBasisFunction::MLBasisFunctionType *type);
	void setShapeInfo(std::shared_ptr<MLShapeInfo> shapeInfo);

	void log();
private:
	Eigen::VectorXd m_center_;
	std::shared_ptr<MLShapeInfo> m_shapeInfo_;
	std::unordered_map<MLBasisFunctionKey, std::shared_ptr<MLBasisFunction>> m_basisFunctions;

};


#endif MUSCLEMASS_SRC_MLPRECOMPUTEDSAMPLE_H_