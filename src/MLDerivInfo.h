#pragma once

#ifndef MUSCLEMASS_SRC_MLDERIVINFO_H_
#define MUSCLEMASS_SRC_MLDERIVINFO_H_

class MLParametricShape;
class MLError;
class MLShapeInfo;
class MLFunctionTestShapeInfo;


#include <Eigen\Dense>
#include <vector>

class MLDerivInfo {
public:
	MLDerivInfo();

	static MLError newFromWeightedSum(const std::vector<std::pair<double, MLShapeInfo*>> weightedSamples, MLDerivInfo **result);
	virtual MLError add(double weight, MLShapeInfo *other) = 0;
	virtual void log() = 0;
};

class MLFunctionTestDerivInfo : public MLDerivInfo {
public:
	MLFunctionTestDerivInfo(double val);
	MLFunctionTestDerivInfo(MLFunctionTestDerivInfo* other, double weight);
	MLError add(double weight, MLShapeInfo *other);
	void log();

	double val;
};




#endif // MUSCLEMASS_SRC_MLDERIVINFO_H_