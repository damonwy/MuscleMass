#pragma once
#ifndef MUSCLEMASS_SRC_MLSHAPEINFO_H_
#define MUSCLEMASS_SRC_MLSHAPEINFO_H_

#include <Eigen\Dense>
#include <vector>
#include <string>

#include "MLCommon.h"

class MLError;

class MLShapeInfo {
public:
	MLShapeInfo();
	virtual ~MLShapeInfo() {};
	virtual MLError computeDifference(MLShapeInfo *other, double *result) = 0;
	virtual MLError computeErrorVector(MLShapeInfo *other, Eigen::VectorXd *result) = 0;
	virtual MLError add(double weight, MLShapeInfo *other) = 0;
	virtual MLError getHomeophicMLShapeInfo(MLShapeInfo *source, MLShapeInfo **result) = 0;
	virtual MLError addBoundaryConditions(bool useForceAngle) = 0;
	virtual void addNewPrecomputedPhysics(std::string name, std::vector<double> &values) = 0;
	virtual void clearAllPhysics() = 0;
	virtual void log() = 0;
};


class MLFunctionTestShapeInfo: public MLShapeInfo {
public:
	MLFunctionTestShapeInfo(double val);
	
	MLError computeDifference(MLShapeInfo *other, double *result);
	MLError computeErrorVector(MLShapeInfo *other, Eigen::VectorXd *result);
	MLError add(double weight, MLShapeInfo *other);
	MLError getHomeophicMLShapeInfo(MLShapeInfo *source, MLShapeInfo **result);
	MLError addBoundaryConditions(bool useForceAngle);
	void addNewPrecomputedPhysics(std::string name, std::vector<double> &values);
	void clearAllPhysics() {}
	void log();

	double val;
};


class MLJointSpaceShapeInfo : public MLShapeInfo {
public:
	MLJointSpaceShapeInfo();
	MLJointSpaceShapeInfo(const nlohmann::json &elem);
	~MLJointSpaceShapeInfo();
	MLJointSpaceShapeInfo(MLJointSpaceShapeInfo* other, double weight);
	MLError computeDifference(MLShapeInfo *other, double *result);
	MLError computeErrorVector(MLShapeInfo *other, Eigen::VectorXd *result);
	MLError add(double weight, MLShapeInfo *other);
	void log();
	void clear();

	Eigen::VectorXd m_js_vals; // inertial matrix, etc..
};


#endif // MUSCLEMASS_SRC_MLSHAPEINFO_H_