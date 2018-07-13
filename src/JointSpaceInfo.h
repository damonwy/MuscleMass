#pragma once
#ifndef MUSCLEMASS_SRC_JOINTSPACEINFO_H_
#define MUSCLEMASS_SRC_JOINTSPACEINFO_H_

#include <Eigen\Dense>
#include <vector>
#include <string>

class MLError;

class JointSpaceInfo {
public:
	JointSpaceInfo();
	virtual ~JointSpaceInfo() {};
	virtual MLError computeDifference(JointSpaceInfo *other, double *result) = 0;
	virtual MLError computeErrorVector(JointSpaceInfo *other, Eigen::VectorXd *result) = 0;
	virtual MLError add(double weight, JointSpaceInfo *other) = 0;
	virtual MLError getHomeophicShapeInfo(JointSpaceInfo *source, JointSpaceInfo **result) = 0;
	virtual MLError addBoundaryConditions(bool useForceAngle) = 0;
	virtual void addNewPrecomputedPhysics(std::string name, std::vector<double> &values) = 0;
	virtual void clearAllPhysics() = 0;
	virtual void log() = 0;
};


class FunctionTestJointSpaceInfoJointSpaceInfo: public JointSpaceInfo {
public:
	FunctionTestJointSpaceInfoJointSpaceInfo(double val);
	
	MLError computeDifference(JointSpaceInfo *other, double *result);
	MLError computeErrorVector(JointSpaceInfo *other, Eigen::VectorXd *result);
	MLError add(double weight, JointSpaceInfo *other);
	MLError getHomeophicShapeInfo(JointSpaceInfo *source, JointSpaceInfo **result);
	MLError addBoundaryConditions(bool useForceAngle);
	void addNewPrecomputedPhysics(std::string name, std::vector<double> &values);
	void clearAllPhysics() {}
	void log();

	double val;
};


#endif // MUSCLEMASS_SRC_ML_ERROR_