#include "MLError.h"
#include "JointSpaceInfo.h"
#include "time.h"
#include <iostream>
#include <iomanip>
# include <cstdlib>
# include <iomanip>

using namespace std;

JointSpaceInfo::JointSpaceInfo() {

}


FunctionTestJointSpaceInfoJointSpaceInfo::FunctionTestJointSpaceInfoJointSpaceInfo(double val) : JointSpaceInfo()
{
	this->val = val;
}


MLError FunctionTestJointSpaceInfoJointSpaceInfo::computeDifference(JointSpaceInfo *other, double *result)
{
	FunctionTestJointSpaceInfoJointSpaceInfo* otherTest = dynamic_cast<FunctionTestJointSpaceInfoJointSpaceInfo*>(other);
	if (!otherTest)
	{
		return MLError("error in dynamic casting");
	}
	*result = abs(val - otherTest->val);

	return MLError();
}


MLError FunctionTestJointSpaceInfoJointSpaceInfo::computeErrorVector(JointSpaceInfo *other, Eigen::VectorXd *result)
{
	FunctionTestJointSpaceInfoJointSpaceInfo* otherTest = dynamic_cast<FunctionTestJointSpaceInfoJointSpaceInfo*>(other);
	if (!otherTest)
	{
		return MLError("error in dynamic casting");
	}
	*result = abs(val - otherTest->val) * Eigen::VectorXd::Ones(1);

	return MLError();
}



MLError FunctionTestJointSpaceInfoJointSpaceInfo::add(double weight, JointSpaceInfo *other)
{
	FunctionTestJointSpaceInfoJointSpaceInfo* otherTest = dynamic_cast<FunctionTestJointSpaceInfoJointSpaceInfo*>(other);
	if (!otherTest)
	{
		return MLError("error in dynamic casting");
	}
	val += weight*otherTest->val;

	return MLError();
}

MLError FunctionTestJointSpaceInfoJointSpaceInfo::getHomeophicShapeInfo(JointSpaceInfo *source, JointSpaceInfo **result)
{
	FunctionTestJointSpaceInfoJointSpaceInfo* sourceTest = dynamic_cast<FunctionTestJointSpaceInfoJointSpaceInfo*>(source);
	if (!sourceTest)
	{
		return MLError("error in dynamic casting");
	}
	*result = new FunctionTestJointSpaceInfoJointSpaceInfo(this->val);

	return MLError();
}

MLError FunctionTestJointSpaceInfoJointSpaceInfo::addBoundaryConditions(bool useForceAngle)
{
	return MLError();
}


void FunctionTestJointSpaceInfoJointSpaceInfo::addNewPrecomputedPhysics(std::string name, std::vector<double> &values)
{
}


void FunctionTestJointSpaceInfoJointSpaceInfo::log()
{
	cout << "val = " << val << endl;
}
