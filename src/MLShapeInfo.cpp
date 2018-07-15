#include "MLShapeInfo.h"

#include "time.h"
#include <iostream>
# include <cstdlib>
# include <iomanip>

#include "MLError.h"

using namespace std;

MLShapeInfo::MLShapeInfo() {

}


MLFunctionTestShapeInfo::MLFunctionTestShapeInfo(double val) : MLShapeInfo()
{
	this->val = val;
}


MLError MLFunctionTestShapeInfo::computeDifference(MLShapeInfo *other, double *result)
{
	MLFunctionTestShapeInfo* otherTest = dynamic_cast<MLFunctionTestShapeInfo*>(other);
	if (!otherTest)
	{
		return MLError("error in dynamic casting");
	}
	*result = abs(val - otherTest->val);

	return MLError();
}


MLError MLFunctionTestShapeInfo::computeErrorVector(MLShapeInfo *other, Eigen::VectorXd *result)
{
	MLFunctionTestShapeInfo* otherTest = dynamic_cast<MLFunctionTestShapeInfo*>(other);
	if (!otherTest)
	{
		return MLError("error in dynamic casting");
	}
	*result = abs(val - otherTest->val) * Eigen::VectorXd::Ones(1);

	return MLError();
}



MLError MLFunctionTestShapeInfo::add(double weight, MLShapeInfo *other)
{
	MLFunctionTestShapeInfo* otherTest = dynamic_cast<MLFunctionTestShapeInfo*>(other);
	if (!otherTest)
	{
		return MLError("error in dynamic casting");
	}
	val += weight*otherTest->val;

	return MLError();
}

MLError MLFunctionTestShapeInfo::getHomeophicMLShapeInfo(MLShapeInfo *source, MLShapeInfo **result)
{
	MLFunctionTestShapeInfo* sourceTest = dynamic_cast<MLFunctionTestShapeInfo*>(source);
	if (!sourceTest)
	{
		return MLError("error in dynamic casting");
	}
	*result = new MLFunctionTestShapeInfo(this->val);

	return MLError();
}

MLError MLFunctionTestShapeInfo::addBoundaryConditions(bool useForceAngle)
{
	return MLError();
}


void MLFunctionTestShapeInfo::addNewPrecomputedPhysics(std::string name, std::vector<double> &values)
{
}


void MLFunctionTestShapeInfo::log()
{
	cout << "val = " << val << endl;
}

//---------------------------------------------------

MLJointSpaceShapeInfo::MLJointSpaceShapeInfo() {
}

MLJointSpaceShapeInfo::MLJointSpaceShapeInfo(const nlohmann::json &elem) {


}

MLJointSpaceShapeInfo::~MLJointSpaceShapeInfo() {
	clear();
}


MLJointSpaceShapeInfo::MLJointSpaceShapeInfo(MLJointSpaceShapeInfo* other, double weight) {
	MLJointSpaceShapeInfo* otherJS = dynamic_cast<MLJointSpaceShapeInfo*>(other);
	m_js_vals = weight * otherJS->m_js_vals;
}


MLError MLJointSpaceShapeInfo::computeDifference(MLShapeInfo *other, double *result) {
	MLJointSpaceShapeInfo* otherJS = dynamic_cast<MLJointSpaceShapeInfo*>(other);
	if (!otherJS) {
		return MLError("error in dynamic casting");
	}
	if (m_js_vals.size() != otherJS->m_js_vals.size()) {
		return MLError("cannot add joint spaces with different dimemsions");
	}
	
	double diff = (m_js_vals - otherJS->m_js_vals).norm();
	*result = diff;
	return MLError();
}

MLError MLJointSpaceShapeInfo::computeErrorVector(MLShapeInfo *other, Eigen::VectorXd *result) {
	// NO NEED?
}


MLError MLJointSpaceShapeInfo::add(double weight, MLShapeInfo *other) {
	MLJointSpaceShapeInfo* otherJS = dynamic_cast<MLJointSpaceShapeInfo*>(other);
	if (!otherJS) {
		return MLError("error in dynamic casting");
	}
	if (m_js_vals.size() != otherJS->m_js_vals.size()) {
		return MLError("cannot add joint spaces with different dimemsions");
	}
	m_js_vals += weight * otherJS->m_js_vals;
	return MLError();
}


void MLJointSpaceShapeInfo::log() {
	cout << "vals = " << endl;
	for (int i = 0; i < (int)m_js_vals.size(); ++i) {
		cout << m_js_vals(i) << endl;
	}
}

void MLJointSpaceShapeInfo::clear() {
	
}
