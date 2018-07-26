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
	this->m_val = val;
}

MLError MLFunctionTestShapeInfo::computeDifference(shared_ptr<MLShapeInfo> other, double *result)
{
	auto otherTest = std::dynamic_pointer_cast<MLFunctionTestShapeInfo>(other);
	if (!otherTest)
	{
		return MLError("error in dynamic casting");
	}
	*result = abs(m_val - otherTest->m_val);

	return MLError();
}


MLError MLFunctionTestShapeInfo::computeErrorVector(shared_ptr<MLShapeInfo> other, Eigen::VectorXd *result)
{
	auto otherTest = std::dynamic_pointer_cast<MLFunctionTestShapeInfo>(other);
	if (!otherTest)
	{
		return MLError("error in dynamic casting");
	}
	*result = abs(m_val - otherTest->m_val) * Eigen::VectorXd::Ones(1);

	return MLError();
}

MLError MLFunctionTestShapeInfo::add(double weight, std::shared_ptr<MLShapeInfo> other)
{
	auto otherTest = std::dynamic_pointer_cast<MLFunctionTestShapeInfo>(other);
	
	if (!otherTest)
	{
		return MLError("error in dynamic casting");
	}
	m_val += weight*otherTest->m_val;

	return MLError();
}

MLError MLFunctionTestShapeInfo::getHomeophicMLShapeInfo(std::shared_ptr<MLShapeInfo> source, MLShapeInfo **result)
{
	auto sourceTest = std::dynamic_pointer_cast<MLFunctionTestShapeInfo>(source);
	if (!sourceTest)
	{
		return MLError("error in dynamic casting");
	}
	*result = new MLFunctionTestShapeInfo(this->m_val);

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
	cout << "val = " << m_val << endl;
}

//---------------------------------------------------

MLJointSpaceShapeInfo::MLJointSpaceShapeInfo() {
}

MLJointSpaceShapeInfo::MLJointSpaceShapeInfo(const nlohmann::json &elem) {


}

MLJointSpaceShapeInfo::~MLJointSpaceShapeInfo() {
	clear();
}


MLJointSpaceShapeInfo::MLJointSpaceShapeInfo(std::shared_ptr<MLShapeInfo> other, double weight) {

	auto otherJS = std::dynamic_pointer_cast<MLJointSpaceShapeInfo>(other);
	m_js_vals = weight * otherJS->m_js_vals;
}


MLError MLJointSpaceShapeInfo::computeDifference(std::shared_ptr<MLShapeInfo> other, double *result) {
	auto otherJS = std::dynamic_pointer_cast<MLJointSpaceShapeInfo>(other);
	
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

MLError MLJointSpaceShapeInfo::computeErrorVector(std::shared_ptr<MLShapeInfo> other, Eigen::VectorXd *result) {
	// NO NEED?
	return MLError();
}


MLError MLJointSpaceShapeInfo::add(double weight, std::shared_ptr<MLShapeInfo> other) {

	auto otherJS = std::dynamic_pointer_cast<MLJointSpaceShapeInfo>(other);
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
