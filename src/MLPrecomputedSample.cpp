#include "MLPrecomputedSample.h"

#include <iostream>
#include "MLError.h"
#include "MLBasisFunction.h"
#include "MLParametricShape.h"

void GroupFunctionInfo::add(std::shared_ptr<MLPrecomputedSample> sample, double weight) {
	if (!m_initialized) {
		m_centerSum = sample->getCenter() * weight;
		m_initialized = true;
	}
	else {
		m_centerSum += sample->getCenter() * weight;
	}
	m_weightSum += weight;
}

MLPrecomputedSample::MLPrecomputedSample(const Eigen::VectorXd center, 
	std::shared_ptr<MLBasisFunction> basisFunction, 
	std::shared_ptr<MLShapeInfo> shapeInfo) {
	m_center_ = center;
	addBasisFunction(basisFunction);
	m_shapeInfo_ = shapeInfo;

}

void MLPrecomputedSample::addBasisFunction(std::shared_ptr<MLBasisFunction> basisFunction) {

}

std::shared_ptr<MLShapeInfo> MLPrecomputedSample::getShapeInfo() {
	return m_shapeInfo_;
}

void MLPrecomputedSample::setShapeInfo(std::shared_ptr<MLShapeInfo> shapeInfo) {
	m_shapeInfo_ = shapeInfo;
}

void MLPrecomputedSample::log() {
	std::cout << "Precomputed sample:" << std::endl;
	std::cout << "center : " << m_center_.transpose() << std::endl;
	for (auto basisFunction : m_basisFunctions_)
	{
		basisFunction.second->log();
	}
}