#include "MLBasisFunction.h"

#include <iostream>

#include "MLParametricShape.h"
#include "MLError.h"

#define EPSILON 0.000001

MLBasisFunction::MLBasisFunction()
{

}

//double MLBasisFunction::getWeight() {
//	return m_weight_;
//}
//
//void MLBasisFunction::setWeightToZero() {
//	m_weight_ = 0.0;
//}
//
//void MLBasisFunction::setWeight(double weight) {
//	m_weight_ = weight;
//}
//
//void MLBasisFunction::addWeight(double val) {
//	m_weight_ += val;
//}
//
//bool MLBasisFunction::isEqual(std::shared_ptr<MLBasisFunction> other) {
//	if (getType() != other->getType()) {
//		return false;
//	}
//
//	if (getCenter() != other->getCenter()) {
//		return false;
//	}
//
//	if (getSupport() != other->getSupport()) {
//		return false;
//	}
//	return true;
//}
//
//bool MLBasisFunction::checkRegionOverlap(const std::vector<double> &uncoveredRegion) {
//	int nParams = m_center_.size();
//	for (int i = 0; i < nParams; i++)
//	{
//		double minBasisFunction = m_center_[i] - m_support_[i];
//		double maxBasisFunction = m_center_[i] + m_support_[i];
//		double minRegion = uncoveredRegion[2 * i];
//		double maxRegion = uncoveredRegion[2 * i + 1];
//		if ((maxRegion - minRegion) < EPSILON)
//		{
//			//return false;
//		}
//		double allowedError = m_support_[i] / 100.0;
//		if ((maxBasisFunction < (minRegion + allowedError)) || (minBasisFunction >(maxRegion - allowedError)))
//		{
//			return false;
//		}
//	}
//	return true;
//}