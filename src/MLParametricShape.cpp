#include "MLParametricShape.h"

#include "MLShapeInfo.h"
#include "MLError.h"
#include "MLDerivInfo.h"

MLParametricShape::MLParametricShape() {

}

MLError MLParametricShape::getMidPoint(std::vector<double> *result) {
	for (int i = 0; i < getNParams(); i++) {
		double minVal = getMinRange(i);
		double maxVal = getMaxRange(i);
		if (maxVal < minVal){
			return(MLError("ranges are invalid"));
		}
		result->push_back(0.5 * (minVal + maxVal));
	}
	return MLError();
}

MLError MLParametricShape::validateParameters(const std::vector<double> &parameters) const {
	if (parameters.size() != getNParams()) {
		return MLError("number of parameters is incorrect");
	}

	for (size_t i = 0; i < parameters.size(); i++) {
		double minVal = getMinRange(i);
		double maxVal = getMaxRange(i);
		if ((parameters[i] < minVal) || (parameters[i] > maxVal)) {
			return (MLError("ranges are invalid"));
		}
	}
	return MLError();
}

MLError MLParametricShape::getRelativeParamVal(int i, double alpha, double *result) {
	if ((i >= getNParams()) || (alpha < 0) || (alpha > 1)){
		return MLError("parameter out of range");
	} 

	*result = getMinRange(i) * (1 - alpha) + getMaxRange(i) * alpha;
	
	return MLError();
}

std::vector<double> MLParametricShape::mapToStandardHypercube(const std::vector<double> & params) {
	std::vector<double> result(params.size());
	for (int i = 0; i < params.size(); i++) {
		double alpha = (params[i] - getMinRange(i)) / (getMaxRange(i) - getMinRange(i));
		result[i] = alpha * 2 - 1;
	}

	return result;
}

std::vector<double> MLParametricShape::mapFromStandardHypercube(const std::vector<double> &unitParams) {
	std::vector<double> result(unitParams.size());
	for (int i = 0; i < unitParams.size(); i++) {
		double alpha = (unitParams[i] + 1.0) / 2.0;
		result[i] = getMinRange(i) * (1 - alpha) + getMaxRange(i) * alpha;
	}

	return result;
}

void MLParametricShape::getRanges(std::vector<double> * ranges) {
	for (int i = 0; i < getNParams(); i++) {
		ranges->push_back(getMinRange(i));
		ranges->push_back(getMaxRange(i));
	}
}
