#include "MLDerivInfo.h"

#include <iostream>

#include "MLParametricShape.h"
#include "MLError.h"
#include "MLShapeInfo.h"

MLDerivInfo::MLDerivInfo() {
}

MLError MLDerivInfo::newFromWeightedSum(const std::vector<std::pair<double, MLShapeInfo*>> weightedSamples, MLDerivInfo **result) {
	MLError err;

	if (weightedSamples.size() == 0){
		return MLError("cannot create new from empty list of weighted samples");
	}

	if (dynamic_cast<MLFunctionTestShapeInfo*>(weightedSamples[0].second)) {
		MLFunctionTestDerivInfo* derivInfo = new MLFunctionTestDerivInfo(0);
		for (auto weightedSample : weightedSamples) {
			MLErrorReturn(derivInfo->add(weightedSample.first, weightedSample.second));
		}
		*result = derivInfo;
		return err;
	}

	return MLError("Shape Info type not specified");
}

//--------------------------------------------------------------------------------

MLFunctionTestDerivInfo::MLFunctionTestDerivInfo(double val) : MLDerivInfo() {
	this->val = val;
}

MLError MLFunctionTestDerivInfo::add(double weight, MLShapeInfo *other) {
	MLFunctionTestDerivInfo* otherTest = dynamic_cast<MLFunctionTestDerivInfo*>(other);
	if (!otherTest) {
		return MLError("error in dynamic casting");
	}
	val += weight*otherTest->val;
	
	return MLError();
}

void MLFunctionTestDerivInfo::log() {
	std::cout << "val = " << val << std::endl;
}

//--------------------------------------------------------------------------------

