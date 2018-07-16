#pragma once
#ifndef MUSCLEMASS_SRC_MLWRAP_H_
#define MUSCLEMASS_SRC_MLWRAP_H_

#include "MLCommon.h"
#include "MLObject.h"

class MLBody;
class MLError;
class MLPoint;

class MLWrap :public MLObject {
public:
	explicit MLWrap();
	virtual ~MLWrap();

	virtual MLError load(const const nlohmann::json &elem, MLWorld *world, const std::string &folder);
	virtual void save(std::ofstream &ofs);

	virtual void init();
	virtual void draw();

	const std::vector<MLPoint *> * getSamples() const { return &m_samples; }

protected:
	std::vector<MLPoint *> m_samples;	// sample points along the wrapping
};

#endif // MUSCLEMASS_SRC_MLWRAP_H_