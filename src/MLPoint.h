#pragma once
#ifndef MUSCLEMASS_SRC_MLPOINT_H_
#define MUSCLEMASS_SRC_MLPOINT_H_

#include "MLCommon.h"
#include "MLObject.h"

class MLPoint : public MLObject{
public:
	explicit MLPoint();
	virtual ~MLPoint();

	virtual MLError load(const const nlohmann::json &elem, MLWorld *world, const std::string &folder);
	virtual void save(std::ofstream &ofs);

	virtual void init();
	virtual void draw();



};

#endif // MUSCLEMASS_SRC_MLPOINT_H_
