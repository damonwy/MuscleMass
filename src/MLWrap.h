#pragma once
#ifndef MUSCLEMASS_SRC_MLWRAP_H_
#define MUSCLEMASS_SRC_MLWRAP_H_

#include "MLCommon.h"
#include "MLObject.h"

class MLBody;
class MLError;

class MLWrap :public MLObject {
public:
	explicit MLWrap();
	virtual ~MLWrap();

	virtual MLError load(const const nlohmann::json &elem, MLWorld *world, const std::string &folder);
	virtual void save(std::ofstream &ofs);

	virtual void init();
	virtual void draw();

protected:
	
};

#endif // MUSCLEMASS_SRC_MLWRAP_H_