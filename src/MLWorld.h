#pragma once
#ifndef MUSCLEMASS_SRC_MLWORLD_H_
#define MUSCLEMASS_SRC_MLWORLD_H_


#include <map>
#include <vector>

#include "MLCommon.h"

class MLComp;
class MLBody;
class MLConstraint;
class MLError;


class MLWorld {
public:
	explicit MLWorld();
	virtual ~MLWorld();

	virtual MLError load(const nlohmann::json &elem, const std::string &folder);
	virtual void save(std::ofstream &ofs);

	void addBody(std::shared_ptr<MLBody> body);
	void addComp(std::shared_ptr<MLComp> comp);
	void addConstraint(std::shared_ptr<MLConstraint> constraint);
	

	void init();
	void step();
	void draw();


protected:
	// Actual objects created
	std::vector<std::shared_ptr<MLBody>> m_bodies;
	std::vector<std::shared_ptr<MLComp>> m_comps;
	std::vector<std::shared_ptr<MLConstraint>> m_constraints;
	
	typedef std::map<int, std::shared_ptr<MLComp>> MapCompUID;
	typedef std::map<std::string, std::shared_ptr<MLBody>> MapBodyName;
	typedef std::map<int, std::shared_ptr<MLBody>> MapBodyUID;


	MapCompUID m_compUID;
	MapBodyName m_bodyName;
	MapBodyUID m_bodyUID;

};

#endif // MUSCLEMASS_SRC_MLWORLD_H_
