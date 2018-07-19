#pragma once
#ifndef MUSCLEMASS_SRC_MLASSEMBLER_H_
#define MUSCLEMASS_SRC_MLASSEMBLER_H_

#include <Eigen/Dense>
#include <vector>
#include <map>

#include "MLCommon.h"

class MLBody;
class MLConstraintJoint;
class Scene;
class MLWorld;

class MLAssembler {
public:
	explicit MLAssembler();
	virtual ~MLAssembler();

	void clear();

	void setNumRows(const std::string &rowType, unsigned int numRows);
	unsigned int getNumRows(const std::string &rowType) const;

	Eigen::VectorXd * getNewVector(const std::string &rowType);
	Eigen::MatrixXd * getNewMatrix(const std::string &rowType, const std::string &colType, const std::string &matType);

protected:
	typedef std::map<std::string, unsigned int> MapRow;
	typedef std::map<std::string, std::vector<Eigen::VectorXd *>> MapVec;
	typedef std::map<std::pair<std::string, std::string>, std::vector<Eigen::MatrixXd *>> MapMat;

	MapRow m_numRows;
	MapVec m_vecs;
	MapMat m_mats;



};

#endif // MUSCLEMASS_SRC_MLASSEMBLER_H_
