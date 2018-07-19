#pragma once
#ifndef MUSCLEMASS_SRC_MLBODY_H_
#define MUSCLEMASS_SRC_MLBODY_H_

#include <vector>

#include <Eigen/Dense>

#include "MLCommon.h"
#include "MLObject.h"


class MLRigid;
class MLWrench;
class MatrixStack;
class Program;

class MLBody : public MLObject {
public:
	explicit MLBody();
	virtual ~MLBody();

	virtual MLError load(const nlohmann::json &elem, MLWorld *world, const std::string &folder);
	virtual void save(std::ofstream &ofs);

	const Eigen::Matrix4d * getTransfromCached() const { return &m_mat; }
	const Matrix6d * getDiagMassMatrix() const { return m_massMatrix; }
	void setExternalForce(MLWrench *fext) { m_fext = fext; }
	MLWrench * getExternalForce() { return m_fext; }

	virtual void init();
	virtual void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog);
	void update();

	std::string getInfoStr() const;

	void setTransfromInit(Eigen::Matrix4d *E);

private:
	Eigen::Matrix4d m_mat;
	MLWrench *m_fext;
	Matrix6d *m_massMatrix;

};

#endif // MUSCLEMASS_SRC_MLBODY_H_
