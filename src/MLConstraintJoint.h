#pragma once
#ifndef MUSCLEMASS_SRC_MLCONSTRAINTJOINT_H_
#define MUSCLEMASS_SRC_MLCONSTRAINTJOINT_H_

#include <vector>

#include "MLConstraint.h"

class MLBody;
class MLWorld;

class MLConstraintJoint : public MLConstraint
{
public:
	explicit MLConstraintJoint();
	virtual ~MLConstraintJoint();

	const MLBody * getBodyA() const { return m_bodyA; }
	const MLBody * getBodyB() const { return m_bodyB; }
	const Eigen::Matrix4d * getMatAtoJ0() const { return &m_aToJ0; }
	const Eigen::Matrix4d * getMatBtoJ0() const { return &m_bToJ0; }
	int getMatBtoJ(Eigen::Matrix4d *bToJ) const;

	void setPressure(double pressure) { m_pressure = pressure; }
	double getPressure() const { return m_pressure; }

	virtual MLError load(const nlohmann::json &elem, MLWorld *world, const std::string &folder);
	virtual void save(std::ofstream &ofs);

	virtual void createVirtualContacts();

	virtual void init();
	virtual void draw();
	virtual unsigned int getNumRowsBilateral() const;
	virtual unsigned int getNumColsDof() const;

	virtual int updateCache(const Eigen::VectorXd *q);

	virtual int evalBilateral(
		Eigen::VectorXd *g,
		Eigen::MatrixXd *G,
		Eigen::VectorXd *g_,
		Eigen::MatrixXd *G_,
		const Eigen::VectorXd *q);

	virtual int detectCollisions(
		std::vector<MLContact *> */*collisions*/,
		const Eigen::VectorXd */*q*/);

	virtual bool hasPressure() const { return true; }

protected:
	const MLBody *m_bodyA;
	const MLBody *m_bodyB;
	std::vector<unsigned int> m_rows;

	Eigen::Matrix4d m_aToJ0;
	Eigen::Matrix4d m_bToJ0;

	// Cached terms
	Eigen::Matrix4d m_aToW;
	Eigen::Matrix4d m_bToW;

	double m_pressure;

};



#endif // MUSCLEMASS_SRC_MLCONSTRAINTJOINT_H_
