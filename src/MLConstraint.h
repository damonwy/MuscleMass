#pragma once
#ifndef MUSCLEMASS_SRC_MLCONSTRAINT_H_
#define MUSCLEMASS_SRC_MLCONSTRAINT_H_

#include "MLCommon.h"
#include "MLObject.h"

class MLContact;

class MLConstraint : public MLObject {
public:
	explicit MLConstraint();
	virtual ~MLConstraint();

	virtual void createVirtualContacts() = 0;

	virtual void init();
	virtual void draw();

	virtual unsigned int getNumColsDof() const = 0;
	virtual unsigned int getNumRowsBilateral() const = 0;
	unsigned int getNumRowsContact() const;
	unsigned int getNumRowsFriction() const;

	void setIndexBilateral(unsigned int *indexBilateral);
	void setIndexContact(unsigned int *indexContact);
	void setIndexFriction(unsigned int *indexFriction);

	virtual int updateCache(const Eigen::VectorXd *q);

	virtual int evalBilateral(
		Eigen::VectorXd *g,
		Eigen::MatrixXd *G,
		Eigen::VectorXd *g_,
		Eigen::MatrixXd *G_,
		const Eigen::VectorXd *q) = 0;

	/*virtual int evalContact(
		Eigen::VectorXd *n,
		Eigen::MatrixXd *N,
		Eigen::VectorXd *n_,
		Eigen::MatrixXd *N_,
		Eigen::VectorXd *P_,
		const Eigen::VectorXd *q);*/

	/*virtual int evalFriction(
		Eigen::MatrixXd *D,
		const Eigen::VectorXd *q,
		const Eigen::VectorXd *v);

	int findActiveContacts(
		std::vector<unsigned int> *ib,
		std::vector<unsigned int> *bActive,
		std::vector<unsigned int> *aActive,
		std::vector<double> *muac,
		const Eigen::VectorXd *a);*/

	virtual int detectCollisions(
		std::vector<MLContact *> */*collisions*/,
		const Eigen::VectorXd */*q*/) { return 1; }

	virtual bool hasPressure() const { return false; }

	unsigned int getIndexBilateral() const { return m_indexBilateral; }
	unsigned int getIndexContact() const { return m_indexContact; }
	unsigned int getIndexFriction() const { return m_indexFriction; }

	std::string getColLabel() const;
	std::string getRowLabel(char prefix) const;

	void setMu(double mu) { m_mu = mu; }
	double getMu() const { return m_mu; }

protected:
	unsigned int m_indexBilateral; // G
	unsigned int m_indexContact;   // N
	unsigned int m_indexFriction;  // D (both active and inactive)
	unsigned int m_indexFriction_; // D_ (only active)
	double m_mu;
	std::vector<MLContact *> m_contacts;

};

#endif // MUSCLEMASS_SRC_MLCONSTRAINT_H_