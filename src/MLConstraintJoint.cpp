#include "MLConstraintJoint.h"

#include "MLBody.h"
#include "MLWorld.h"
#include "MLContact.h"
#include "MLError.h"

MLConstraintJoint::MLConstraintJoint() :
	m_bodyA(NULL),
	m_bodyB(NULL),
	m_pressure(0.0) 
{

}

MLConstraintJoint::~MLConstraintJoint() {

}

void MLConstraintJoint::init() {
	MLConstraint::init();
}

void MLConstraintJoint::draw() {
	MLConstraint::draw();
}

MLError MLConstraintJoint::load(const nlohmann::json &elem, MLWorld *world, const std::string &folder) {
	return MLError();
}

void MLConstraintJoint::save(std::ofstream &ofs) {

}

void MLConstraintJoint::createVirtualContacts() {

}

int MLConstraintJoint::updateCache(const Eigen::VectorXd *q) {
	return 0;
}


int MLConstraintJoint::evalBilateral(
	Eigen::VectorXd *g,
	Eigen::MatrixXd *G,
	Eigen::VectorXd *g_,
	Eigen::MatrixXd *G_,
	const Eigen::VectorXd *q) {
	return 0;

}

int MLConstraintJoint::getMatBtoJ(Eigen::Matrix4d *bToJ) const {
	return 0;
}


unsigned int MLConstraintJoint::getNumRowsBilateral() const {
	return 0;
}

unsigned int MLConstraintJoint::getNumColsDof() const {
	return 0;
}

int MLConstraintJoint::detectCollisions(
	std::vector<MLContact *> */*collisions*/,
	const Eigen::VectorXd */*q*/) {
	return 0;
}