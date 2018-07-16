#include "MLConstraintJoint.h"



#include "MLBody.h"
#include "MLWorld.h"
#include "MLContact.h"

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

