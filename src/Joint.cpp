#include "Joint.h"
#include "Rigid.h"
#include <iostream>

using namespace std;
using namespace Eigen;

Joint::Joint():
dtheta(0.0), theta(0.0), theta_0(0.0), thetadot(0.0), thetaddot(0.0), min_theta(0.0), max_theta(0.0){
	E_P_J.setIdentity();
	E_C_J.setIdentity();
	E_C_J_0.setIdentity();
	E_P_J_0.setIdentity();
};

Joint::Joint(Matrix4d _E_P_J, Matrix4d _E_C_J, double _min_theta, double _max_theta) :
	dtheta(0.0), thetadot(0.0), thetaddot(0.0), min_theta(_min_theta), max_theta(_max_theta),
	E_C_J(_E_C_J), E_C_J_0(_E_C_J), E_P_J(_E_P_J), E_P_J_0(_E_P_J){
	this->theta_0 = 0.0;
	this->theta = 0.0;
}

VectorXd Joint::getThetaVector(vector<shared_ptr<Joint>> joints) {
	VectorXd thetalist;
	thetalist.resize((int)joints.size());
	thetalist.setZero();
	for (int i = 0; i < (int)joints.size(); ++i) {
		thetalist(i) = joints[i]->getTheta();
	}
	return thetalist;
}

VectorXd Joint::getThetadotVector(vector<shared_ptr<Joint>> joints) {
	VectorXd thetadotlist;
	thetadotlist.resize((int)joints.size());
	thetadotlist.setZero();
	for (int i = 0; i < (int)joints.size(); ++i) {
		thetadotlist(i) = joints[i]->getThetadot();
	}
	return thetadotlist;
}

void Joint::setThetadotVector(vector <shared_ptr<Joint>> joints, VectorXd thetadotlist) {	
	for (int i = 0; i < (int)joints.size(); ++i) {
		joints[i]->setThetadot(thetadotlist(i));
	}
}

void Joint::reset() {
	E_C_J = E_C_J_0;
	E_P_J = E_P_J_0;
	//theta = theta_0;
	//thetadot = 0.0;
	//thetaddot = 0.0;
}

Joint::~Joint() {

}