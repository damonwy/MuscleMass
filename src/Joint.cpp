#include "Joint.h"
#include <iostream>

using namespace std;
using namespace Eigen;

Joint::Joint():
dtheta(0.0), theta(0.0), theta_0(0.0), min_theta(0.0), max_theta(0.0){
	E_P_J.setIdentity();
	E_C_J.setIdentity();
	E_C_J_0.setIdentity();
	E_P_J_0.setIdentity();
};

Joint::Joint(Matrix4d _E_P_J, Matrix4d _E_C_J, double _dtheta, double _min_theta, double _max_theta) :
	dtheta(_dtheta), min_theta(_min_theta), max_theta(_max_theta),
	E_C_J(_E_C_J), E_C_J_0(_E_C_J), E_P_J(_E_P_J), E_P_J_0(_E_P_J){
	this->theta_0 = 0.0;
	this->theta = 0.0;

}


void Joint::reset() {
	E_C_J = E_C_J_0;
	E_P_J = E_P_J_0;
	theta = theta_0;
}

Joint::~Joint() {

}