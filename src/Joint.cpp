#include "Joint.h"
#include <iostream>

using namespace std;
using namespace Eigen;

Joint::Joint() {
	E_P_J.setZero();
	E_C_J.setZero();
	theta = 0.0;
};

Joint::~Joint() {

}

Eigen::Matrix4d Joint::getE_P_J() const {
	return this->E_P_J;
}


Eigen::Matrix4d Joint::getE_C_J() const {
	return this->E_C_J;
}

double Joint::getTheta() const {
	return theta;
}

void Joint::setE_C_J(Matrix4d _E_C_J) {
	this->E_C_J = _E_C_J;
}

void Joint::setE_P_J(Matrix4d _E_P_J) {
	this->E_P_J = _E_P_J;
}

void Joint::setTheta(double _theta) {
	this->theta = _theta;
}

void Joint::setE_C_J_0(Matrix4d _E_C_J_0) {
	this->E_C_J_0 = _E_C_J_0;
}

void Joint::setE_P_J_0(Eigen::Matrix4d _E_P_J_0) {
	this->E_P_J_0 = _E_P_J_0;
}

void Joint::setTheta_0(double _theta_0) {
	this->theta_0 = _theta_0;
}

void Joint::reset() {
	E_C_J = E_C_J_0;
	E_P_J = E_P_J_0;
	theta = theta_0;
}