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
	return E_P_J;
}


Eigen::Matrix4d Joint::getE_C_J() const {
	return E_C_J;
}

double Joint::getTheta() const {
	return theta;
}

void Joint::setE_C_J(Matrix4d _E_C_J) {
	E_C_J = _E_C_J;
}

void Joint::setE_P_J(Matrix4d _E_P_J) {
	E_P_J = _E_P_J;
}

void Joint::setTheta(double _theta) {
	theta = _theta;
}