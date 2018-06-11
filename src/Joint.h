#pragma once
#ifndef MUSCLEMASS_SRC_JOINT_H_
#define MUSCLEMASS_SRC_JOINT_H_

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <iostream>
class Joint {
private:
	Eigen::Matrix4d E_P_J;	// Where the joint is wrt parent
	Eigen::Matrix4d E_C_J;	// Where the joint is wrt children
	double dtheta;			// Current rotation angle
	double theta;			// Current joint angle
	double max_theta;		
	double min_theta;

	Eigen::Matrix4d E_P_J_0;
	Eigen::Matrix4d E_C_J_0;
	double theta_0;
	
public:
	Eigen::Matrix4d getE_P_J() const { return this->E_P_J; }
	Eigen::Matrix4d getE_C_J() const { return this->E_C_J; }
	double getTheta() const { return this->theta; }
	double getDTheta() const { return this->dtheta; }
	double getMaxTheta() const { return this->max_theta; }
	double getMinTheta() const { return this->min_theta; }

	void setE_C_J(Eigen::Matrix4d _E_C_J) { this->E_C_J = _E_C_J; }
	void setE_P_J(Eigen::Matrix4d _E_P_J) { this->E_P_J = _E_P_J; }
	void setDTheta(double _dtheta) { this->dtheta = _dtheta; this->theta += _dtheta;}
	void setTheta(double _theta) { this->theta = _theta; }
	void setE_C_J_0(Eigen::Matrix4d _E_C_J_0) { this->E_C_J_0 = _E_C_J_0; }
	void setE_P_J_0(Eigen::Matrix4d _E_P_J_0) { this->E_P_J_0 = _E_P_J_0; }
	void setTheta_0(double _theta_0) { this->theta_0 = _theta_0; }
	void reset();
	Joint();
	Joint(Eigen::Matrix4d _E_P_J, Eigen::Matrix4d _E_C_J, double _min_theta, double _max_theta);

	~Joint();
};

#endif // MUSCLEMASS_SRC_JOINT_H_