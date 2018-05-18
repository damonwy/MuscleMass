#pragma once
#ifndef MUSCLEMASS_SRC_JOINT_H_
#define MUSCLEMASS_SRC_JOINT_H_

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Joint {
private:
	Eigen::Matrix4d E_P_J;	// Where the joint is wrt parent
	Eigen::Matrix4d E_C_J;	// Where the joint is wrt children
	double theta;			// Current joint angle

	//int type;		// the type of joint
	//int hinge_type;
	
public:
	Eigen::Matrix4d getE_P_J() const;
	Eigen::Matrix4d getE_C_J() const;
	double getTheta() const;

	void setE_C_J(Eigen::Matrix4d _E_C_J);
	void setE_P_J(Eigen::Matrix4d _E_P_J);
	void setTheta(double _theta);

	Joint();
	~Joint();
};

#endif // MUSCLEMASS_SRC_JOINT_H_