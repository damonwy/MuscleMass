#pragma once
#ifndef MUSCLEMASS_SRC_JOINT_H_
#define MUSCLEMASS_SRC_JOINT_H_

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

struct Joint {
	
	Eigen::Matrix4d E_P_J;	// Where the joint is wrt parent
	Eigen::Matrix4d E_C_J;	// Where the joint is wrt children
	double theta;			// Current joint angle

	//int type;		// the type of joint
	//int hinge_type;
	Joint();
};

#endif // MUSCLEMASS_SRC_JOINT_H_