#include "Solver.h"
#include "Rigid.h"
#include "Joint.h"

#include <iostream>

using namespace std;
using namespace Eigen;

Solver::Solver(vector< shared_ptr<Rigid> > _boxes) {
	this->boxes = _boxes;
	double n = 6 * (int)boxes.size() + 6 + 5 * (int)boxes.size();
	A.resize(n, n);
	x.resize(n);
	b.resize(n);

	A.setZero();
	x.setZero();
	b.setZero();
}

void Solver::step(double h) {
	// Solve linear system
	int j = 0;
	for (int i = 0; i < (int)boxes.size(); i++) {
		auto box = boxes[i];

		A.block<6, 6>(6 * i, 6 * i) = box->getMassMatrix();
		b.segment<6>(6 * i) = box->getMassMatrix() * box->getTwist() + h * box->getForce();

		//
		if (i == 0) {
			Matrix6d I;
			I.setIdentity();
			A.block<6, 6>(6 * (int)boxes.size(), 0) = I;
			A.block<6, 6>(0, 6 * (int)boxes.size()) = I;
		}
		else {
			auto joint = box->getJoint();
			Matrix6d Ad_J_P = Rigid::adjoint(joint->getE_P_J().inverse());
			Matrix6d Ad_J_C = -Rigid::adjoint(joint->getE_C_J().inverse());
			int id_P = box->getParent()->getIndex();
			int id_C = box->getIndex();

			Matrix5x6d G_J_P;
			G_J_P.block<2, 6>(0, 0) = Ad_J_P.block<2, 6>(0, 0);
			G_J_P.block<3, 6>(2, 0) = Ad_J_P.block<3, 6>(3, 0);

			Matrix5x6d G_J_C;
			G_J_C.block<2, 6>(0, 0) = Ad_J_C.block<2, 6>(0, 0);
			G_J_C.block<3, 6>(2, 0) = Ad_J_C.block<3, 6>(3, 0);

			A.block<5, 6>(6 * (int)boxes.size() + 6 + 5 * j, 6 * id_P) = G_J_P;
			A.block<6, 5>(6 * id_P, 6 * (int)boxes.size() + 6 + 5 * j) = G_J_P.transpose();

			A.block<5, 6>(6 * (int)boxes.size() + 6 + 5 * j, 6 * id_C) = G_J_C;
			A.block<6, 5>(6 * id_C, 6 * (int)boxes.size() + 6 + 5 * j) = G_J_C.transpose();

			j++;
		}
	}

	x = A.ldlt().solve(b);

	// Update boxes
	for (int i = 0; i < (int)boxes.size(); i++) {
		auto box = boxes[i];
		box->setTwist(x.segment<6>(6 * i));
	}
}

Solver::~Solver() {

}