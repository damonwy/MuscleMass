#include "Solver.h"
#include "Rigid.h"
#include "Joint.h"

#include <iostream>

using namespace std;
using namespace Eigen;

Solver::Solver(vector< shared_ptr<Rigid> > _boxes, bool _isReduced) {
	this->boxes = _boxes;
	this->isReduced = _isReduced;
	double n;

	if (isReduced) {
		double m = 6 * (int)boxes.size();
		n = 6 + 1 * ((int)boxes.size() - 1);
		M.resize(m, m);
		J.resize(m, n);
		f.resize(m);
		M.setZero();
		J.setZero();
		f.setZero();
	}
	else {
		n = 6 * (int)boxes.size() + 6 + 5 * ((int)boxes.size()-1);
	}

	A.resize(n, n);
	x.resize(n);
	b.resize(n);

	A.setZero();
	x.setZero();
	b.setZero();
}

void Solver::step(double h) {
	// Solve linear system
	if (isReduced) {
		int j = 0;
		for (int i = 0; i < (int)boxes.size(); i++) {
			auto box = boxes[i];

			M.block<6, 6>(6 * i, 6 * i) = box->getMassMatrix();
			f.segment<6>(6 * i) = box->getMassMatrix() * box->getTwist() + h * box->getForce();

			if (i == 0) {
				Matrix6d I;
				I.setIdentity();
				J.block<6, 6>(0, 0) = I;
			}
			else {
				auto joint = box->getJoint();
				Matrix6d Ad_C_W = Rigid::adjoint(box->getE().inverse());
				Matrix6d Ad_W_P = Rigid::adjoint(box->getParent()->getE());
				Matrix6d Ad_C_P = Ad_C_W * Ad_W_P;


				Matrix6d Ad_J_P = Rigid::adjoint(joint->getE_P_J().inverse());
				Matrix6d Ad_J_C = -Rigid::adjoint(joint->getE_C_J().inverse());
				int id_P = box->getParent()->getIndex();
				int id_C = box->getIndex();

			}

		}

		A = J.transpose() * M * J;
		b = J.transpose() * f;
		x = A.ldlt().solve(b);

		// Update Boxes
		VectorXd phi = J.transpose() * x;
		for (int i = 0; i < (int)boxes.size(); i++) {
			auto box = boxes[i];
			box->setTwist(phi.segment<6>(6 * i));
		}
	}
	else {
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
}

Solver::~Solver() {

}