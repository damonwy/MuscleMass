#include "Solver.h"
#include "Rigid.h"
#include "Joint.h"

#include <iostream>

using namespace std;
using namespace Eigen;

Solver::Solver(vector< shared_ptr<Rigid> > _boxes) {
	this->boxes = _boxes;

	A.resize(6 * (int)boxes.size(), 6 * (int)boxes.size());
	x.resize(6 * (int)boxes.size());
	b.resize(6 * (int)boxes.size());

	A.setZero();
	x.setZero();
	b.setZero();

}

void Solver::step(double h) {
	// Solve linear system
	for (int i = 0; i < (int)boxes.size(); i++) {
		auto box = boxes[i];

		A.block<6, 6>(6 * i, 6 * i) = box->getMassMatrix();
		b.segment<6>(6 * i) = box->getMassMatrix() * box->getTwist() + h * box->getForce();
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