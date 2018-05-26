#include <iostream>

#include "Scene.h"
#include "Particle.h"
#include "Shape.h"
#include "Program.h"
#include "Rigid.h"
#include "Solver.h"
#include "Joint.h"
#include "MatlabDebug.h"

using namespace std;
using namespace Eigen;

Scene::Scene() :
	t(0.0),
	h(1e-2),
	step_i(0),
	grav(0.0, 0.0, 0.0),
	t_start(0.0),
	t_stop(100.0),
	n_step(10),
	K(0.0),
	V(0.0)
{
	theta_list.resize(n_step);
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{
	// Units: meters, kilograms, seconds
	h = 1e-2;
	grav << 0.0, -9.8, 0.0;
	bool isReduced = true;
	double scale = 1.0;
	double mass = 8.0;
	Integrator time_integrator = SYMPLECTIC;

	// Init boxes
boxShape = make_shared<Shape>();
boxShape->loadMesh(RESOURCE_DIR + "box2.obj");

Matrix3d R;
//R.setIdentity();
R << 0, -1, 0,
1, 0, 0,
0, 0, 1;
Vector3d p = Vector3d(0.0, 0.0, 0.0);
Vector3d dimension = Vector3d(1.0, 4.0, 1.0);

auto box0 = make_shared<Rigid>(boxShape, R, p, dimension, scale, mass, isReduced);
box0->setIndex(0);
boxes.push_back(box0);

p += Vector3d(4.0, 0.0, 0.0);
auto box1 = make_shared<Rigid>(boxShape, R, p, dimension, scale, mass, isReduced);
box1->setIndex(1);
box1->setParent(box0);
//box0->addChild(box1);
boxes.push_back(box1);

p += Vector3d(4.0, 0.0, 0.0);
auto box2 = make_shared<Rigid>(boxShape, R, p, dimension, scale, mass, isReduced);
box2->setIndex(2);
box2->setParent(box1);
//box1->addChild(box2);
boxes.push_back(box2);

// Init joints
auto joint1 = box1->getJoint();
Matrix4d E_P_J;
E_P_J.setIdentity();
E_P_J(1, 3) = -2.0;
Matrix4d E_C_J = box1->getE().inverse() * box0->getE() * E_P_J;

joint1->setE_P_J_0(E_P_J);
joint1->setE_C_J_0(E_C_J);
joint1->setTheta_0(0.0);
joint1->reset();

auto joint2 = box2->getJoint();
E_C_J = box2->getE().inverse() * box1->getE() * E_P_J;

joint2->setE_P_J_0(E_P_J);
joint2->setE_C_J_0(E_C_J);
joint2->setTheta_0(0.0);
joint2->reset();

y.push_back(0.0);
y.push_back(0.0);
y.push_back(0.0);
y.push_back(0.0);
y.resize(2 * (boxes.size() - 1));

yp.push_back(0.0);
yp.push_back(0.0);
yp.push_back(0.0);
yp.push_back(0.0);
yp.resize(2 * (boxes.size() - 1));

// Init solver	
solver = make_shared<Solver>(boxes, isReduced, time_integrator);
}

void Scene::init()
{
	boxShape->init();
}

void Scene::tare()
{
	for (int i = 0; i < (int)boxes.size(); ++i) {
		boxes[i]->tare();
	}
}

void Scene::reset()
{
	t = 0.0;
	for (int i = 0; i < (int)boxes.size(); ++i) {
		boxes[i]->reset();
	}
}

void Scene::step()
{
	if (solver->time_integrator == RKF45) {

		vector<double> tmp = solver->solve(&y[0], &yp[0], (int)2 * (boxes.size() - 1), t, t + h);

		y.clear();

		for (int i = 0; i < (int)2 * (boxes.size() - 1); i++) {
			y.push_back(tmp[i]);
		}
		VectorXd thetadot;
		thetadot.resize(boxes.size() - 1);

		for (int i = 0; i < boxes.size() - 1; i++){
			thetadot(i) = y[boxes.size() - 1 + i];
		}

		for (int i = 0; i < (int)boxes.size(); i++) {
			auto box = boxes[i];
			if (i != 0) {
				// Update joint angles
				box->setJointAngle(y[i-1], true);	
			}
		}

		MatrixXd J = solver->getJ_twist_thetadot();
		MatrixXd JJ = J.block(0, solver->n - solver->num_joints, solver->m, solver->num_joints);

		VectorXd phi = JJ * thetadot;
		for (int i = 0; i < (int)boxes.size(); i++) {
			auto box = boxes[i];
			if (i != 0) {
				box->setTwist(phi.segment<6>(6 * i));
			}
		}
		// solve_once

		//if (step_i == 0) {
		//	double y[2] = { 0.0, 0.0 };
		//	double yp[2];
		//	theta_list = solver->solve_once(y, yp, 2, t_start, t_stop, n_step);
		//	for (int i = 0; i < (int)boxes.size(); i++) {
		//		auto box = boxes[i];
		//		if (i != 0) {
		//			// Update joint angles
		//			box->setJointAngle(theta_list(step_i), true);
		//		}
		//	}
		//}
		//else {
		//	for (int i = 0; i < (int)boxes.size(); i++) {
		//		auto box = boxes[i];
		//		if (i != 0) {
		//		// Update joint angles
		//		box->setJointAngle(theta_list(step_i), true);	
		//		}
		//	}
		//}

	}
	else if (solver->time_integrator == SYMPLECTIC) {
		solver->step(h);
		for (int i = 0; i < (int)boxes.size(); ++i) {
			boxes[i]->step(h);
			cout << boxes[i]->getTwist() << endl;
		}
	}

	t += h;
	step_i += 1;
	computeEnergy();
	saveData(10000);
}


void Scene::saveData(int num_steps) {
	// Save data and plot in MATLAB
	Kvec.push_back(K);
	Vvec.push_back(V);
	Tvec.push_back(step_i);

	if (step_i == num_steps) {
		cout << "finished" << endl;
		VectorXd Kv;
		VectorXd Vv;
		VectorXd Tv;
		Kv.resize(Kvec.size());
		Vv.resize(Vvec.size());
		Tv.resize(Tvec.size());

		for (int i = 0; i < Kvec.size(); i++) {
			Kv(i) = Kvec[i];
			Vv(i) = Vvec[i];
			Tv(i) = Tvec[i];
		}

		vec_to_file(Kv, "K");
		vec_to_file(Vv, "V");
		vec_to_file(Tv, "T");
	}
}

void Scene::computeEnergy() {
	K = 0.0;
	V = 0.0;

	for (int i = 0; i < (int)boxes.size(); ++i) {
		double vi = boxes[i]->m * grav.transpose() * boxes[i]->getP();
		V += vi;
		double ki = 0.5 * boxes[i]->getTwist().transpose() * boxes[i]->getMassMatrix() * boxes[i]->getTwist();
		K += ki;
	}
}


void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	for (int i = 0; i < (int)boxes.size(); ++i) {
		boxes[i]->draw(MV, prog);
	}
}
