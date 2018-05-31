#include <iostream>
#include <fstream>

#include <json.hpp>
#include "Scene.h"
#include "Particle.h"
#include "Shape.h"
#include "Program.h"
#include "Rigid.h"
#include "WrapSphere.h"
#include "WrapCylinder.h"
#include "WrapDoubleCylinder.h"
#include "Solver.h"
#include "Joint.h"
#include "MatlabDebug.h"
#include "Vector.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;


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
	//read a JSON file
	ifstream i(RESOURCE_DIR + "input.json");
	i >> js;
	i.close();

	// Units: meters, kilograms, seconds
	h = js["h"];
	grav << 0.0, -9.8, 0.0;
	Integrator time_integrator = RKF45;

	// Init boxes
	boxShape = make_shared<Shape>();
	boxShape->loadMesh(RESOURCE_DIR + "box2.obj");

	Matrix3d R;
	//R.setIdentity();
	R << 0, -1, 0,
		1, 0, 0,
		0, 0, 1;

	Vector3d p;
	p.setZero();

	Vector3d dimension = Vector3d(1.0, 4.0, 1.0);
	double mass = js["mass"];
	double scale = js["scale"];

	auto box0 = make_shared<Rigid>(boxShape, R, p, dimension, scale, mass, js["isReduced"]);
	box0->setIndex(0);
	boxes.push_back(box0);

	p += Vector3d(4.0, 0.0, 0.0);
	auto box1 = make_shared<Rigid>(boxShape, R, p, dimension, scale, mass, js["isReduced"]);
	box1->setIndex(1);
	box1->setParent(box0);
	//box0->addChild(box1);
	boxes.push_back(box1);

	p += Vector3d(4.0, 0.0, 0.0);
	auto box2 = make_shared<Rigid>(boxShape, R, p, dimension, scale, mass, js["isReduced"]);
	box2->setIndex(2);
	box2->setParent(box1);
	//box1->addChild(box2);
	boxes.push_back(box2);

	R.setIdentity();
	p.setZero();
	/*auto box3 = make_shared<Rigid>(boxShape, R, p, dimension, scale, mass, js["isReduced"]);
	box3->setIndex(3);
	boxes.push_back(box3);*/

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

	// Init ODE params
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

	// Init Particles
	sphereShape = make_shared<Shape>();
	sphereShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
	auto wc_p = make_shared<Particle>(sphereShape);
	wc_p->x0 = Vector3d(1.0, 2.0, 0.0);
	wc_p->r = 0.1;
	wc_p->update(box1->getE());
	box1->addPoint(wc_p);

	auto wc_s = make_shared<Particle>(sphereShape);
	wc_s->x0 = Vector3d(1.0, -2.0, 0.0);
	wc_s->r = 0.1;
	wc_s->update(box2->getE());
	box2->addPoint(wc_s);

	double cylinder_radius = js["cylinder_radius"];
	auto wc_o = make_shared<Particle>(sphereShape);
	wc_o->x0 << cylinder_radius + box2->getDimension()(0), -0.5 * box2->getDimension()(1) + 1.0, 0.0;
	wc_o->r = 0.1;
	wc_o->update(box2->getE());
	box2->addPoint(wc_o);

	auto wdc_u = make_shared<Particle>(sphereShape);
	wdc_u->x0 << cylinder_radius + box1->getDimension()(0), -0.5 * box1->getDimension()(1) + 1.0, 0.0;
	wdc_u->r = 0.1;
	wdc_u->update(box1->getE());
	box1->addPoint(wdc_u);

	auto wdc_v = make_shared<Particle>(sphereShape);
	wdc_v->x0 << cylinder_radius + box2->getDimension()(0), -0.5 * box2->getDimension()(1) + 1.0, 0.0;
	wdc_v->r = 0.1;
	wdc_v->update(box2->getE());
	box2->addPoint(wdc_v);

	auto wc_z = make_shared<Vector>();
	wc_z->dir0 << 0.0, 0.0, 1.0;
	wc_z->dir = wc_z->dir0;
	wc_z->setP(wc_o);
	wc_z->update(box2->getE());

	auto wdc_z_u = make_shared<Vector>();
	wdc_z_u->dir0 << 0.0, 0.0, -1.0;
	wdc_z_u->dir = wdc_z_u->dir0;
	wdc_z_u->setP(wdc_u);
	wdc_z_u->update(box1->getE());

	auto wdc_z_v = make_shared<Vector>();
	wdc_z_v->dir0 << 0.0, 0.0, 1.0;
	wdc_z_v->dir = wdc_z_v->dir0;
	wdc_z_v->setP(wdc_v);
	wdc_z_v->update(box2->getE());

	// Init WrapCylinder
	cylinderShape = make_shared<Shape>();
	cylinderShape->loadMesh(RESOURCE_DIR + "cylinder2.obj");
	
	Vector3d O = Vector3d(cylinder_radius + box2->getDimension()(0), -0.5 * box2->getDimension()(1) + 1.0, 0.0);
	
	R << 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, -1.0, 0.0;

	Matrix4d Ewrap;
	Ewrap.setIdentity();
	Ewrap.block<3, 3>(0, 0) = R;
	Ewrap.block<3, 1>(0, 3) = O;

	int num_points_on_arc = js["num_points_on_arc"];

	auto wrap_cylinder0 = make_shared<WrapCylinder>(cylinderShape, O, R, cylinder_radius, num_points_on_arc);
	wrap_cylinder0->setE(box2->getE() * Ewrap);
	wrap_cylinder0->setP(wc_p);
	wrap_cylinder0->setS(wc_s);
	wrap_cylinder0->setO(wc_o);
	wrap_cylinder0->setZ(wc_z);

	box2->addCylinder(wrap_cylinder0);

	// Init WrapDoubleCylinder
	auto wrap_doublecylinder = make_shared<WrapDoubleCylinder>(cylinderShape, cylinder_radius, cylinder_radius, num_points_on_arc);
	wrap_doublecylinder->setP(wc_p);
	wrap_doublecylinder->setS(wc_s);
	wrap_doublecylinder->setU(wdc_u);
	wrap_doublecylinder->setV(wdc_v);
	wrap_doublecylinder->setZ_U(wdc_z_u);
	wrap_doublecylinder->setZ_V(wdc_z_v);
	wrap_doublecylinder->setE_U(box1->getE() * Ewrap);
	wrap_doublecylinder->setE_P_U(Ewrap);
	wrap_doublecylinder->setE_V(box2->getE() * Ewrap);
	wrap_doublecylinder->setE_P_V(Ewrap);
	wrap_doublecylinder->setParent_U(box1);
	wrap_doublecylinder->setParent_V(box2);

	box2->addDoubleCylinder(wrap_doublecylinder); // Only add double cylinder to the latest updated rigid body, so that all the positions are updated
	box2->setCylinderStatus(js["isCylinder"]);
	box2->setDoubleCylinderStatus(js["isDoubleCylinder"]);

	auto ws_p = make_shared<Particle>(sphereShape);
	ws_p->x0 = Vector3d(1.0, 2.0, 0.0);
	ws_p->r = 0.1;
	ws_p->update(box1->getE());
	box1->addPoint(ws_p);

	auto ws_s = make_shared<Particle>(sphereShape);
	ws_s->x0 = Vector3d(1.0, -2.0, 0.0);
	ws_s->r = 0.1;
	ws_s->update(box2->getE());
	box2->addPoint(ws_s);

	auto ws_o = make_shared<Particle>(sphereShape);
	ws_o->x0 = Vector3d(cylinder_radius + box2->getDimension()(0), -0.5 * box2->getDimension()(1) + 1.0, 0.0);
	ws_o->r = 0.1;
	ws_o->update(box2->getE());
	box2->addPoint(ws_o);

	// Init WrapSphere
	auto wrap_sphere = make_shared<WrapSphere>(cylinder_radius, num_points_on_arc);

	wrap_sphere->setP(ws_p);
	wrap_sphere->setS(ws_s);
	wrap_sphere->setO(ws_o);
	
	wrap_sphere->setParent(box2);
	box2->addSphere(wrap_sphere);

	// Init solver	
	solver = make_shared<Solver>(boxes, js["isReduced"], time_integrator);
}

void Scene::init()
{
	boxShape->init();
	cylinderShape->init();
	sphereShape->init();
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
		}
	}

	t += h;
	step_i += 1;
	computeEnergy();

	int plot_steps = js["plot_steps"];
	saveData(plot_steps);
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

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> prog2, shared_ptr<MatrixStack> P) const
{
	for (int i = 0; i < (int)boxes.size(); ++i) {
		boxes[i]->draw(MV, prog, prog2, P);
	}
}
