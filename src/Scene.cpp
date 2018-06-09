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

#include "Joint.h"
#include "MatlabDebug.h"
#include "Vector.h"
#include "JsonEigen.h"
#include "Spring.h"
#include "SymplecticIntegrator.h"
#include "RKF45Integrator.h"
#include "Solver.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

#include <unsupported/Eigen/MatrixFunctions> // TODO: avoid using this later, write a func instead

Scene::Scene() :
	t(0.0),
	h(1e-2),
	step_i(0),
	grav(0.0, 0.0, 0.0),
	t_start(0.0),
	t_stop(100.0),
	n_step(10000),
	K(0.0),
	V(0.0)
{
	theta_list.resize(n_step * 2);
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{	
	// Init shapes
	boxShape = make_shared<Shape>();
	boxShape->loadMesh(RESOURCE_DIR + "box2.obj");
	sphereShape = make_shared<Shape>();
	sphereShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
	cylinderShape = make_shared<Shape>();
	cylinderShape->loadMesh(RESOURCE_DIR + "cylinder2.obj");

	//read a JSON file
	ifstream i(RESOURCE_DIR + "input.json");
	i >> js;
	i.close();

	// Units: meters, kilograms, seconds
	h = js["h"];
	Eigen::from_json(js["grav"], grav);
	time_integrator = SYMPLECTIC;

	// Init boxes	
	auto box0 = addBox(js["Rz"], js["p0"], js["dimension"], js["scale"], js["mass"], boxShape, js["isReduced"], 0);	// No parent
	auto box1 = addBox(js["Rz"], js["p1"], js["dimension"], js["scale"], js["mass"], boxShape, js["isReduced"], 1, box0);
	auto box2 = addBox(js["Rz"], js["p2"], js["dimension"], js["scale"], js["mass"], boxShape, js["isReduced"], 2, box1);

	// Init joints
	addJoint(box0, box1, js["E_P_J_1"], js["min_theta_1"], js["max_theta_1"]);
	addJoint(box1, box2, js["E_P_J_2"], js["min_theta_2"], js["max_theta_2"]);

	// Init Particles
	auto s_p = make_shared<Particle>(sphereShape);
	from_json(js["s_p_x0"], s_p->x0);
	s_p->r = js["s_p_r"];
	s_p->update(box0->getE());
	box0->addPoint(s_p);

	// Attached to springs
	auto s_s = make_shared<Particle>(sphereShape);
	from_json(js["s_s_x0"], s_s->x0);
	s_s->r = js["s_s_r"];
	s_s->update(box1->getE());
	box1->addPoint(s_s);

	// Attached to wrapCylinder
	auto wc_p = make_shared<Particle>(sphereShape);
	from_json(js["wc_p_x0"], wc_p->x0);
	wc_p->r = js["wc_p_r"];
	wc_p->update(box1->getE());
	box1->addPoint(wc_p);

	auto wc_s = make_shared<Particle>(sphereShape);
	from_json(js["wc_s_x0"], wc_s->x0);
	wc_s->r = js["wc_s_r"];
	wc_s->update(box2->getE());
	box2->addPoint(wc_s);

	double cylinder_radius = js["cylinder_radius"];
	auto wc_o = make_shared<Particle>(sphereShape);
	wc_o->x0 << cylinder_radius + box2->getDimension()(0), -0.5 * box2->getDimension()(1) + 1.0, 0.0;
	wc_o->r = js["wc_o_r"];
	wc_o->update(box2->getE());
	box2->addPoint(wc_o);

	// Attached to wrapDoubleCylinder
	auto wdc_u = make_shared<Particle>(sphereShape);
	wdc_u->x0 << cylinder_radius + box1->getDimension()(0), -0.5 * box1->getDimension()(1) + 1.0, 0.0;
	wdc_u->r = js["wdc_u_r"];
	wdc_u->update(box1->getE());
	box1->addPoint(wdc_u);

	auto wdc_v = make_shared<Particle>(sphereShape);
	wdc_v->x0 << cylinder_radius + box2->getDimension()(0), -0.5 * box2->getDimension()(1) + 1.0, 0.0;
	wdc_v->r = js["wdc_v_r"];
	wdc_v->update(box2->getE());
	box2->addPoint(wdc_v);

	auto wc_z = make_shared<Vector>();
	from_json(js["wc_z_dir0"], wc_z->dir0);
	wc_z->dir = wc_z->dir0;
	wc_z->setP(wc_o);
	wc_z->update(box2->getE());

	auto wdc_z_u = make_shared<Vector>();
	from_json(js["wdc_z_u_dir0"], wdc_z_u->dir0);
	wdc_z_u->dir = wdc_z_u->dir0;
	wdc_z_u->setP(wdc_u);
	wdc_z_u->update(box1->getE());

	auto wdc_z_v = make_shared<Vector>();
	from_json(js["wdc_z_v_dir0"], wdc_z_v->dir0);
	wdc_z_v->dir = wdc_z_v->dir0;
	wdc_z_v->setP(wdc_v);
	wdc_z_v->update(box2->getE());

	// Init WrapCylinder
	Vector3d O = Vector3d(cylinder_radius + box2->getDimension()(0), -0.5 * box2->getDimension()(1) + 1.0, 0.0);

	Matrix3d R;
	from_json(js["Rx"], R);
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
	box2->setSphereStatus(js["isSphere"]);

	auto ws_p = make_shared<Particle>(sphereShape);
	from_json(js["ws_p_x0"], ws_p->x0);
	ws_p->r = js["ws_p_r"];
	ws_p->update(box1->getE());
	box1->addPoint(ws_p);

	auto ws_s = make_shared<Particle>(sphereShape);
	from_json(js["ws_s_x0"], ws_s->x0);
	ws_s->r = js["ws_s_r"];
	ws_s->update(box2->getE());
	box2->addPoint(ws_s);

	auto ws_o = make_shared<Particle>(sphereShape);
	ws_o->x0 = Vector3d(cylinder_radius + box2->getDimension()(0), -0.5 * box2->getDimension()(1) + 1.0, 0.0);
	ws_o->r = js["ws_o_r"];
	ws_o->update(box2->getE());
	box2->addPoint(ws_o);

	// Init WrapSphere
	auto wrap_sphere = make_shared<WrapSphere>(cylinder_radius, num_points_on_arc);

	wrap_sphere->setP(ws_p);
	wrap_sphere->setS(ws_s);
	wrap_sphere->setO(ws_o);

	wrap_sphere->setParent(box2);
	box2->addSphere(wrap_sphere);

	// Init Spring
	double spring_mass = js["spring_mass"];
	if (js["isSpring"]) {
		auto spring = make_shared<Spring>(s_p, s_s, spring_mass);
		springs.push_back(spring);
	}

	if (time_integrator == SYMPLECTIC) {
		symplectic_solver = make_shared<SymplecticIntegrator>(boxes, springs, js["isReduced"]);
	}
	else if (time_integrator == RKF45) {
		rkf45_solver = make_shared<RKF45Integrator>(boxes, springs, js["isReduced"]);
	}
}

shared_ptr<Joint> Scene::addJoint(shared_ptr<Rigid> parent, shared_ptr<Rigid> child, json jE_P_J, json jmin_theta, json jmax_theta) {
	Matrix4d E_P_J;
	Eigen::from_json(jE_P_J, E_P_J);
	Matrix4d E_C_J = child->getE().inverse() * parent->getE() * E_P_J;
	auto joint = make_shared<Joint>(E_P_J, E_C_J, jmin_theta.get<double>() / 180.0 * PI, jmax_theta.get<double>() / 180.0 * PI);
	child->setJoint(joint);
	return joint;
}

shared_ptr<Rigid> Scene::addBox(json _R, json _p, json _dimension, json _scale, json _mass, const shared_ptr<Shape> shape, json _isReduced, int id, shared_ptr<Rigid> parent) {
	Matrix3d R;	
	Vector3d p, dimension;
	Eigen::from_json(_R, R);
	Eigen::from_json(_p, p);
	Eigen::from_json(_dimension, dimension);
	double scale = _scale;
	double mass = _mass;
	auto box = make_shared<Rigid>(shape, R, p, dimension, scale, mass, _isReduced);
	box->setIndex(id);
	box->setParent(parent);
	boxes.push_back(box);
	return box;
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
	if (time_integrator == SYMPLECTIC) {
		symplectic_solver->step(h);
		for (int i = 0; i < (int)boxes.size(); ++i) {
			boxes[i]->step(h);

			if (i != 0) {
				this->phi.segment<6>(6 * (i - 1)) = boxes[i]->getTwist();
			}
		}
		MatrixXd jj = symplectic_solver->getJ_twist_thetadot();
	}
	else if (time_integrator == RKF45) {
		// nothing

	}else{
		// nothing
	}

	t += h;
	step_i += 1;
	computeEnergy();

	if (js["isPlotEnergy"]) {
		int plot_steps = js["plot_steps"];
		saveData(plot_steps);
	}
}

void Scene::saveData(int num_steps) {
	// Save data and plot in MATLAB
	Kvec.push_back(K);
	Vvec.push_back(V);
	Tvec.push_back(step_i);
	Twist_vec.push_back(phi);

	if (step_i == num_steps) {
		cout << "finished" << endl;
		VectorXd Kv;
		VectorXd Vv;
		VectorXd Tv;
		VectorXd Twistv;

		Kv.resize(Kvec.size());
		Vv.resize(Vvec.size());
		Tv.resize(Tvec.size());
		Twistv.resize(12 * Twist_vec.size());

		for (int i = 0; i < Kvec.size(); i++) {
			Kv(i) = Kvec[i];
			Vv(i) = Vvec[i];
			Tv(i) = Tvec[i];
			Twistv.segment<12>(12 * i) = Twist_vec[i];
		}

		vec_to_file(Kv, "K");
		vec_to_file(Vv, "V");
		vec_to_file(Tv, "T");
		vec_to_file(Twistv, "Twistv");
	}
}

void Scene::computeEnergy() {
	K = 0.0;
	V = 0.0;

	for (int i = 0; i < (int)boxes.size(); ++i) {
		double vi = boxes[i]->m * grav.transpose() * boxes[i]->getP();

		//cout << "pos: " << boxes[i]->getP()(1) << endl;
		V += vi;
		double ki = 0.5 * boxes[i]->getTwist().transpose() * boxes[i]->getMassMatrix() * boxes[i]->getTwist();
		//cout << "twist: " << boxes[i]->getTwist() << endl;
		K += ki;
	}
	if (step_i == 1) {
		V0 = V;
		K0 = K;
	}
	//cout << V0 - V << endl;
	//cout << K - K0 << endl;
	//cout << V0 - V + K - K0 << endl;

}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> prog2, shared_ptr<MatrixStack> P) const
{
	for (int i = 0; i < (int)boxes.size(); ++i) {
		boxes[i]->draw(MV, prog, prog2, P);
	}

	for (int i = 0; i < (int)springs.size(); ++i) {
		springs[i]->draw(MV, prog, prog2, P);
	}
}
