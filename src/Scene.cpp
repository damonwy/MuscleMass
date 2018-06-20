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
	auto joint1 = addJoint(box0, box1, js["E_P_J_1"], js["min_theta_1"], js["max_theta_1"]);
	auto joint2 = addJoint(box1, box2, js["E_P_J_2"], js["min_theta_2"], js["max_theta_2"]);

	// Init Spring
	auto spring0 = addSpring(js["s_p_x0"], box1, js["s_s_x0"], box2, js["spring_mass"]);

	// Init WrapCylinder
	auto wc0 = addWrapCylinder(js["wc_p_x0"], box1, js["wc_s_x0"], box2, js["wc_o_x0"], box2, js["cylinder_radius"], js["wc_z_dir0"]);
	
	// Init WrapDoubleCylinder
	auto wdc0 = addWrapDoubleCylinder(js["wdc_p_x0"], box1, js["wdc_s_x0"], box2, js["wdc_u_x0"], box1, js["wdc_v_x0"], box2, js["cylinder_radius"], js["cylinder_radius"], js["wdc_z_u_dir0"], js["wdc_z_v_dir0"]);

	// Init WrapSphere
	auto ws0 = addSphere(js["ws_p_x0"], box1, js["ws_s_x0"], box2, js["ws_o_x0"], box2, js["ws_o_r"]);

	box2->setCylinderStatus(js["isCylinder"]);
	box2->setDoubleCylinderStatus(js["isDoubleCylinder"]);
	box2->setSphereStatus(js["isSphere"]);

	if (time_integrator == SYMPLECTIC) {
		symplectic_solver = make_shared<SymplecticIntegrator>(boxes, joints, springs, js["isReduced"], js["num_samples_on_muscle"], js["grav"], js["epsilon"]);
	}
	else if (time_integrator == RKF45) {
		rkf45_solver = make_shared<RKF45Integrator>(boxes, springs, js["isReduced"]);
	}

	for (int i = 0; i < (int)springs.size(); ++i) {
		springs[i]->step(joints);
	}
}

shared_ptr<Joint> Scene::addJoint(shared_ptr<Rigid> parent, shared_ptr<Rigid> child, json jE_P_J, json jmin_theta, json jmax_theta) {
	Matrix4d E_P_J;
	Eigen::from_json(jE_P_J, E_P_J);
	Matrix4d E_C_J = child->getE().inverse() * parent->getE() * E_P_J;
	auto joint = make_shared<Joint>(E_P_J, E_C_J, jmin_theta.get<double>() / 180.0 * PI, jmax_theta.get<double>() / 180.0 * PI);
	child->setJoint(joint);
	joint->setChild(child);
	joint->setParent(parent);
	joints.push_back(joint);
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
	auto box = make_shared<Rigid>(shape, R, p, dimension, scale, mass, _isReduced, js["grav"]);
	box->setIndex(id);
	box->setParent(parent);
	boxes.push_back(box);
	return box;
}

shared_ptr<Spring> Scene::addSpring(json jp_x, shared_ptr<Rigid> p_parent, json js_x, shared_ptr<Rigid> s_parent, json jmass) {
	auto p = make_shared<Particle>(sphereShape);
	from_json(jp_x, p->x0);
	p->r = js["particle_r"];
	p->setParent(p_parent);
	p->update(p_parent->getE());
	p_parent->addPoint(p);

	auto s = make_shared<Particle>(sphereShape);
	from_json(js_x, s->x0);
	s->r = js["particle_r"];
	s->setParent(s_parent);
	s->update(s_parent->getE());
	s_parent->addPoint(s);

	// Init Spring
	double spring_mass = jmass;
	auto spring = make_shared<Spring>(p, s, spring_mass, js["num_samples_on_muscle"], grav, js["epsilon"], js["isReduced"], js["stiffness"]);
	if (js["isSpring"]) {	
		springs.push_back(spring);		
	}
	return spring;
}

shared_ptr<WrapSphere> Scene::addSphere(json jp_x, shared_ptr<Rigid> p_parent, json js_x, shared_ptr<Rigid> s_parent, json jo_x, shared_ptr<Rigid> o_parent, json jradius) {
	auto ws_p = make_shared<Particle>(sphereShape);
	from_json(jp_x, ws_p->x0);
	ws_p->r = js["particle_r"];
	ws_p->update(p_parent->getE());
	p_parent->addPoint(ws_p);

	auto ws_s = make_shared<Particle>(sphereShape);
	from_json(js_x, ws_s->x0);
	ws_s->r = js["particle_r"];
	ws_s->update(s_parent->getE());
	s_parent->addPoint(ws_s);

	double s_radius = jradius;
	Vector3d o_x;
	from_json(jo_x, o_x);
	auto ws_o = make_shared<Particle>(sphereShape);
	ws_o->x0 = Vector3d(s_radius + o_parent->getDimension()(0), -0.5 * o_parent->getDimension()(1) + 1.0, 0.0);
	ws_o->x0 += o_x;
	ws_o->r = s_radius;
	ws_o->update(o_parent->getE());
	o_parent->addPoint(ws_o);

	// Init WrapSphere
	auto wrap_sphere = make_shared<WrapSphere>(s_radius, js["num_samples_on_muscle"]);

	wrap_sphere->setP(ws_p);
	wrap_sphere->setS(ws_s);
	wrap_sphere->setO(ws_o);

	wrap_sphere->setParent(s_parent);
	s_parent->addSphere(wrap_sphere);
	wrap_spheres.push_back(wrap_sphere);
	return wrap_sphere;
}

shared_ptr<WrapCylinder> Scene::addWrapCylinder(json jp_x, shared_ptr<Rigid> p_parent, json js_x, shared_ptr<Rigid> s_parent, json jo_x, shared_ptr<Rigid> o_parent, json jradius, json jzdir) {
	auto wc_p = make_shared<Particle>(sphereShape);
	from_json(jp_x, wc_p->x0);
	wc_p->r = js["particle_r"];
	wc_p->setParent(p_parent);
	wc_p->update(p_parent->getE());
	p_parent->addPoint(wc_p);

	auto wc_s = make_shared<Particle>(sphereShape);
	from_json(js_x, wc_s->x0);
	wc_s->r = js["particle_r"];
	wc_s->setParent(s_parent);
	wc_s->update(s_parent->getE());
	s_parent->addPoint(wc_s);

	double cylinder_radius = jradius;
	auto wc_o = make_shared<Particle>(sphereShape);
	Vector3d o_x;
	from_json(jo_x, o_x);
	wc_o->x0 << cylinder_radius + o_parent->getDimension()(0), -0.5 * o_parent->getDimension()(1) + 1.0, 0.0;
	wc_o->x0 += o_x;
	wc_o->r = js["particle_r"];
	wc_o->setParent(o_parent);
	wc_o->update(o_parent->getE());
	o_parent->addPoint(wc_o);

	int num_points_on_arc = js["num_points_on_arc"];
	Matrix3d R;
	from_json(js["Rx"], R);
	Matrix4d Ewrap;
	Ewrap.setIdentity();
	Ewrap.block<3, 3>(0, 0) = R;
	Ewrap.block<3, 1>(0, 3) = wc_o->x0;

	auto wc_z = make_shared<Vector>();
	from_json(jzdir, wc_z->dir0);
	wc_z->dir = wc_z->dir0;
	wc_z->setP(wc_o);
	wc_z->update(o_parent->getE());

	auto wrap_cylinder = make_shared<WrapCylinder>(cylinderShape, wc_o->x0, R, cylinder_radius, num_points_on_arc);
	wrap_cylinder->setE(o_parent->getE() * Ewrap);
	wrap_cylinder->setP(wc_p);
	wrap_cylinder->setS(wc_s);
	wrap_cylinder->setO(wc_o);
	wrap_cylinder->setZ(wc_z);

	o_parent->addCylinder(wrap_cylinder);

	wrap_cylinders.push_back(wrap_cylinder);
	return wrap_cylinder;
}

shared_ptr<WrapDoubleCylinder> Scene::addWrapDoubleCylinder(json jp_x, shared_ptr<Rigid> p_parent, json js_x, shared_ptr<Rigid> s_parent, json ju_x, shared_ptr<Rigid> u_parent, json jv_x, shared_ptr<Rigid> v_parent, json juradius, json jvradius, json jzudir, json jzvdir) {
	auto wdc_p = make_shared<Particle>(sphereShape);
	from_json(jp_x, wdc_p->x0);
	wdc_p->r = js["particle_r"];
	wdc_p->setParent(p_parent);
	wdc_p->update(p_parent->getE());
	p_parent->addPoint(wdc_p);

	auto wdc_s = make_shared<Particle>(sphereShape);
	from_json(js_x, wdc_s->x0);
	wdc_s->r = js["particle_r"];
	wdc_s->setParent(s_parent);
	wdc_s->update(s_parent->getE());
	s_parent->addPoint(wdc_s);

	double u_radius = juradius;
	double v_radius = jvradius;

	auto wdc_u = make_shared<Particle>(sphereShape);
	Vector3d u_x;
	from_json(ju_x, u_x);
	wdc_u->x0 << u_radius + u_parent->getDimension()(0), -0.5 * u_parent->getDimension()(1) + 1.0, 0.0;
	wdc_u->x0 += u_x;
	wdc_u->r = js["particle_r"];
	wdc_u->update(u_parent->getE());
	u_parent->addPoint(wdc_u);

	auto wdc_v = make_shared<Particle>(sphereShape);
	Vector3d v_x;
	from_json(jv_x, v_x);
	wdc_v->x0 << v_radius + v_parent->getDimension()(0), -0.5 * v_parent->getDimension()(1) + 1.0, 0.0;
	wdc_v->x0 += v_x;
	wdc_v->r = js["particle_r"];
	wdc_v->update(v_parent->getE());
	v_parent->addPoint(wdc_v);

	auto wdc_z_u = make_shared<Vector>();
	from_json(jzudir, wdc_z_u->dir0);
	wdc_z_u->dir = wdc_z_u->dir0;
	wdc_z_u->setP(wdc_u);
	wdc_z_u->update(u_parent->getE());

	auto wdc_z_v = make_shared<Vector>();
	from_json(jzvdir, wdc_z_v->dir0);
	wdc_z_v->dir = wdc_z_v->dir0;
	wdc_z_v->setP(wdc_v);
	wdc_z_v->update(v_parent->getE());

	Matrix3d R;
	from_json(js["Rx"], R);
	Matrix4d Ewrapu, Ewrapv;
	Ewrapu.setIdentity();
	Ewrapu.block<3, 3>(0, 0) = R;
	Ewrapu.block<3, 1>(0, 3) = wdc_u->x0;

	Ewrapv.setIdentity();
	Ewrapv.block<3, 3>(0, 0) = R;
	Ewrapv.block<3, 1>(0, 3) = wdc_v->x0;

	//// Init WrapDoubleCylinder
	auto wrap_doublecylinder = make_shared<WrapDoubleCylinder>(cylinderShape, juradius, jvradius, js["num_points_on_arc"]);
	wrap_doublecylinder->setP(wdc_p);
	wrap_doublecylinder->setS(wdc_s);
	wrap_doublecylinder->setU(wdc_u);
	wrap_doublecylinder->setV(wdc_v);
	wrap_doublecylinder->setZ_U(wdc_z_u);
	wrap_doublecylinder->setZ_V(wdc_z_v);
	wrap_doublecylinder->setE_U(u_parent->getE() * Ewrapu);
	wrap_doublecylinder->setE_P_U(Ewrapu);
	wrap_doublecylinder->setE_V(v_parent->getE() * Ewrapv);
	wrap_doublecylinder->setE_P_V(Ewrapv);
	wrap_doublecylinder->setParent_U(u_parent);
	wrap_doublecylinder->setParent_V(v_parent);
	wrap_doublecylinders.push_back(wrap_doublecylinder);

	s_parent->addDoubleCylinder(wrap_doublecylinder);// Only add double cylinder to the latest updated rigid body, so that all the positions are updated
	return wrap_doublecylinder;
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
				//this->phi.segment<6>(6 * (i - 1)) = boxes[i]->getTwist();
			}
		}
		//MatrixXd jj = symplectic_solver->getJ_twist_thetadot();
		for (int i = 0; i < (int)springs.size(); ++i) {
			springs[i]->step(joints);
		}
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

	// Rigid Body:
	for (int i = 0; i < (int)boxes.size(); ++i) {
		double vi = boxes[i]->getPotentialEnergy();
		double ki = boxes[i]->getKineticEnergy();
		V += vi;
		K += ki;
	}	

	// Spring:
	for (int i = 0; i < (int)springs.size(); ++i) {
		double vi = springs[i]->getPotentialEnergy();
		double ki = springs[i]->getKineticEnergy();
		V += vi;
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

	symplectic_solver->draw(MV, prog2, P);
}
