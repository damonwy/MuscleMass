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
	//read a JSON file
	ifstream i(RESOURCE_DIR + "input.json");
	i >> js;
	i.close();

	// Units: meters, kilograms, seconds
	h = js["h"];
	Eigen::from_json(js["grav"], grav);
	time_integrator = SYMPLECTIC;

	// Init boxes
	boxShape = make_shared<Shape>();
	boxShape->loadMesh(RESOURCE_DIR + "box2.obj");

	Matrix3d R;
	Eigen::from_json(js["Rz"], R);
	//R.setIdentity();

	Vector3d p;
	p.setZero();

	Vector3d dimension = js["dimension"];
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
	//p -= Vector3d(2.0, 2.0, 0.0);
	//from_json(js["Rx"], R);
	//R.setIdentity();
	auto box2 = make_shared<Rigid>(boxShape, R, p, dimension, scale, mass, js["isReduced"]);
	box2->setIndex(2);
	box2->setParent(box1);
	//box1->addChild(box2);

	if (js["isCylinder"] || js["isDoubleCylinder"] || js["isSpring"] || js["isSphere"]) {
		// If I want to show wrapping, I need box2
		boxes.push_back(box2);
	}

	R.setIdentity();
	p.setZero();
	/*auto box3 = make_shared<Rigid>(boxShape, R, p, dimension, scale, mass, js["isReduced"]);
	box3->setIndex(3);
	boxes.push_back(box3);*/

	// Init joints
	
	Matrix4d E_P_J;
	E_P_J.setIdentity();
	E_P_J(1, 3) = -2.0;
	Matrix4d E_C_J = box1->getE().inverse() * box0->getE() * E_P_J;

	auto joint1 = make_shared<Joint>(E_P_J, E_C_J, double(js["theta_1"]) / 180.0 * PI, double(js["min_theta_1"]) / 180.0 * PI, double(js["max_theta_1"]) / 180.0 * PI);
	box1->setJoint(joint1);

	E_C_J = box2->getE().inverse() * box1->getE() * E_P_J;
	
	auto joint2 = make_shared<Joint>(E_P_J, E_C_J, double(js["theta_2"]) / 180.0 * PI, double(js["min_theta_2"]) / 180.0 * PI, double(js["max_theta_2"]) / 180.0 * PI);
	box2->setJoint(joint2);

	// Init ODE params
	if (js["isReduced"]) {
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
	}
	else {
		Matrix4d E1 = box1->getE();
		Matrix4d E2 = box2->getE();
		VectorXd e1 = Rigid::unbracket6(E1.log());
		VectorXd e2 = Rigid::unbracket6(E2.log());
		for (int i = 0; i < 6; i++) {
			y.push_back(e1(i));
			yp.push_back(0.0);
		}
		for (int i = 0; i < 6; i++) {
			y.push_back(e2(i));
			yp.push_back(0.0);
		}
		for (int i = 0; i < 12; i++) {
			y.push_back(0.0);
			yp.push_back(0.0);
		}

		y.resize(12 * (boxes.size() - 1));
		yp.resize(12 * (boxes.size() - 1));
	}
	

	// Init Particles
	sphereShape = make_shared<Shape>();
	sphereShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
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
	cylinderShape = make_shared<Shape>();
	cylinderShape->loadMesh(RESOURCE_DIR + "cylinder2.obj");
	
	Vector3d O = Vector3d(cylinder_radius + box2->getDimension()(0), -0.5 * box2->getDimension()(1) + 1.0, 0.0);
	
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

	// Init solver	
	//solver = make_shared<Solver>(boxes, springs, js["isReduced"], time_integrator);

	if (time_integrator == SYMPLECTIC) {
		symplectic_solver = make_shared<SymplecticIntegrator>(boxes, springs, js["isReduced"]);
	}
	else if (time_integrator == RKF45) {
		rkf45_solver = make_shared<RKF45Integrator>(boxes, springs, js["isReduced"]);
	}

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
//	if (solver->time_integrator == RKF45) {
//		if (js["isReduced"]) {
//			MatrixXd J = solver->getJ_twist_thetadot();
//			//cout << J << endl;
//			MatrixXd JJ = J.block(0, solver->n - solver->num_joints, solver->m, solver->num_joints);
//
//			for (int i = 0; i < (int)2 * (boxes.size() - 1); i++) {
//				
//				cout << "yinput" << y[i] << endl;
//			}
//
//			vector<double> tmp = solver->solve(&y[0], &yp[0], (int)2 * (boxes.size() - 1), t, t + h);
//			///vector<double> tmp = solver->solve(&y[0], &yp[0], (int)2 * (boxes.size() - 1), 0, t);
//
//			y.clear();
//
//			for (int i = 0; i < (int)2 * (boxes.size() - 1); i++) {
//				y.push_back(tmp[i]);
//				cout <<"youtput"<< tmp[i] << endl;
//			}
//
//			VectorXd thetadot;
//			thetadot.resize(boxes.size() - 1);
//
//			for (int i = 0; i < boxes.size() - 1; i++) {
//				thetadot(i) = y[boxes.size() - 1 + i];
//				//cout << thetadot(i) << endl;
//			}
//
//			for (int i = 0; i < (int)boxes.size(); i++) {
//				auto box = boxes[i];
//				if (i != 0) {
//					// Update joint angles
//					box->setJointAngle(y[i - 1], true);
//					//cout << box->getP()(1) << endl;
//				}
//			}
//
//			VectorXd phi = JJ * thetadot;
//			for (int i = 0; i < (int)boxes.size(); i++) {
//				auto box = boxes[i];
//				if (i != 0) {
//					box->setTwist(phi.segment<6>(6 * i));
//				}
//			}
//
//			//this->phi = phi.segment<12>(6);
//
//			// solve_once
//
//			//if (step_i == 0) {
//			//	double y[4] = { 0.0, 0.0, 0.0, 0.0 };
//			//	double yp[4];
//			//	theta_list = solver->solve_once(y, yp, 4, t_start, t_stop, n_step);
//			//	for (int i = 0; i < (int)boxes.size(); i++) {
//			//		auto box = boxes[i];
//			//		if (i != 0) {
//			//			// Update joint angles
//			//			box->setJointAngle(theta_list(2 * step_i + i - 1), true);
//			//			cout << theta_list(2 * step_i + i - 1) << endl;
//			//		}
//			//	}
//			//}
//			//else {
//			//	for (int i = 0; i < (int)boxes.size(); i++) {
//			//		auto box = boxes[i];
//			//		if (i != 0) {
//			//		// Update joint angles
//			//		box->setJointAngle(theta_list(2 * step_i + i - 1), true);	
//			//		cout << theta_list(step_i) << endl;
//			//		}
//			//	}
//			//}
//		}
//		else {
//			// Maximal Coordinate
//			vector<double> tmp = solver->solve(&y[0], &yp[0], (int)12 * (boxes.size() - 1), t, t + h);
//			y.clear();
//
//			for (int i = 0; i < (int)12 * (boxes.size() - 1); i++) {
//				y.push_back(tmp[i]);
//				//cout << tmp[i] << endl;
//			}
//
//			Vector6d e1, e2, phi1, phi2;
//			// Only need twist to update boxes
//			
//			for (int i = 12; i < 18; i++) {
//				phi1(i - 12) = y[i];
//			}
//			for (int i = 18; i < 24; i++) {
//				phi2(i - 18) = y[i];
//			}
//
//			// Update twists
//			for (int i = 0; i < (int)boxes.size(); i++) {
//				// Update E
//				boxes[i]->step(h);
//			}
//			boxes[1]->setTwist(phi1);
//			boxes[2]->setTwist(phi2);
//			
//			Matrix4d E1 = boxes[1]->getE();
//			Matrix4d E2 = boxes[2]->getE();
//			e1 = Rigid::unbracket6(E1.log());
//			e2 = Rigid::unbracket6(E2.log());
//
//			for (int i = 0; i < 6; i++) {
//				y[i] = e1(i);
//				y[6 + i] = e2(i);
//			}
//			//cout << "phi2:" << endl << phi2 << endl;
//		}		
//	}
//	else if (solver->time_integrator == SYMPLECTIC) {
//		solver->step(h);
//		for (int i = 0; i < (int)boxes.size(); ++i) {
//			boxes[i]->step(h);
//			
//			if (i != 0) {
//				this->phi.segment<6>(6 * (i - 1)) = boxes[i]->getTwist();
//				//cout << boxes[i]->getP()(1) << endl;
//
//			}
//		}
//		MatrixXd jj = solver->getJ_twist_thetadot();
//		///cout << jj << endl;
//	}
//
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
	//cout << "V" << V << endl;
	//cout << "K" << K << endl;
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
