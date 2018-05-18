#include <iostream>

#include "Scene.h"
#include "Particle.h"
#include "Shape.h"
#include "Program.h"
#include "Rigid.h"
#include "Solver.h"
#include "Joint.h"

using namespace std;
using namespace Eigen;

Scene::Scene() :
	t(0.0),
	h(1e-2),
	grav(0.0, 0.0, 0.0)
{
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{
	// Units: meters, kilograms, seconds
	h = 1e-3;	
	grav << 0.0, -9.8, 0.0;
	bool isReduced = true;
	double scale = 1.0;
	double mass = 8.0; 

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

	joint1->setE_P_J(E_P_J);
	joint1->setE_C_J(E_C_J);

	auto joint2 = box2->getJoint();
	E_C_J = box2->getE().inverse() * box1->getE() * E_P_J;

	joint2->setE_P_J(E_P_J);
	joint2->setE_C_J(E_C_J);


	// Init solver
	solver = make_shared<Solver>(boxes, isReduced);

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
	t += h;
	solver->step(h);

	for (int i = 0; i < (int)boxes.size(); ++i) {
		boxes[i]->step(h);
	}
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	for (int i = 0; i < (int)boxes.size(); ++i) {
		boxes[i]->draw(MV, prog);
	}
}
