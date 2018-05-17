#include <iostream>

#include "Scene.h"
#include "Particle.h"
#include "Shape.h"
#include "Program.h"
#include "Rigid.h"
#include "Solver.h"

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
	double scale = 1.0;
	double mass = 8.0; 

	auto box = make_shared<Rigid>(boxShape, R, p, dimension, scale, mass);
	boxes.push_back(box);

	// Init solver
	solver = make_shared<Solver>(boxes);

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
