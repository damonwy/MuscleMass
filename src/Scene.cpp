#include <iostream>

#include "Scene.h"
#include "Particle.h"
#include "Shape.h"
#include "Program.h"
#include "Rigid.h"

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
	h = 1e-2;	
	grav << 0.0, -9.8, 0.0;

	boxShape = make_shared<Shape>();
	boxShape->loadMesh(RESOURCE_DIR + "box2.obj");

	auto box = make_shared<Rigid>(boxShape);
	boxes.push_back(box);
	// init box params
	box->r = 1.0;
	box->x << 0.0, 0.0, 0.0;

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
