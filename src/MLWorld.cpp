#include "MLWorld.h"

#include <fstream>
#include <iomanip>

#include "MLError.h"
#include "MLBody.h"

#include <iostream>
#include <fstream>

#include "MLComp.h"
#include "MLConstraint.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

MLWorld::MLWorld() {

}

MLWorld::~MLWorld() {

}

MLError MLWorld::load(const json &elem, const std::string &folder) {
	return MLError();
}

void MLWorld::save(std::ofstream &ofs) {

}

void MLWorld::addBody(std::shared_ptr<MLBody> body) {

}

void MLWorld::addComp(std::shared_ptr<MLComp> comp) {
	m_comps.push_back(comp);
	m_compUID[comp->getUID()] = comp;
}

void MLWorld::addConstraint(std::shared_ptr<MLConstraint> constraint) {

}

void MLWorld::init() {
	for (int i = 0; i < (int)m_bodies.size(); ++i) {
		m_bodies[i]->init();
	}
	for (int i = 0; i < (int)m_constraints.size(); ++i) {
		m_constraints[i]->init();
	}
	
}

void MLWorld::step() {

}

void MLWorld::draw() {

}
