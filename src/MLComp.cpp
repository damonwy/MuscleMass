#include "MLComp.h"

#include <iostream>
#include <fstream>

#include "MLError.h"
#include "MLBody.h"
#include "MLWorld.h"

using json = nlohmann::json;

MLComp::MLComp(std::shared_ptr<MLBody> parent) {
	m_parent = parent;
}

MLComp::~MLComp() {

}

MLError MLComp::load(const nlohmann::json &elem, std::shared_ptr<MLWorld> world, const std::string &folder) {


}

void MLComp::save(std::ofstream &ofs) {
}

void MLComp::init() {

}

void MLComp::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog) {

	// draw box
}

MLError MLComp::getTransformWorld(Eigen::Matrix4d &sToWord, const Eigen::Matrix4d *pToWorld) const {
	sToWord = *pToWorld * m_cToP;
	return MLError();
}