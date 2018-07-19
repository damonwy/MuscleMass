#include "MLCompCylinder.h"

#include <iostream>
#include <fstream>

#include "MatrixStack.h"
#include "Program.h"

#include "MLError.h"
#include "MLBody.h"
#include "MLWorld.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

MLCompCylinder::MLCompCylinder(std::shared_ptr<MLBody> parent) : MLComp(parent) {

}

MLCompCylinder::~MLCompCylinder() {

}

MLError MLCompCylinder::load(const json &elem, std::shared_ptr<MLWorld> world, const std::string &folder) {
	auto error = MLComp::load(elem, world, folder);
	try {
		m_radius = elem.at("radius").get<double>();
		m_height = elem.at("height").get<double>();
	}
	catch (json::exception& e) {
		// output exception information
		std::cout << "message: " << e.what() << '\n' << "exception id: " << e.id << std::endl;
		return MLError("error parsing object");
	}

	return MLError();
}

void MLCompCylinder::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P) const {

}

void MLCompCylinder::save(std::ofstream &ofs) {

}