#include "MLComp.h"

#include <iostream>
#include <fstream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "MLError.h"
#include "MLBody.h"
#include "MLWorld.h"
#include "Shape.h"
#include "MatrixStack.h"
#include "Program.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

MLComp::MLComp(std::shared_ptr<MLBody> parent) {
	m_parent = parent;
}

MLComp::~MLComp() {

}

MLError MLComp::load(const nlohmann::json &elem, std::shared_ptr<MLWorld> world, const std::string &folder) {
	auto error = MLObject::load(elem, world, folder);
	try {
		// Add pointer to this component to world
		world->addComp(shared_from_this());
		// Read in cToW

		// get *mToW


	}
	catch (json::exception& e) {
		// output exception information
		std::cout << "message: " << e.what() << '\n' << "exception id: " << e.id << std::endl;
		return MLError("error parsing object");
	}

	// read shape file in folder
	std::string objFile = folder + "/" + m_name + ".obj";
	m_shape->loadMesh(objFile);
	return MLError();

}

void MLComp::save(std::ofstream &ofs) {
}

void MLComp::init() {

}

void MLComp::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P) {

	prog->bind();
	if (m_shape) {
		Eigen::Matrix4d cToW;
		const Eigen::Matrix4d *pToW = m_parent->getTransfromCached();
		cToW = *pToW * m_cToP;
		
		glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
		MV->pushMatrix();
		glUniform3f(prog->getUniform("lightPos1"), 1.0, 1.0, 1.0);
		glUniform1f(prog->getUniform("intensity_1"), 0.8);
		glUniform3f(prog->getUniform("lightPos2"), -1.0, 1.0, 1.0);
		glUniform1f(prog->getUniform("intensity_2"), 0.2);
		glUniform1f(prog->getUniform("s"), 200);
		glUniform3f(prog->getUniform("ka"), 0.2, 0.2, 0.2);
		glUniform3f(prog->getUniform("kd"), 0.8, 0.7, 0.7);
		glUniform3f(prog->getUniform("ks"), 1.0, 0.9, 0.8);
		MV->multMatrix(eigen_to_glm(cToW));
	
		//MV->scale(r);
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_shape->draw(prog);
		MV->popMatrix();
	}
	prog->unbind();
	// draw box
}

MLError MLComp::getTransformWorld(Eigen::Matrix4d &sToWord, const Eigen::Matrix4d *pToWorld) const {
	sToWord = *pToWorld * m_cToP;
	return MLError();
}