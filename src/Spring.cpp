#include "Spring.h"

#include <iostream>
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Particle.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"

using namespace std;
using namespace Eigen;

Spring::Spring(shared_ptr<Particle> p0, shared_ptr<Particle> p1, double _mass) :
	E(1.0), mass(_mass)
{
	assert(p0);
	assert(p1);
	assert(p0 != p1);
	this->p0 = p0;
	this->p1 = p1;
	Vector3d x0 = p0->x;
	Vector3d x1 = p1->x;
	Vector3d dx = x1 - x0;
	L = dx.norm();
	this->p0_b = x0;
	this->p1_b = x1;
	assert(L > 0.0);
}

double Spring::computeLength() {
	return (p0->x - p1->x).norm();
}

void Spring::setPosBeforePert() {
	this->p0_b = p0->x_temp;
	this->p1_b = p1->x_temp;
}

void Spring::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const {

	prog->bind();
	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 0.0, 0.0).data());
	glUniform3fv(prog->getUniform("kdBack"), 1, Vector3f(1.0, 1.0, 0.0).data());

	// Draw P, Spoints
	this->p0->draw(MV, prog);
	this->p1->draw(MV, prog);

	prog->unbind();

	// Draw wrapping
	prog2->bind();
	glUniformMatrix4fv(prog2->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog2->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	MV->pushMatrix();
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	glColor3f(0.0, 0.0, 0.0); // black
	glLineWidth(3);
	Vector3d p = p0->x;
	Vector3d s = p1->x;

	glBegin(GL_LINE_STRIP);
	glVertex3f(p(0), p(1), p(2));
	glVertex3f(s(0), s(1), s(2));
	glEnd();
	MV->popMatrix();
	prog2->unbind();
}

Spring::~Spring()
{
	
}
