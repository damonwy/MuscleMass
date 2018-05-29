#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Vector.h"
#include "Program.h"
#include "MatrixStack.h"
#include "Particle.h"

using namespace std;
using namespace Eigen;

Vector::Vector() :
	dir(0.0, 0.0, 0.0),
	dir0(0.0, 0.0, 0.0),
	fixed(false)
{

}

Vector::~Vector()
{
}

void Vector::reset()
{
	dir = dir0;
}

void Vector::update(Matrix4d E) {
	Vector4d pos;
	pos.segment<3>(0) = this->dir0;
	pos(3) = 0.0;
	pos = E * pos;
	this->dir = pos.segment<3>(0);
}

void Vector::setP(shared_ptr<Particle> _p) {
	this->p = _p;
}

void Vector::draw(shared_ptr<MatrixStack> MV, shared_ptr<MatrixStack> P, const shared_ptr<Program> prog) const
{
	if (p) {
		prog->bind();
		glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		MV->pushMatrix();
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		glColor3f(0.0, 0.0, 0.0); // black
		glLineWidth(3);
		glBegin(GL_LINES);
		Vector3f p0 = p->x.cast<float>();
		glVertex3f(p0(0), p0(1), p0(2));
		Vector3f p1 = (p0 + this->dir.cast<float>());
		glVertex3f(p1(0), p1(1), p1(2));
		glEnd();
		MV->popMatrix();
		prog->unbind();
	}
}
