#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Rigid.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"

using namespace std;

Rigid::Rigid() :
	r(1.0),
	m(1.0),
	i(-1),
	x(0.0, 0.0, 0.0),
	v(0.0, 0.0, 0.0),
	fixed(true)
{

}

Rigid::Rigid(const shared_ptr<Shape> s) :
	r(1.0),
	m(1.0),
	i(-1),
	x(0.0, 0.0, 0.0),
	v(0.0, 0.0, 0.0),
	fixed(true),
	normal(0.0, 0.0, 0.0),
	box(s)
{

}

Rigid::~Rigid()
{
}

void Rigid::tare()
{
	x0 = x;
	v0 = v;
}

void Rigid::reset()
{
	x = x0;
	v = v0;
}

void Rigid::step(double h) {

}

void Rigid::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	if (box) {
		MV->pushMatrix();
		MV->translate(x(0), x(1), x(2));
		MV->scale(r);
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		box->draw(prog);
		MV->popMatrix();
	}
}
