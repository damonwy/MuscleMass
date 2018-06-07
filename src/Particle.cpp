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

Particle::Particle() :
	r(1.0),
	m(1.0),
	i(-1),
	x(0.0, 0.0, 0.0),
	v(0.0, 0.0, 0.0),
	fixed(false)
{
	
}

Particle::Particle(const shared_ptr<Shape> s) :
	r(1.0),
	m(1.0),
	i(-1),
	x(0.0, 0.0, 0.0),
	v(0.0, 0.0, 0.0),
	fixed(false),
	normal(0.0,0.0,0.0),
	sphere(s)
{
	
}

Particle::~Particle()
{
}

void Particle::tare()
{
	x0 = x;
	v0 = v;
}

void Particle::reset()
{
	x = x0;
	v = v0;
}

void Particle::update(Matrix4d E) {
	Vector4d pos;
	pos.segment<3>(0) = this->x0;
	pos(3) = 1.0;

	pos = E * pos;
	this->x = pos.segment<3>(0);
}

void Particle::updateTemp(Matrix4d E) {
	Vector4d pos;
	pos.segment<3>(0) = this->x0;
	pos(3) = 1.0;

	pos = E * pos;
	this->x_temp = pos.segment<3>(0);
}

void Particle::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	if(sphere) {
		MV->pushMatrix();
		MV->translate(x(0), x(1), x(2));
		MV->scale(r);
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		sphere->draw(prog);
		MV->popMatrix();
	}
}
