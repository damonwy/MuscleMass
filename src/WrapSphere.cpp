#include "WrapSphere.h"
#include <iostream>
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "WrapCylinder.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"
#include "Rigid.h"
#include "Particle.h"
#include "Vector.h"

using namespace std;
using namespace Eigen;

void WrapSphere::compute()
{
	Eigen::Vector3d OS = this->point_S - this->point_O;
	OS = OS / OS.norm();
	Eigen::Vector3d OP = this->point_P - this->point_O;
	OP = OP / OP.norm();
	Eigen::Vector3d N = OP.cross(OS);
	N = N / N.norm();

	this->M << OS.transpose(), N.cross(OS).transpose(), N.transpose();
	//std::cout << this->M << std::endl;

	Eigen::Vector3d p = this->M * (this->point_P - this->point_O);
	Eigen::Vector3d s = this->M * (this->point_S - this->point_O);

	double denom_q = p(0)*p(0) + p(1)*p(1);
	double denom_t = s(0)*s(0) + s(1)*s(1);
	double R = this->radius;

	if ((denom_q - R*R < 0.0) || (denom_t - R*R < 0.0))
	{
		this->status = inside_radius;

	}

	double root_q = sqrt(denom_q - R*R);
	double root_t = sqrt(denom_t - R*R);

	Eigen::Vector3d q(0.0, 0.0, 0.0);
	Eigen::Vector3d t(0.0, 0.0, 0.0);
	q(0) = (p(0) * R*R + R * p(1) * root_q) / denom_q;
	q(1) = (p(1) * R*R - R * p(0) * root_q) / denom_q;
	t(0) = (s(0) * R*R - R * s(1) * root_t) / denom_t;
	t(1) = (s(1) * R*R + R * s(0) * root_t) / denom_t;

	if (R * (q(0) * t(1) - q(1) * t(0)) > 0.0)
	{
		this->status = no_wrap;
	}

	this->status = wrap;
	this->point_q = q;
	this->point_t = t;

	//std::cout << q << std::endl << t << std::endl;

	Eigen::Vector3d Q = this->M.transpose() * q + this->point_O;
	Eigen::Vector3d T = this->M.transpose() * t + this->point_O;

	//  std::cout << Q.transpose() << std::endl << T.transpose() << std::endl;

	this->path_length = R * acos(1.0 - 0.5 *
		((Q(0) - T(0)) * (Q(0) - T(0))
			+ (Q(1) - T(1)) * (Q(1) - T(1))) / (R*R));
}

Eigen::MatrixXd WrapSphere::getPoints(int num_points)
{
	double theta_q = atan(this->point_q(1) / this->point_q(0));
	if (this->point_q(0) < 0.0)
		theta_q += PI;

	double theta_t = atan(this->point_t(1) / this->point_t(0));
	if (this->point_t(0) < 0.0)
		theta_t += PI;

	Eigen::MatrixXd points(3, num_points + 1);

	double theta_s, theta_e;

	if (theta_q < theta_t)
	{
		theta_s = theta_q; theta_e = theta_t;
	}
	else
	{
		theta_s = theta_t; theta_e = theta_q;
	}

	if (theta_e - theta_s > theta_s + 2 * PI - theta_e)
	{
		double tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2 * PI;
	}

	int col = 0;
	for (double i = theta_s; i <= theta_e + 0.001;
		i += (theta_e - theta_s) / num_points)
	{
		Eigen::Vector3d point = this->radius * this->M.transpose() *
			Eigen::Vector3d(cos(i), sin(i), 0.0) + this->point_O;
		points.col(col++) = point;
	}

	return points;
}

void WrapSphere::reset() {

}

void WrapSphere::step() {
	this->point_P = P->x;
	this->point_S = S->x;
	this->point_O = O->x;
	compute();
	if (this->status == wrap) {
		arc_points = getPoints(num_points);
		// for heon
		/*cout << "P" << point_P << endl;
		cout << "S" << point_S << endl;
		cout << "O" << point_O << endl;*/
	}
}

void WrapSphere::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const {

	prog->bind();
	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 0.0, 0.0).data());
	glUniform3fv(prog->getUniform("kdBack"), 1, Vector3f(1.0, 1.0, 0.0).data());

	// Draw P, Spoints
	this->P->draw(MV, prog);
	this->S->draw(MV, prog);
	
	// Draw Sphere
	this->O->draw(MV, prog);
	prog->unbind();

	// Draw wrapping
	prog2->bind();
	glUniformMatrix4fv(prog2->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog2->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	MV->pushMatrix();
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	glColor3f(0.0, 0.0, 0.0); // black
	glLineWidth(3);
	glBegin(GL_LINE_STRIP);
	glVertex3f(this->point_S(0), this->point_S(1), this->point_S(2));

	if (this->status == wrap) {
		for (int i = 0; i < this->arc_points.cols(); i++) {
			Vector3f p = this->arc_points.block<3, 1>(0, i).cast<float>();
			glVertex3f(p(0), p(1), p(2));
		}
	}

	glVertex3f(this->point_P(0), this->point_P(1), this->point_P(2));
	glEnd();
	MV->popMatrix();
	prog2->unbind();
}

void WrapSphere::setP(shared_ptr<Particle> _P) {
	this->P = _P;
}

void WrapSphere::setS(shared_ptr<Particle> _S) {
	this->S = _S;
}

void WrapSphere::setO(shared_ptr<Particle> _O) {
	this->O = _O;
	this->O->r = this->radius;
}

void WrapSphere::setParent(shared_ptr<Rigid> _parent) {
	this->parent = _parent;
}