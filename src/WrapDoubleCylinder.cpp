#include <iostream>
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "WrapDoubleCylinder.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"
#include "Rigid.h"
#include "Particle.h"
#include "Vector.h"

using namespace std;
using namespace Eigen;

void WrapDoubleCylinder::compute()
{
	// compute Matrix U and V
	Eigen::Vector3d OP = this->point_P - this->point_U;
	OP = OP / OP.norm();
	Eigen::Vector3d vec_Z_U = vec_z_U / vec_z_U.norm();
	Eigen::Vector3d vec_X_U = vec_Z_U.cross(OP);
	vec_X_U = vec_X_U / vec_X_U.norm();
	Eigen::Vector3d vec_Y_U = vec_Z_U.cross(vec_X_U);
	vec_Y_U = vec_Y_U / vec_Y_U.norm();

	Eigen::Vector3d OS = this->point_S - this->point_V;
	OS = OS / OS.norm();
	Eigen::Vector3d vec_Z_V = vec_z_V / vec_z_V.norm();
	Eigen::Vector3d vec_X_V = vec_Z_V.cross(OS);
	vec_X_V = vec_X_V / vec_X_V.norm();
	Eigen::Vector3d vec_Y_V = vec_Z_V.cross(vec_X_V);
	vec_Y_V = vec_Y_V / vec_Y_V.norm();

	this->M_U.resize(3, 3);
	this->M_U.row(0) = vec_X_U.transpose();
	this->M_U.row(1) = vec_Y_U.transpose();
	this->M_U.row(2) = vec_Z_U.transpose();

	this->M_V.resize(3, 3);
	this->M_V.row(0) = vec_X_V.transpose();
	this->M_V.row(1) = vec_Y_V.transpose();
	this->M_V.row(2) = vec_Z_V.transpose();

	// step 1: compute H and T
	Eigen::Vector3d pv = this->M_V * (this->point_P - this->point_V);
	Eigen::Vector3d sv = this->M_V * (this->point_S - this->point_V);

	double denom_h = pv(0)*pv(0) + pv(1)*pv(1);
	double denom_t = sv(0)*sv(0) + sv(1)*sv(1);
	double Rv = this->radius_V;

	double root_h = sqrt(denom_h - Rv*Rv);
	double root_t = sqrt(denom_t - Rv*Rv);

	Eigen::Vector3d h(0.0, 0.0, 0.0);
	Eigen::Vector3d t(0.0, 0.0, 0.0);
	h(0) = (pv(0) * Rv*Rv + Rv * pv(1) * root_h) / denom_h;
	h(1) = (pv(1) * Rv*Rv - Rv * pv(0) * root_h) / denom_h;
	t(0) = (sv(0) * Rv*Rv - Rv * sv(1) * root_t) / denom_t;
	t(1) = (sv(1) * Rv*Rv + Rv * sv(0) * root_t) / denom_t;

	this->status = wrap;

	std::complex<double> ht_i = 1.0 - 0.5 *
		((h(0) - t(0)) * (h(0) - t(0))
			+ (h(1) - t(1)) * (h(1) - t(1))) / (Rv*Rv);
	std::complex<double> ph_i = 1.0 - 0.5 *
		((pv(0) - h(0)) * (pv(0) - h(0))
			+ (pv(1) - h(1)) * (pv(1) - h(1))) / (Rv*Rv);
	std::complex<double> ts_i = 1.0 - 0.5 *
		((t(0) - sv(0)) * (t(0) - sv(0))
			+ (t(1) - sv(1)) * (t(1) - sv(1))) / (Rv*Rv);

	double ht_xy = abs(Rv * acos(ht_i));
	double ph_xy = abs(Rv * acos(ph_i));
	double ts_xy = abs(Rv * acos(ts_i));

	h(2) = pv(2) + (sv(2) - pv(2)) * ph_xy / (ph_xy + ht_xy + ts_xy);
	t(2) = sv(2) - (sv(2) - pv(2)) * ts_xy / (ph_xy + ht_xy + ts_xy);

	Eigen::Vector3d H = this->M_V.transpose() * h + this->point_V;
	Eigen::Vector3d T = this->M_V.transpose() * t + this->point_V;
	Eigen::Vector3d H0 = H;

	Eigen::Vector3d q(0.0, 0.0, 0.0);
	Eigen::Vector3d g(0.0, 0.0, 0.0);
	Eigen::Vector3d Q, G;

	double len = 0.0;

	for (int i = 0; i < 30; i++)
	{
		len = 0.0;

		// step 2: compute Q and G
		Eigen::Vector3d pu = this->M_U * (this->point_P - this->point_U);
		Eigen::Vector3d hu = this->M_U * (H - this->point_U);

		double denom_q = pu(0)*pu(0) + pu(1)*pu(1);
		double denom_g = hu(0)*hu(0) + hu(1)*hu(1);
		double Ru = -this->radius_U;

		double root_q = sqrt(denom_q - Ru*Ru);
		double root_g = sqrt(denom_g - Ru*Ru);

		q(0) = (pu(0) * Ru*Ru + Ru * pu(1) * root_q) / denom_q;
		q(1) = (pu(1) * Ru*Ru - Ru * pu(0) * root_q) / denom_q;
		g(0) = (hu(0) * Ru*Ru - Ru * hu(1) * root_g) / denom_g;
		g(1) = (hu(1) * Ru*Ru + Ru * hu(0) * root_g) / denom_g;

		std::complex<double> qg_i = 1.0 - 0.5 *
			((q(0) - g(0)) * (q(0) - g(0))
				+ (q(1) - g(1)) * (q(1) - g(1))) / (Ru*Ru);
		std::complex<double> pq_i = 1.0 - 0.5 *
			((pu(0) - q(0)) * (pu(0) - q(0))
				+ (pu(1) - q(1)) * (pu(1) - q(1))) / (Ru*Ru);
		std::complex<double> gh_i = 1.0 - 0.5 *
			((g(0) - hu(0)) * (g(0) - hu(0))
				+ (g(1) - hu(1)) * (g(1) - hu(1))) / (Ru*Ru);

		double qg_xy = abs(Rv * acos(qg_i));
		double pq_xy = abs(Rv * acos(pq_i));
		double gh_xy = abs(Rv * acos(gh_i));
		len += qg_xy;

		q(2) = pu(2) + (hu(2) - pu(2)) * pq_xy / (pq_xy + qg_xy + gh_xy);
		g(2) = hu(2) - (hu(2) - pu(2)) * gh_xy / (pq_xy + qg_xy + gh_xy);

		Q = this->M_U.transpose() * q + this->point_U;
		G = this->M_U.transpose() * g + this->point_U;

		// step 3: compute H based on G and T
		Eigen::Vector3d gv = this->M_V * (G - this->point_V);

		double denom_h = gv(0)*gv(0) + gv(1)*gv(1);
		double root_h = sqrt(denom_h - Rv*Rv);

		h = Eigen::Vector3d(0.0, 0.0, 0.0);
		h(0) = (gv(0) * Rv*Rv + Rv * gv(1) * root_h) / denom_h;
		h(1) = (gv(1) * Rv*Rv - Rv * gv(0) * root_h) / denom_h;

		std::complex<double> ht_i = 1.0 - 0.5 *
			((h(0) - t(0)) * (h(0) - t(0))
				+ (h(1) - t(1)) * (h(1) - t(1))) / (Rv*Rv);
		gh_i = 1.0 - 0.5 *
			((gv(0) - h(0)) * (gv(0) - h(0))
				+ (gv(1) - h(1)) * (gv(1) - h(1))) / (Rv*Rv);
		std::complex<double> ts_i = 1.0 - 0.5 *
			((t(0) - sv(0)) * (t(0) - sv(0))
				+ (t(1) - sv(1)) * (t(1) - sv(1))) / (Rv*Rv);

		double ht_xy = abs(Rv * acos(ht_i));
		gh_xy = abs(Rv * acos(gh_i));
		double ts_xy = abs(Rv * acos(ts_i));
		len += ht_xy;

		h(2) = gv(2) + (sv(2) - gv(2)) * gh_xy / (gh_xy + ht_xy + ts_xy);
		t(2) = sv(2) - (sv(2) - gv(2)) * ts_xy / (gh_xy + ht_xy + ts_xy);

		H = this->M_V.transpose() * h + this->point_V;
		T = this->M_V.transpose() * t + this->point_V;

		len += (G - H).norm();

		double dist = (H - H0).norm();
		if (dist == 0) break;

		H0 = H;
	}

	this->path_length = len;
	this->point_q = q;
	this->point_g = g;
	this->point_h = h;
	this->point_t = t;
	/*
	std::cout << Q.transpose() << std::endl << G.transpose() << std::endl
	<< H.transpose() << std::endl << T.transpose() << std::endl;
	*/
}

Eigen::MatrixXd WrapDoubleCylinder::getPoints(int num_points)
{
	double theta_q = atan(this->point_q(1) / this->point_q(0));
	if (this->point_q(0) < 0.0)
		theta_q += PI;

	double theta_g = atan(this->point_g(1) / this->point_g(0));
	if (this->point_g(0) < 0.0)
		theta_g += PI;

	double theta_h = atan(this->point_h(1) / this->point_h(0));
	if (this->point_h(0) < 0.0)
		theta_h += PI;

	double theta_t = atan(this->point_t(1) / this->point_t(0));
	if (this->point_t(0) < 0.0)
		theta_t += PI;

	Eigen::MatrixXd points(3, 3 * num_points + 1);

	double theta_s, theta_e, z_s, z_e;

	// q to g
	if (theta_q < theta_g)
	{
		theta_s = theta_q; theta_e = theta_g;
		z_s = this->point_q(2); z_e = this->point_g(2);
	}
	else
	{
		theta_s = theta_g; theta_e = theta_q;
		z_s = this->point_g(2); z_e = this->point_q(2);
	}

	if (theta_e - theta_s > theta_s + 2 * PI - theta_e)
	{
		double tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2 * PI;
		tmp = z_s; z_s = z_e; z_e = tmp;
	}

	int col = 0;
	double z_i = z_s, dz = (z_e - z_s) / num_points;
	for (double i = theta_s; i <= theta_e + 0.001;
		i += (theta_e - theta_s) / num_points)
	{	
		if (col == num_points + 1) {
			break;
		}
		Eigen::Vector3d point = this->M_U.transpose() *
			Eigen::Vector3d(this->radius_U * cos(i),
				this->radius_U * sin(i), z_i) +
			this->point_U;
		z_i += dz;
		points.col(col++) = point;
	}

	// g to h
	Eigen::Vector3d G = this->M_U.transpose() * this->point_g + this->point_U;
	Eigen::Vector3d H = this->M_V.transpose() * this->point_h + this->point_V;
	Eigen::Vector3d diff = H - G;

	for (int i = 1; i < num_points; i++)
	{
		Eigen::Vector3d point = G + diff / num_points * i;
		points.col(col++) = point;
	}

	// h to t
	if (theta_h < theta_t)
	{
		theta_s = theta_h; theta_e = theta_t;
		z_s = this->point_h(2); z_e = this->point_t(2);
	}
	else
	{
		theta_s = theta_t; theta_e = theta_h;
		z_s = this->point_t(2); z_e = this->point_h(2);
	}

	if (theta_e - theta_s > theta_s + 2 * PI - theta_e)
	{
		double tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2 * PI;
		tmp = z_s; z_s = z_e; z_e = tmp;
	}

	z_i = z_s;
	dz = (z_e - z_s) / num_points;

	int tt = col;
	col = 3 * num_points ;
	
	for (double i = theta_s; i <= theta_e + 0.001;
		i += (theta_e - theta_s) / num_points)
	{	
		if (tt == 3 * num_points + 1) {
			break;
		}
		Eigen::Vector3d point = this->M_V.transpose() *
			Eigen::Vector3d(this->radius_V * cos(i),
				this->radius_V * sin(i), z_i) +
			this->point_V;
		z_i += dz;
		points.col(col--) = point;
		tt++;
	}

	return points;
}

void WrapDoubleCylinder::step() {
	point_P = P->x;
	point_U = U->x;
	point_V = V->x;
	point_S = S->x;
	vec_z_U = z_U->dir;
	vec_z_V = z_V->dir;
	compute();

	if (this->status == wrap) {
		this->arc_points = getPoints(this->num_points);
		// for heon 
		/*cout << "P:" << endl<< P->x << endl << endl;
		cout << "U:" << endl << U->x << endl << endl;
		cout << "V:" << endl << V->x << endl << endl;
		cout << "S:" << endl << S->x << endl << endl;
		cout << "zu:" << endl << z_U->dir << endl << endl;
		cout << "zv:" << endl << z_V->dir << endl << endl;*/
	}
}

void WrapDoubleCylinder::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> prog2, shared_ptr<MatrixStack> P) const {
	// Draw double cylinder
	prog->bind();

	if (cylinder_shape) {
		// Draw U
		glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 0.0, 0.0).data());
		glUniform3fv(prog->getUniform("kdBack"), 1, Vector3f(1.0, 1.0, 0.0).data());
		MV->pushMatrix();
		Vector3d x = getp_U();
		MV->translate(x(0), x(1), x(2));
		// Decompose R into 3 Euler angles
		Matrix3d R = getR_U();
		double theta_x = atan2(R(2, 1), R(2, 2));
		double theta_y = atan2(-R(2, 0), sqrt(pow(R(2, 1), 2) + pow(R(2, 2), 2)));
		double theta_z = atan2(R(1, 0), R(0, 0));
		MV->rotate(theta_z, 0.0f, 0.0f, 1.0f);
		MV->rotate(theta_y, 0.0f, 1.0f, 0.0f);
		MV->rotate(theta_x, 1.0f, 0.0f, 0.0f);
		MV->scale(this->radius_U);
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		cylinder_shape->draw(prog);
		MV->popMatrix();
		
		// Draw V
		MV->pushMatrix();
		x = this->getp_V();
		MV->translate(x(0), x(1), x(2));
		// Decompose R into 3 Euler angles
		R = this->getR_V();
		theta_x = atan2(R(2, 1), R(2, 2));
		theta_y = atan2(-R(2, 0), sqrt(pow(R(2, 1), 2) + pow(R(2, 2), 2)));
		theta_z = atan2(R(1, 0), R(0, 0));
		MV->rotate(theta_z, 0.0f, 0.0f, 1.0f);
		MV->rotate(theta_y, 0.0f, 1.0f, 0.0f);
		MV->rotate(theta_x, 1.0f, 0.0f, 0.0f);
		MV->scale(this->radius_V);
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		cylinder_shape->draw(prog);
		MV->popMatrix();	
	}

	// Draw P, S, U, V points
	this->P->draw(MV, prog);
	this->S->draw(MV, prog);
	this->U->draw(MV, prog);
	this->V->draw(MV, prog);

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
	glVertex3f(this->point_P(0), this->point_P(1), this->point_P(2));
	//todo
	if (this->status == wrap) {
		for (int i = 0; i < this->arc_points.cols(); i++) {
			Vector3f p = this->arc_points.block<3, 1>(0, i).cast<float>();
			glVertex3f(p(0), p(1), p(2));
		}
	}
	glVertex3f(this->point_S(0), this->point_S(1), this->point_S(2));
	glEnd();

	// Draw z axis
	z_U->draw(MV, P, prog2);
	z_V->draw(MV, P, prog2);

	MV->popMatrix();
	prog2->unbind();
}

// get
Vector3d WrapDoubleCylinder::getp_U() const {
	return this->E_W_U.block<3, 1>(0, 3);
}

Matrix3d WrapDoubleCylinder::getR_U() const {
	return this->E_W_U.block<3, 3>(0, 0);
}

Vector3d WrapDoubleCylinder::getp_V() const {
	return this->E_W_V.block<3, 1>(0, 3);
}

Matrix3d WrapDoubleCylinder::getR_V() const {
	return this->E_W_V.block<3, 3>(0, 0);
}

// set

void WrapDoubleCylinder::setE_U(Matrix4d E) {
	this->E_W_U = E;
}

void WrapDoubleCylinder::setE_V(Matrix4d E) {
	this->E_W_V = E;
}

void WrapDoubleCylinder::setP(shared_ptr<Particle> _P) {
	this->P = _P;
	this->point_P = P->x;
}

void WrapDoubleCylinder::setS(shared_ptr<Particle> _S) {
	this->S = _S;
	this->point_S = S->x;
}

void WrapDoubleCylinder::setU(shared_ptr<Particle> _U) {
	this->U = _U;
	this->point_U = U->x;
}

void WrapDoubleCylinder::setV(shared_ptr<Particle> _V) {
	this->V = _V;
	this->point_V = V->x;
}

void WrapDoubleCylinder::setZ_U(shared_ptr<Vector> _z_U) {
	this->z_U = _z_U;
	this->vec_z_U = this->z_U->dir;
}

void WrapDoubleCylinder::setZ_V(shared_ptr<Vector> _z_V) {
	this->z_V = _z_V;
	this->vec_z_V = this->z_V->dir;
}

void WrapDoubleCylinder::setNumPoints(int _num_points) {
	this->num_points = _num_points;
}

void WrapDoubleCylinder::setParent_U(std::shared_ptr<Rigid> _parent_U) {
	this->parent_U = _parent_U;
}

void WrapDoubleCylinder::setParent_V(std::shared_ptr<Rigid> _parent_V) {
	this->parent_V = _parent_V;
}

void WrapDoubleCylinder::setE_P_U(Eigen::Matrix4d E) {
	this->E_P_U = E;
}

void WrapDoubleCylinder::setE_P_V(Eigen::Matrix4d E) {
	this->E_P_V = E;
}


Matrix4d WrapDoubleCylinder::getE_U() const {
	return this->E_W_U;
}

Matrix4d WrapDoubleCylinder::getE_V() const {
	return this->E_W_V;
}

Matrix4d WrapDoubleCylinder::getE_P_U() const {
	return this->E_P_U;
}

Matrix4d WrapDoubleCylinder::getE_P_V() const {
	return this->E_P_V;
}

shared_ptr<Rigid> WrapDoubleCylinder::getParent_U() const {
	return this->parent_U;
}

shared_ptr<Rigid> WrapDoubleCylinder::getParent_V() const {
	return this->parent_V;
}