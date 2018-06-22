#include "Spring.h"

#include <iostream>
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Particle.h"
#include "Joint.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"
#include "Rigid.h"
#include <unsupported/Eigen/MatrixFunctions> // TODO: avoid using this later, write a func instead

using namespace std;
using namespace Eigen;

Spring::Spring(shared_ptr<Particle> p0, shared_ptr<Particle> p1, double _mass, int num_samples, Vector3d _grav, double _epsilon, bool _isReduced, double _stiffness) :
	E(_stiffness), mass(_mass), grav(_grav), epsilon(_epsilon), isReduced(_isReduced)
{
	assert(p0);
	assert(p1);
	assert(p0 != p1);
	this->p0 = p0;
	this->p1 = p1;
	Vector3d x0 = p0->x;
	Vector3d x1 = p1->x;
	Vector3d dx = x1 - x0;
	this->L = dx.norm();
	this->l = L;
	this->phi_box.setZero();

	assert(L > 0.0);
	box_id(0) = p0->getParent()->getIndex();
	box_id(1) = p1->getParent()->getIndex();

	Matrix3x12d J;
	J.setZero();
	
	double dm = this->mass / num_samples;

	for (int i = 0; i < num_samples; ++i) {
		auto sample = make_shared<Particle>();
		double s = i / double(num_samples - 1);
		//cout << s << endl;
		sample->s = s;
		sample->x = (1 - s)*p0->x + s*p1->x;
		sample->m = dm;
		if (!isReduced) {
			sample->setJacobianMatrix(J);
		}
		samples.push_back(sample);
	}
}

void Spring::step(vector<shared_ptr<Joint>> joints) {
	computeLength();
	updateSamplesPosition();
	updateSamplesJacobian(joints);
	computeEnergy();
}

double Spring::computeLength() {
	this->l = (p0->x - p1->x).norm();
	return this->l;
}

void Spring::updateSamplesPosition() {
	for (int i = 0; i < (int)this->samples.size(); ++i) {
		auto sample = samples[i];
		double s = sample->s;
		sample->x = (1 - s)*p0->x + s*p1->x;
	}
}

void Spring::updateSamplesJacobian(vector<shared_ptr<Joint>> joints) {
	

	if (isReduced) {
		//// Reduced Coordinate
		//int num_joints = (int)joints.size();

		//// Resize Jacobian
		//MatrixXd J(3, num_joints);
		//J.setZero();
		//for (int i = 0; i < (int)samples.size(); ++i) {		
		//	samples[i]->setJacobianMatrix(J);
		//}

		//VectorXd pert(num_joints);
		//VectorXd thetalist(num_joints);
		//thetadotlist.resize(num_joints);
		//thetalist = Joint::getThetaVector(joints);

		//// Used to compute the velocity and kinetic energy of samples
		//thetadotlist = Joint::getThetadotVector(joints);

		//// For each component of thetalist add a relative small perturbation
		//for (int ii = 0; ii < num_joints; ++ii) {
		//	// Change ii-th component
		//	pert = thetalist;
		//	pert(ii) += epsilon;

		//	// Compute new configuration
		//	/*for (int jj = 0; jj < num_joints; ++jj) {
		//		auto joint = joints[jj];
		//		auto box = joint->getChild();
		//		box->setJointAngle(pert(jj), false);
		//		box->updateTempPoints();
		//	}*/
		//	//joints[ii]->getChild()->setJointAngle(pert(ii), false);

		//	// Fill in the ii-th column of Jacobian
		//	for (int isample = 0; isample < (int)samples.size(); ++isample) {
		//		auto sample = samples[isample];
		//		double s = sample->s;

		//		Vector3d p_nopert = sample->x;
		//		Vector3d p_pert = (1 - s) * p0->getTempPos() + s * p1->getTempPos();

		//		Vector3d diff = (p_pert - p_nopert) / epsilon;
		//		
		//		// To test if the Jacobian is correct
		//		int nth_joint = 0;
		//		if (nth_joint == 1) {
		//			auto debug = make_shared<Particle>();
		//			debug->x = p_pert;
		//			this->debug_points.push_back(debug);
		//		}
		//					
		//		sample->setJacobianMatrixCol(diff, ii);
		//	}
		//}

		MatrixXd J(3, 2);
		J.setZero();
		for (int i = 0; i < (int)samples.size(); ++i) {		
			samples[i]->setJacobianMatrix(J);
		}

		auto b0 = p0->getParent();
		auto b1 = p1->getParent();

		Vector2d theta;
		theta(0) = b0->getAngle();
		theta(1) = b1->getAngle();

		// Used to compute velocity of samples and energy
		thetadot(0) = b0->getThetadot();
		thetadot(1) = b1->getThetadot();
		Vector2d pert;
		// For two angles add a small perturbation
		pert = theta;
		pert(0) += epsilon;

		b0->setSingleJointAngle(pert(0));
		b1->setJointAngle(pert(1), false);
		b0->updateTempPoints();
		b1->updateTempPoints();

		for (int isample = 0; isample < (int)samples.size(); ++isample) {
			auto sample = samples[isample];
			double s = sample->s;

			Vector3d p_nopert = sample->x;
			Vector3d p_pert = (1 - s) * p0->getTempPos() + s * p1->getTempPos();
			//Vector3d p_pert = (1 - s) * p0->getTempPos() + s * p1->x;
			Vector3d diff = (p_pert - sample->x) / epsilon;
			sample->setJacobianMatrixCol(diff, 0);
		}

		pert = theta;
		pert(1) += epsilon;

		b1->setSingleJointAngle(pert(1));
		b1->updateTempPoints();

		for (int isample = 0; isample < (int)samples.size(); ++isample) {
			auto sample = samples[isample];
			double s = sample->s;

			Vector3d p_nopert = sample->x;
			Vector3d p_pert = (1 - s) * p0->x + s * p1->getTempPos();

			Vector3d diff = (p_pert - sample->x) / epsilon;
			sample->setJacobianMatrixCol(diff, 1);
		}
	}
	else {
		// Maximal Coordinate
		auto b0 = p0->getParent();
		auto b1 = p1->getParent();

		// Used to compute velocity of samples and energy
		phi_box.segment<6>(0) = b0->getTwist();
		phi_box.segment<6>(6) = b1->getTwist();

		Vector6d pert;

		// For each component of phi(i = 0, 1, 2..,11) add a relative small perturbation
		for (int ii = 0; ii < 6; ++ii) {
			// Change ii-th component
			pert.setZero();
			pert(ii) += epsilon;

			// Compute new configuration
			Matrix4d E_pert = b0->getE() * Rigid::bracket6(pert).exp();
			b0->setEtemp(E_pert);
			b0->updateTempPoints();

			for (int isample = 0; isample < (int)samples.size(); ++isample) {
				auto sample = samples[isample];
				double s = sample->s;

				Vector3d p_nopert = sample->x;
				Vector3d p_pert = (1 - s) * p0->getTempPos() + s * p1->x;

				Vector3d diff = (p_pert - sample->x) / epsilon;
				sample->setJacobianMatrixCol(diff, ii);
			}
		}

		for (int ii = 0; ii < 6; ++ii) {
			// Change ii-th component
			pert.setZero();
			pert(ii) += epsilon;

			// Compute new configuration
			Matrix4d E_pert = b1->getE() * Rigid::bracket6(pert).exp();
			b1->setEtemp(E_pert);
			b1->updateTempPoints();

			// Fill in the ii-th column of Jacobian
			for (int isample = 0; isample < (int)samples.size(); ++isample) {
				auto sample = samples[isample];
				double s = sample->s;

				Vector3d p_nopert = sample->x;
				Vector3d p_pert = (1 - s) * p0->x + s * p1->getTempPos();
				
				/*if (ii == 0) {
					auto debug = make_shared<Particle>();
					debug->x = p_pert;
					debug_points.push_back(debug);
				}*/
				

				Vector3d diff = (p_pert - sample->x) / epsilon;
				sample->setJacobianMatrixCol(diff, ii + 6);
			}
		}
	}
}

void Spring::computeEnergy() {
	this->V = 0.0;
	this->K = 0.0;
	double V_ii, K_ii;

	for (int isample = 0; isample < (int)samples.size(); ++isample) {
		auto sample = samples[isample];

		// Update the energy	
		V_ii = sample->computePotentialEnergy(grav);
		if (isReduced) {
			//K_ii = sample->computeKineticEnergy(thetadotlist);
			K_ii = sample->computeKineticEnergy(thetadot);
		}
		else {
			K_ii = sample->computeKineticEnergy(phi_box);
		}		

		this->V += V_ii;
		this->K += K_ii;
	}
}

MatrixXd Spring::computeMassMatrix(vector<shared_ptr<Spring> > springs, int num_boxes, bool isReduced) {
	if (isReduced) {
		int n = num_boxes - 1;
		MatrixXd M_s(n, n);
		M_s.setZero();

		for (int i = 0; i < (int)springs.size(); ++i) {
			auto spring = springs[i];
			auto samples = spring->getSamples();

			// Sum up the inertia matrix of all the sample points
			for (int isample = 0; isample < (int)samples.size(); ++isample) {
				auto sample = samples[isample];
				MatrixXd J_ii = sample->getJacobianMatrix();
				MatrixXd M_ii = J_ii.transpose() * J_ii * sample->m;
			
				// Fill in Mass matrix
				M_s += M_ii;
			}
		}
		return M_s;
	}
	else {
		int n = 6 * num_boxes;
		MatrixXd M_s(n, n);
		M_s.setZero();

		for (int i = 0; i < (int)springs.size(); ++i) {
			auto spring = springs[i];
			auto samples = spring->getSamples();
			Vector2d box_id = spring->getBoxID();

			// Sum up the inertia matrix of all the sample points
			for (int isample = 0; isample < (int)samples.size(); ++isample) {
				auto sample = samples[isample];
				Matrix3x12d J_ii = sample->getJacobianMatrix();
				Matrix12d M_ii = J_ii.transpose() * J_ii * sample->m;

				// Fill in Mass matrix
				M_s.block<6, 6>(6 * box_id(0), 6 * box_id(0)) += M_ii.block<6, 6>(0, 0);
				M_s.block<6, 6>(6 * box_id(0), 6 * box_id(1)) += M_ii.block<6, 6>(0, 6);
				M_s.block<6, 6>(6 * box_id(1), 6 * box_id(0)) += M_ii.block<6, 6>(6, 0);
				M_s.block<6, 6>(6 * box_id(1), 6 * box_id(1)) += M_ii.block<6, 6>(6, 6);
			}
		}
		return M_s;
	}	
}

VectorXd Spring::computeGravity(vector<shared_ptr<Spring> > springs, int num_boxes, bool isReduced) {
	if (isReduced) {
		int n = num_boxes - 1;
		VectorXd b_s(n);
		b_s.setZero();
		for (int i = 0; i < (int)springs.size(); ++i) {
			auto spring = springs[i];
			auto samples = spring->getSamples();
			
			for (int isample = 0; isample < (int)samples.size(); ++isample) {
				auto sample = samples[isample];
				MatrixXd J_ii = sample->getJacobianMatrix();

				// Add gravtiy force for each sample to b vector
				VectorXd f_ii = J_ii.transpose() * sample->m * spring->grav;
				b_s += f_ii;
			}
		}
		return b_s;
	}
	else {
		int n = 6 * num_boxes;
		VectorXd b_s(n);
		b_s.setZero();

		for (int i = 0; i < (int)springs.size(); ++i) {
			auto spring = springs[i];
			auto samples = spring->getSamples();
			Vector2d box_id = spring->getBoxID();

			for (int isample = 0; isample < (int)samples.size(); ++isample) {
				auto sample = samples[isample];
				Matrix3x12d J_ii = sample->getJacobianMatrix();

				// Add gravtiy force for each sample to b vector
				Vector12d f_ii = J_ii.transpose() * sample->m * spring->grav;

				b_s.segment<6>(6 * box_id(0)) += f_ii.segment<6>(0);
				b_s.segment<6>(6 * box_id(1)) += f_ii.segment<6>(6);
			}
		}
		return b_s;
	}
 }


void Spring::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const {

	prog->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	MV->pushMatrix();
	glUniform3f(prog->getUniform("lightPos1"), 1.0, 1.0, 1.0);
	glUniform1f(prog->getUniform("intensity_1"), 0.8);
	glUniform3f(prog->getUniform("lightPos2"), -1.0, 1.0, 1.0);
	glUniform1f(prog->getUniform("intensity_2"), 0.2);
	glUniform1f(prog->getUniform("s"), 200);
	glUniform3f(prog->getUniform("ka"), 0.9, 0.5, 0.2);
	glUniform3f(prog->getUniform("kd"), 0, 0, 1);
	glUniform3f(prog->getUniform("ks"), 0, 1.0, 0);

	// Draw P, Spoints
	this->p0->draw(MV, prog);
	this->p1->draw(MV, prog);

	prog->unbind();

	// Draw wrapping
	prog2->bind();
	glUniformMatrix4fv(prog2->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog2->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	MV->pushMatrix();
	glColor3f(0.1, 0.5, 0.6); 
	glLineWidth(4);
	Vector3d p = p0->x;
	Vector3d s = p1->x;

	glBegin(GL_LINE_STRIP);
	glVertex3f(p(0), p(1), p(2));
	glVertex3f(s(0), s(1), s(2));
	glEnd();

	glPointSize(3.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < (int)debug_points.size(); ++i) {
		Vector3f p = debug_points[i]->x.cast<float>();
		glVertex3f(p(0), p(1), p(2));
	}
	glEnd();

	MV->popMatrix();
	prog2->unbind();
}

Spring::~Spring()
{
	
}
