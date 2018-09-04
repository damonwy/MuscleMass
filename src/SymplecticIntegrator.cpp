#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#define _USE_MATH_DEFINES
#include <math.h>

#include "SymplecticIntegrator.h"

#include "Rigid.h"
#include "Joint.h"
#include "MatlabDebug.h"
#include "Spring.h"
#include "Particle.h"
#include "QuadProgMosek.h"
#include "Scene.h"
#include "Program.h"
#include "MatrixStack.h"
#include "MLShapeInfo.h"
#include "MLError.h"

#include <iostream>
#include <iomanip>

# include <cstdlib>

# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;
using namespace Eigen;
#include <unsupported/Eigen/MatrixFunctions> // TODO: avoid using this later, write a func instead
double inf = numeric_limits<double>::infinity();

SymplecticIntegrator::SymplecticIntegrator(vector< shared_ptr<Rigid> > _boxes, vector< shared_ptr<Joint>> _joints, vector< shared_ptr<Spring> > _springs, bool _isReduced, int _num_samples, Vector3d _grav, double _epsilon):
		num_joints(_boxes.size() - 1),
		num_samples(_num_samples),
		boxes(_boxes),
		joints(_joints),
		springs(_springs),
		isReduced(_isReduced),
		epsilon(_epsilon),
		grav(_grav)
{

	if (isReduced) {
		m = 6 * (int)boxes.size();
		n = 6 + 1 * ((int)boxes.size() - 1);

		M.resize(m, m);	
		M.setZero();

		J.resize(m, n);
		J.setZero();

		Jdot.resize(m, n);
		Jdot.setZero();

		f.resize(m);
		f_b.resize(m);
		f_c.resize(m);
		ftest.resize(m);
		
		f.setZero();
		ftest.setZero();
		f_b.setZero();
		f_c.setZero();
	}
	else {
		m = 6 * (int)boxes.size();
		n = 6 * (int)boxes.size() + 6 + 5 * ((int)boxes.size() - 1);
		M.resize(m, m);
		M.setZero();
	}

	A.resize(n, n);
	x.resize(n);
	b.resize(n);

	A.setZero();
	x.setZero();
	b.setZero();
}

MatrixXd SymplecticIntegrator::getGlobalJacobian(VectorXd thetalist) {
	// Transfer reduced coords to maximal coords
	// Assume the first box is fixed and do not include the first Identity matrix

	int n = thetalist.size();
	MatrixXd J(6 * n, n);
	J.setZero();

	// Rotate about Z axis
	Vector6d z;
	z.setZero();
	z(2) = 1.0;

	for (int i = 0; i < n; ++i) {
		auto joint = joints[i];
		double target_theta = thetalist(i);

		Matrix4d R;
		R.setIdentity();
		R.block<2, 2>(0, 0) << cos(target_theta), -sin(target_theta),
			sin(target_theta), cos(target_theta);

		Matrix4d E_C_J_new = joint->getE_C_J_0() * R.inverse();
		Matrix6d Ad_C_J = Rigid::adjoint(E_C_J_new);
		Vector6d Adz_C_J =  Ad_C_J * z;

		if (i != 0) {
			// If it is not the first joint, need to compute the off-diag entries
			Matrix6d Ad_J_P = Rigid::adjoint(joint->getE_P_J().inverse());
			Matrix6d Ad_C_P = Ad_C_J * Ad_J_P;
			J.block(6 * i, 0, 6, n-1) = Ad_C_P * J.block(6 * (i - 1), 0, 6, n - 1);	
		}
		// Compute the diag entries
		J.block<6, 1>(6 * i, i) = Adz_C_J;
	}
	return J;
}


MatrixXd SymplecticIntegrator::getJ_twist_thetadot() {
	J.setZero();
	int j = 6;

	// Rotate about Z axis
	Vector6d z;
	z.setZero();
	z(2) = 1.0;

	for (int i = 0; i < (int)boxes.size(); i++) {
		auto box = boxes[i];

		M.block<6, 6>(6 * i, 6 * i) = box->getMassMatrix();

		if (i == 0) {
			Matrix6d I;
			I.setIdentity();
			J.block<6, 6>(0, 0) = I;
		}
		else if (i == 1) {
			auto joint = box->getJoint();
			Matrix6d Ad_C_J = Rigid::adjoint(joint->getE_C_J());
			Matrix6d Ad_J_P = Rigid::adjoint(joint->getE_P_J().inverse());
			Matrix6d Ad_C_P = Ad_C_J * Ad_J_P;
			Vector6d Adz_C_J = Ad_C_J * z;
			
			J.block<6, 6>(j, 0) = Ad_C_P;
			J.block<6, 1>(j, 5 + i) = Adz_C_J;
			j = j + 6;
		}
		else {
			auto joint = box->getJoint();
			Matrix6d Ad_C_J = Rigid::adjoint(joint->getE_C_J());
			Matrix6d Ad_J_P = Rigid::adjoint(joint->getE_P_J().inverse());
			Vector6d Adz_C_J = Ad_C_J * z;
			Matrix6d Ad_C_P = Ad_C_J * Ad_J_P;

			J.block(j, 0, 6, 5 + i) = Ad_C_P * J.block(j - 6, 0, 6, 5 + i);
			J.block<6, 1>(j, 5 + i) = Adz_C_J;
			j = j + 6;
		}
	}
	return J;
}

Eigen::MatrixXd SymplecticIntegrator::getJdot(VectorXd thetadotlist) {
	Jdot.setZero();
	MatrixXd JJ = J.block(0, n - num_joints, m, num_joints);
	VectorXd Phi = JJ * thetadotlist;
	cout << "Phi:" << Phi << endl;
	for (int i = 1; i < (int)boxes.size(); i++) {
		cout << "i: " << i << endl;
		auto box = boxes[i];
		int ii = (i ) * 6;
		int p = i - 1;
		int pp = (p) * 6;

		auto joint = box->getJoint();
		Matrix4d E_C_J = joint->getE_C_J();
		//cout << E_C_J << endl;
		Matrix4d E_J_P = joint->getE_P_J().inverse();
		Matrix4d E_C_P = E_C_J * E_J_P;

		Matrix6d Ad_C_P = Rigid::adjoint(E_C_P);
		Matrix6d Ad_C_0 = Rigid::adjoint(box->getE().inverse());
		Matrix6d Ad_0_P = Rigid::adjoint(box->getParent()->getE());

		Matrix6d Addot_0_C = Rigid::dAddt(box->getE(), Phi.segment<6>(i*6));

		Matrix6d Addot_0_P = Rigid::dAddt(box->getParent()->getE(), Phi.segment<6>(p*6));
		Matrix6d Addot_C_P = -Ad_C_0 * (Addot_0_C * Ad_C_0 * Ad_0_P - Addot_0_P);

		for (int j = 0; j < i ; j++) {
			Jdot.block<6, 1>(ii, j + 6) = Ad_C_P * Jdot.block<6, 1>(pp, j+6) + Addot_C_P * J.block<6, 1>(pp, j + 6);
			cout << Ad_C_P * Jdot.block<6, 1>(pp, j + 6) + Addot_C_P * J.block<6, 1>(pp, j + 6) << endl;

		}

	}
	return Jdot;

}

void SymplecticIntegrator::step(double h) {
	A.setZero();
	x.setZero();
	b.setZero();
	M.setZero();

	// Test shapeinfo
	//std::shared_ptr<MLFunctionTestShapeInfo> test0 = make_shared<MLFunctionTestShapeInfo>(1.5);
	//std::shared_ptr<MLFunctionTestShapeInfo> test1 = make_shared<MLFunctionTestShapeInfo>(3.5);
	//double result;
	//test0->computeDifference(test1, &result);
	//cout << "result" << result << endl;


	// Solve linear system
	if (isReduced) {
		M.setZero();
		J.setZero();
		f.setZero();
		f_b.setZero();
		f_c.setZero();
		ftest.setZero();

		// Input
		VectorXd thetadotlist = Joint::getThetadotVector(joints);
		VectorXd thetalist = Joint::getThetaVector(joints);

		// Matlab Solution
		double c1 = cos( M_PI / 2.0 + thetalist(0));
		double c2 = cos(thetalist(1));
		double s1 = sin( M_PI / 2.0 + thetalist(0));
		double s2 = sin(thetalist(1));
		double c12 = cos(thetalist(0) + M_PI / 2.0 + thetalist(1));
		double s12 = sin(thetalist(0) + M_PI / 2.0 + thetalist(1));
		double mu1 = 1.0;
		double mu2 = 1.0;
		double mum = 1.0;
		double l = 1.0;
		double r = 0.1;
		double bigm = 4.2;
		double bigl = 4.0;

		Matrix2d I1;
		I1.setZero();
		I1(0, 0) = mu1 / 3.0 * l * l;

		Matrix2d I2;
		I2.setZero();
		I2 << 1 + 3 * l * (l + c2), 1 + 1.5 * l * c2, 1 + 1.5 * l * c2, 1.0;
		I2 *= mu2 / 3.0;

		Matrix2d Im;
		Im.setZero();
		Im << l * l + r * r + 2 * l * r * c2, r * r + l * r * c2, r * r + l * r * c2, r * r;
		Im *= mum / 3.0;

		Matrix2d Iall = bigm * bigl * bigl * (I1 + I2 + Im);

		Matrix2d temp0, temp1;
		temp0 << 1.0, 0.5, 0.5, 0.0;
		temp1 << 2.0, 1.0, 1.0, 0.0;

		Matrix2d dIdtheta2 = -s2 * (mu2 * l * temp0 + mum / 3.0 * l * r * temp1);
		//dIdtheta2 = -s2 * (mu2 * l * temp0 );

		Vector2d fk = -thetadotlist(1) * dIdtheta2 * thetadotlist;
		fk(1) = fk(1) + 0.5 * thetadotlist.transpose() * dIdtheta2 * thetadotlist;
		
		Vector2d fv1, fv2;
		fv1 << -0.5 * mu1 * 9.81 * l * s1, 0.0;
		
		fv2 << 2 * l * s1 + s12, s12;
		fv2 *= -0.5 * mu2 * 9.81;
		
		Vector2d fvm;
		fvm << l * s1 + r * s12, r *s12;
		fvm *= -0.5 * mum * 9.81;
		
		Vector2d fk_matlab = bigm * bigl * bigl * fk;
	

		Vector2d fall = fk_matlab + bigm * bigl * (fv1 + fv2 + fvm);
		Vector2d thetaddot_matlab = Iall.ldlt().solve(fall);
		// End of Matlab solution

		for (int i = 0; i < (int)boxes.size(); i++) {
			auto box = boxes[i];
			//f.segment<6>(6 * i) = box->getMassMatrix() * box->getTwist() + h * box->getForce();
			f.segment<6>(6 * i) = box->getForce();
			f_b.segment<6>(6 * i) = box->getBodyForce();
			f_c.segment<6>(6 * i) = box->getCoriolisForce();
			Vector3d omega = box->getTwist().segment<3>(0);
			Matrix6d twist_bracket;
			twist_bracket.setZero();
			twist_bracket.block<3, 3>(0, 0) = Rigid::bracket3(omega.segment<3>(0));
			twist_bracket.block<3, 3>(3, 3) = Rigid::bracket3(omega.segment<3>(0));
			VectorXd tettt = twist_bracket * box->getMassMatrix() *box->getTwist();
			//cout << "cof" << tettt << endl;

			ftest.segment<6>(6 * i) = box->getMassMatrix() * box->getTwist() + h * box->getForce();
			//cout << box->getTwist() << endl;
		}
		

		J = getJ_twist_thetadot();
		// TODO GET JDOT
		Jdot = getJdot(thetadotlist);
		cout << "Jdot" << endl << Jdot << endl;

		//cout << "Jcomp" << J << endl;
		VectorXd testtest(8);
		testtest.setZero();
		testtest.segment<2>(6) = thetadotlist;
		//cout << "computed phi" << J * testtest << endl;

		MatrixXd Jtest = getGlobalJacobian(thetalist);

		A = J.transpose() * M * J;
		b = J.transpose() * f;
		
		VectorXd b_b = J.transpose() * f_b;
		VectorXd b_c = J.transpose() * f_c;
		//cout << "b_b" << b_b << endl;
		//cout << "b_c" << b_c << endl;
		//cout << "fk" << bigm * bigl * bigl * fk << endl;
		//cout << "fv" << bigm * bigl * (fv1 + fv2) << endl;

		//cout << "fc" << f_c << endl;
		
		// Compute the inertia matrix of spring using finite difference 
		MatrixXd M_s = Spring::computeMassMatrix(springs, (int)boxes.size(), isReduced);

		x.setZero();
		MatrixXd JJ = J.block(0, n - num_joints, m, num_joints);
		//JJ.block<12, 2>(6, 0) = Jtest;
		A = JJ.transpose() * M * JJ;
		VectorXd bnew = JJ.transpose() * f;
		VectorXd bbb = JJ.block<12, 2>(6, 0).transpose() * f_c.segment<12>(6);
		//cout << "bbbfc" << bbb << endl;

		VectorXd bb = J.transpose() * ftest;
		//cout << A - bigm * bigl * bigl * (I1 + I2) << endl;

		A += M_s;
		
		//cout << "test Mass Matrix: " << endl;
		//cout << A - Iall << endl;
		//cout << "A" << A << endl;
		//cout << "Iall:" << Iall << endl;

		
		VectorXd b_s = Spring::computeGravity(springs, (int)boxes.size(), isReduced);
		VectorXd xtest = x;
		//x.segment(6, num_joints) = A.ldlt().solve(b.segment(6, num_joints) + b_s * h + M_s * thetadotlist);	// thetadot
		x.segment(6, num_joints) = A.ldlt().solve(b.segment(6, num_joints) + b_s);	// thetaddot
		xtest.segment(6, num_joints) = A.ldlt().solve(bb.segment(6, num_joints) + b_s * h + M_s * thetadotlist); // thetadot
		//cout << "f" << f << endl;
		//cout << "b" << b.segment(6, num_joints) << endl;
		//cout << "fall" << fall << endl;
		//cout << "Idiff:" << (Iall - A).norm() << endl;

		//cout << b.segment(6, num_joints) << endl;
		//cout << endl << fall - b_s << endl;
		//cout << "fdiff" << endl << (b.segment(6, num_joints) + b_s - fall).norm() << endl;
		//cout << "b_ccc" << b.segment(6, num_joints) + b_s - fall - b_c.segment(6, num_joints) + bigm * bigl * bigl * fk << endl;
		//cout << "fvmdiff" << endl << (b_s - 4.2 * 4.0 * (fvm)).norm() << endl;
		//cout << "accel diff:" << endl << (thetaddot_matlab - x.segment(6, num_joints)).norm() << endl;
		
		Vector2d newthetadotlist = thetadotlist + h * x.segment(6, num_joints);
		newthetadotlist = thetadotlist + h * thetaddot_matlab;

		//cout << thetaddot_matlab << endl;
		//cout << thetadotlist << endl;

		//cout << "velocity diff:" << (xtest.segment(6, num_joints) - newthetadotlist).norm() << endl;
		//cout << "thetadot diff" << endl << (thetaddot_matlab - (x.segment(6, num_joints) - thetadotlist) / h).norm() << endl;

		VectorXd phi = JJ * newthetadotlist;
		//cout << "JJ" << JJ << endl;

		Vector2d thetalistnew;

		// For QP
		VectorXd xl, xu; // lower, upper bound 
		xl.resize(num_joints); 
		xu.resize(num_joints);
		bool isQP = false;

		for (int i = 1; i < (int)boxes.size(); i++) {
			auto box = boxes[i];
			
			// Check if the new joint angle is out of range
			double dtheta = h * x(5 + i);
			double theta_old = box->getJoint()->getTheta();
			double theta_new = theta_old + dtheta;

			double min_theta = box->getJoint()->getMinTheta();
			double max_theta = box->getJoint()->getMaxTheta();
			
			// Two cases need to solve QP 
			if (theta_new > max_theta) {
				// thetadot <= 0.0 
				xl(i - 1) = -inf;
				xu(i - 1) = 0.0;
				isQP = true;
			}
			else if (theta_new < min_theta) {
				// thetadot >= 0.0
				xl(i - 1) = 0.0;
				xu(i - 1) = inf;
				isQP = true;
			}
			else {
				xl(i - 1) = -inf;
				xu(i - 1) = inf;
			}
		}

		isQP = false;//TOCHANGE

		if (isQP) {
			shared_ptr<QuadProgMosek> program_ = make_shared <QuadProgMosek>();
			program_->setParamInt(MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT);
			program_->setParamInt(MSK_IPAR_LOG, 10);
			program_->setParamInt(MSK_IPAR_LOG_FILE, 1);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_DFEAS, 1e-8);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_INFEAS, 1e-10);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_MU_RED, 1e-8);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_NEAR_REL, 1e3);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_PFEAS, 1e-8);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_REL_GAP, 1e-8);

			program_->setNumberOfVariables(num_joints);
			program_->setLowerVariableBound(xl);
			program_->setUpperVariableBound(xu);

			program_->setObjectiveMatrix(A.sparseView());
			program_->setObjectiveVector(b.segment(6, num_joints));

			bool success = program_->solve();
			VectorXd sol = program_->getPrimalSolution();
			VectorXd phi = JJ * sol;

			for (int i = 1; i < (int)boxes.size(); i++) {
				auto box = boxes[i];
				// Update joint angles
				box->setRotationAngle(h * sol(i - 1));
				box->setThetadot(sol(i - 1));
				// Don't forget to update twists as well, we will use it to compute forces
				box->setTwist(phi.segment<6>(6 * i));			
			}
		}
		else{

			// Use the result of KKT
			for (int i = 0; i < (int)boxes.size(); i++) {
				auto box = boxes[i];
				if (i != 0) {
					// Update joint angles
					box->setRotationAngle(h * thetadotlist(i - 1));
					//box->setThetadot(x(5 + i));
					box->setThetadot(newthetadotlist(i - 1));
					//cout << endl << box->getThetadot() << endl;
					//cout << endl << box->getAngle() << endl;
					thetalistnew(i - 1) = box->getAngle();

					// Don't forget to update twists as well, we will use it to compute forces
					box->setTwist(phi.segment<6>(6 * i));
				}
			}
			//MatrixXd Jnew = getGlobalJacobian(thetalistnew);
			//cout << "jNEW" << Jnew << endl;

			//VectorXd phi = Jnew * newthetadotlist;
			for (int i = 0; i < (int)boxes.size(); i++) {
				auto box = boxes[i];
				if (i != 0) {				
					// Don't forget to update twists as well, we will use it to compute forces
					//box->setTwist(phi.segment<6>(6 * (i - 1)));
				}
			}
		}
	}
	else {
		VectorXd boxtwists;
		boxtwists.resize((int)6 * boxes.size());
		boxtwists.setZero();

		// Maximal Coordinate
		int j = 0;
		for (int i = 0; i < (int)boxes.size(); i++) {
			auto box = boxes[i];

			A.block<6, 6>(6 * i, 6 * i) = box->getMassMatrix();
			M.block<6, 6>(6 * i, 6 * i) = box->getMassMatrix();
			boxtwists.segment<6>(6 * i) = box->getTwist();
			b.segment<6>(6 * i) = h * box->getForce();
			if (i == 0) {
				Matrix6d I;
				I.setIdentity();
				A.block<6, 6>(6 * (int)boxes.size(), 0) = I;
				A.block<6, 6>(0, 6 * (int)boxes.size()) = I;
			}
			else {
				auto joint = box->getJoint();
				Matrix6d Ad_J_P = Rigid::adjoint(joint->getE_P_J().inverse());
				Matrix6d Ad_J_C = -Rigid::adjoint(joint->getE_C_J().inverse());
				int id_P = box->getParent()->getIndex();
				int id_C = box->getIndex();

				Matrix5x6d G_J_P;
				G_J_P.block<2, 6>(0, 0) = Ad_J_P.block<2, 6>(0, 0);
				G_J_P.block<3, 6>(2, 0) = Ad_J_P.block<3, 6>(3, 0);

				Matrix5x6d G_J_C;
				G_J_C.block<2, 6>(0, 0) = Ad_J_C.block<2, 6>(0, 0);
				G_J_C.block<3, 6>(2, 0) = Ad_J_C.block<3, 6>(3, 0);

				A.block<5, 6>(6 * (int)boxes.size() + 6 + 5 * j, 6 * id_P) = G_J_P;
				A.block<6, 5>(6 * id_P, 6 * (int)boxes.size() + 6 + 5 * j) = G_J_P.transpose();

				A.block<5, 6>(6 * (int)boxes.size() + 6 + 5 * j, 6 * id_C) = G_J_C;
				A.block<6, 5>(6 * id_C, 6 * (int)boxes.size() + 6 + 5 * j) = G_J_C.transpose();

				j++;
			}
		}
		//cout << "phi_n" << boxtwists << endl;
		double kk = 0.5 * boxtwists.transpose() * M * boxtwists;
		//cout << "kk_n" << kk << endl;

		MatrixXd M_s = Spring::computeMassMatrix(springs, (int)boxes.size(), isReduced);

		// Add mass matrix of spring to A matrix
		A.block(0, 0, 6 * (int)boxes.size(), 6 * (int)boxes.size()) += M_s;

		VectorXd b_s = Spring::computeGravity(springs, (int)boxes.size(), isReduced);
		b.segment(0, 6 * (int)boxes.size()) += A.block(0, 0, 6 * (int)boxes.size(), 6 * (int)boxes.size()) * boxtwists + b_s * h;

		x = A.ldlt().solve(b);
		//cout << x << endl;
		//cout << A << endl;
		//cout << b << endl;
		MatrixXd AA(22, 22);

		AA.setZero();
		AA.block<12, 12>(0, 0) = M.block<12, 12>(0, 0);
		AA.block<5, 6>(12, 0) = A.block<5, 6>(24, 6);
		AA.block<5, 12>(17, 0) = A.block<5, 12>(29, 6);
		AA.block<6, 5>(0, 12) = A.block<5, 6>(24, 6).transpose();
		AA.block<12, 5>(0, 17) = A.block<5, 12>(29, 6).transpose();
		//cout << "AA"<< AA << endl;
		VectorXd xxxx = AA.ldlt().solve(b.segment<22>(6));
		//cout << "bbb" << b.segment<22>(6) << endl;
		//cout << xxxx << endl;

		//cout << "phi_n+1" << x.segment<18>(0) << endl;
		double kkk = 0.5 * x.segment<18>(0).transpose() * M * x.segment<18>(0);
		//cout << "kk_n+1" << kkk << endl;

		// Update boxes
		for (int i = 0; i < (int)boxes.size(); i++) {
			auto box = boxes[i];
			box->setTwist(x.segment<6>(6 * i));
		}
	}
}


void SymplecticIntegrator::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, shared_ptr<MatrixStack> P) const
{	
	prog->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	MV->pushMatrix();
	glColor3f(1.0, 0.0, 0.0); // black
	glPointSize(3.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < (int)debug_points.size(); ++i) {
		Vector3f p = debug_points[i]->x.cast<float>();
		glVertex3f(p(0), p(1), p(2));
	}	
	glEnd();
	MV->popMatrix();
	prog->unbind();
}

SymplecticIntegrator::~SymplecticIntegrator() {

}