#include <iostream>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Cloth.h"
#include "Particle.h"
#include "Spring.h"
#include "MatrixStack.h"
#include "Program.h"
#include "GLSL.h"
#include <stdio.h>
#include <mosek.h>

using namespace std;
using namespace Eigen;
typedef Eigen::Triplet<double> T;

// // prints log output from MOSEK to the terminal
// static void MSKAPI printstr(void *handle, MSKCONST char str[]){
// 	printf("%s", str);
// }


shared_ptr<Spring> createSpring(const shared_ptr<Particle> p0, const shared_ptr<Particle> p1, double E)
{
	auto s = make_shared<Spring>(p0, p1);
	s->E = E;
	Vector3d x0 = p0->x;
	Vector3d x1 = p1->x;
	Vector3d dx = x1 - x0;
	s->L = dx.norm();
	return s;
}

Cloth::Cloth(int rows, int cols,
			 const Vector3d &x00,
			 const Vector3d &x01,
			 const Vector3d &x10,
			 const Vector3d &x11,
			 double mass,
			 double stiffness,
			 const Vector2d &damping)
{
	assert(rows > 1);
	assert(cols > 1);
	assert(mass > 0.0);
	assert(stiffness > 0.0);
	
	this->rows = rows;
	this->cols = cols;
	this->damping = damping;
	
	// Create particles
	n = 0;
	double r = 0.02; // Used for collisions
	int nVerts = rows*cols;
	for(int i = 0; i < rows; ++i) {
		double u = i / (rows - 1.0);
		Vector3d x0 = (1 - u)*x00 + u*x10;
		Vector3d x1 = (1 - u)*x01 + u*x11;
		for(int j = 0; j < cols; ++j) {
			double v = j / (cols - 1.0);
			Vector3d x = (1 - v)*x0 + v*x1;
			auto p = make_shared<Particle>();
			particles.push_back(p);
			p->r = r;
			p->x = x;
			p->v << 0.0, 0.0, 0.0;
			p->m = mass/(nVerts);
			// Pin two particles
			if(i == 0 && (j == 0 || j == cols-1)) {
				p->fixed = true;
				p->i = -1;
			} else {
				p->fixed = false;
				p->i = n;
				n += 3;
			}
		}
	}
	
	// Create x springs
	for(int i = 0; i < rows; ++i) {
		for(int j = 0; j < cols-1; ++j) {
			int k0 = i*cols + j;
			int k1 = k0 + 1;
			springs.push_back(createSpring(particles[k0], particles[k1], stiffness));
		}
	}
	
	// Create y springs
	for(int j = 0; j < cols; ++j) {
		for(int i = 0; i < rows-1; ++i) {
			int k0 = i*cols + j;
			int k1 = k0 + cols;
			springs.push_back(createSpring(particles[k0], particles[k1], stiffness));
		}
	}
	
	// Create shear springs
	for(int i = 0; i < rows-1; ++i) {
		for(int j = 0; j < cols-1; ++j) {
			int k00 = i*cols + j;
			int k10 = k00 + 1;
			int k01 = k00 + cols;
			int k11 = k01 + 1;
			springs.push_back(createSpring(particles[k00], particles[k11], stiffness));
			springs.push_back(createSpring(particles[k10], particles[k01], stiffness));
		}
	}

	// Build system matrices and vectors
	M.resize(n,n);
	K.resize(n,n);
	v.resize(n);
	f.resize(n);

	
	// Build vertex buffers
	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();
	posBuf.resize(nVerts*3);
	norBuf.resize(nVerts*3);
	updatePosNor();
	// Texture coordinates (don't change)
	for(int i = 0; i < rows; ++i) {
		for(int j = 0; j < cols; ++j) {
			texBuf.push_back(i/(rows-1.0));
			texBuf.push_back(j/(cols-1.0));
		}
	}
	// Elements (don't change)
	for(int i = 0; i < rows-1; ++i) {
		for(int j = 0; j < cols; ++j) {
			int k0 = i*cols + j;
			int k1 = k0 + cols;
			// Triangle strip
			eleBuf.push_back(k0);
			eleBuf.push_back(k1);
		}
	}
}

Cloth::~Cloth()
{
}

void Cloth::tare()
{
	for(int k = 0; k < (int)particles.size(); ++k) {
		particles[k]->tare();
	}
}

void Cloth::reset()
{
	for(int k = 0; k < (int)particles.size(); ++k) {
		particles[k]->reset();
	}
	updatePosNor();
}

void Cloth::updatePosNor()
{
	// Position
	for(int i = 0; i < rows; ++i) {
		for(int j = 0; j < cols; ++j) {
			int k = i*cols + j;
			Vector3d x = particles[k]->x;
			posBuf[3*k+0] = x(0);
			posBuf[3*k+1] = x(1);
			posBuf[3*k+2] = x(2);
		}
	}
	// Normal
	for(int i = 0; i < rows; ++i) {
		for(int j = 0; j < cols; ++j) {
			// Each particle has four neighbors
			//
			//      v1
			//     /|\
			// u0 /_|_\ u1
			//    \ | /
			//     \|/
			//      v0
			//
			// Use these four triangles to compute the normal
			int k = i*cols + j;
			int ku0 = k - 1;
			int ku1 = k + 1;
			int kv0 = k - cols;
			int kv1 = k + cols;
			Vector3d x = particles[k]->x;
			Vector3d xu0, xu1, xv0, xv1, dx0, dx1, c;
			Vector3d nor(0.0, 0.0, 0.0);
			int count = 0;
			// Top-right triangle
			if(j != cols-1 && i != rows-1) {
				xu1 = particles[ku1]->x;
				xv1 = particles[kv1]->x;
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			// Top-left triangle
			if(j != 0 && i != rows-1) {
				xu1 = particles[kv1]->x;
				xv1 = particles[ku0]->x;
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			// Bottom-left triangle
			if(j != 0 && i != 0) {
				xu1 = particles[ku0]->x;
				xv1 = particles[kv0]->x;
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			// Bottom-right triangle
			if(j != cols-1 && i != 0) {
				xu1 = particles[kv0]->x;
				xv1 = particles[ku1]->x;
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			nor /= count;
			nor.normalize();
			norBuf[3*k+0] = nor(0);
			norBuf[3*k+1] = nor(1);
			norBuf[3*k+2] = nor(2);
		}
	}
}

void Cloth::step(double h, const Vector3d &grav, const vector< shared_ptr<Particle> > spheres)
{
	M.setZero();
	K.setZero();
	v.setZero();
	f.setZero();

	//
	// IMPLEMENT ME!
	//

	// Sparse Version
	// Construct matrix A and vector b while looping through particles and springs
	vector<T> A_;
	Eigen::VectorXd b;
	b.resize(n);
	b.setZero();

	for(int i=0; i<particles.size(); i++) {
		int idx = particles[i]->i;
		if(idx != -1){
			double mass = particles[i]->m;
			Vector3d vel = particles[i]->v;

			// filling vector v for initial guess
			v.segment<3>(idx) = vel;

			for(int j = 0; j<3; j++){
				A_.push_back(T(idx+j, idx+j, mass));
				b(idx + j) += mass * vel(j) + h * grav(j) * mass;

			}
		}
	}


	for(int i=0; i<springs.size(); i++){

		int idx0 = springs[i]->p0->i;
 		int idx1 = springs[i]->p1->i;

 		Vector3d dx = springs[i]->p1->x - springs[i]->p0->x;

 		// The vector Fs for one spring

 		double l = sqrt(pow(dx(0),2)+pow(dx(1),2)+pow(dx(2),2));
 		Vector3d fs = (springs[i]->E)*(l - springs[i]->L)/l * dx;


 		if(idx0 != -1){
 			for (int t=0; t<3; t++){
 				b(idx0+t) += h * fs(t);
 			}
 		}

 		if(idx1 != -1){
 			for (int t=0; t<3; t++){
 				b(idx1+t) -= h * fs(t);
 			}
 		}
 		

 		// The matrix Ks for one spring

 		Matrix3d xxT = dx * dx.transpose();
 		double xTx = dx.transpose() * dx;
 		Matrix3d I = Matrix3d::Identity();
 		Matrix3d ks = springs[i]->E/pow(l,2) * ( (1-(l-springs[i]->L)/l)*xxT + (l - springs[i]->L)/l * xTx * I);

 		if( idx0 != -1 && idx1 != -1 ) {

 			for(int k = 0; k < 3; k++ ){
 				for(int g = 0; g < 3; g++ ){
 					
 					// val = beta * h^2 * ks(k,g)

 					double val = damping(1) * pow(h, 2) * ks(k,g);

 					// filling the diagonal values of K
 					A_.push_back(T(idx0 + k, idx0 + g, val));
 					A_.push_back(T(idx1 + k, idx1 + g, val));

 					// filling the negated, off-diagonal values of K

 					A_.push_back(T(idx0 + k, idx1 + g, -val));
 					A_.push_back(T(idx1 + k, idx0 + g, -val));
 				}
 			}
 		}
	}

	Eigen::SparseMatrix<double> A(n, n);
	A.setFromTriplets(A_.begin(), A_.end());

	ConjugateGradient< SparseMatrix<double> > cg;
	cg.setMaxIterations(25);
	cg.setTolerance(1e-3);
	cg.compute(A);
	VectorXd x = cg.solveWithGuess(b, v);

	// cout << b << endl;
	// cout << endl;
	// for (int p = 0;  p< particles.size(); p++) {
	// 	int index = particles[p]->i;
	// 	if(index != -1){
	// 		double mass = particles[p]->m;
	// 		Matrix3d A;
	// 		A << mass, 0, 0,
	// 			0, mass, 0,
	// 			0, 0, mass;
	// 		M.block<3,3>(index, index) = A;

	// 		Vector3d B;
	// 		B << particles[p]->v;
	// 		v.segment<3>(index) = B;

	// 		f.segment<3>(index) = mass * grav;
	// 	}

	// }

 // 	for (int i = 0; i<springs.size(); i++){
 // 		int idx0 = springs[i]->p0->i;
 // 		int idx1 = springs[i]->p1->i;
 // 		Vector3d dx = springs[i]->p1->x - springs[i]->p0->x;

 // 		// fs
 // 		double l = sqrt(pow(dx(0),2)+pow(dx(1),2)+pow(dx(2),2));
 // 		Vector3d fs = (springs[i]->E)*(l - springs[i]->L)/l * dx;

 // 		// ks
 // 		Matrix3d xxT = dx * dx.transpose();
 // 		double xTx = dx.transpose() * dx;
 // 		Matrix3d I = Matrix3d::Identity();
 // 		Matrix3d ks = springs[i]->E/pow(l,2) * ( (1-(l-springs[i]->L)/l)*xxT + (l - springs[i]->L)/l * xTx * I);

 // 		int ei = 0;
 // 		int si = 0;
 // 		if(idx0 > idx1){
 // 			ei = idx0;
 // 			si = idx1;
 // 		}else{
 // 			ei = idx1;
 // 			si = idx0;
 // 		}

 // 		if(idx0 !=-1 && idx1 != -1){
 // 			K.block<3,3>(si, si) += ks;
 // 			K.block<3,3>(ei, ei) += ks;
 // 			K.block<3,3>(si, ei) -= ks;
 // 			K.block<3,3>(ei, si) -= ks;
 // 		}
 		


 // 		if(idx0 != -1){

 // 			f.segment<3>(idx0) += fs;
 // 		}

 // 		if(idx1 != -1){
 // 			f.segment<3>(idx1) -= fs;
 // 		}
	
	// }

	 // cout << M << endl;
	 // cout << v << endl;
	// cout << f << endl;
	// cout << (M + damping(1)*pow(h,2)*K)<< endl;
	// cout << M*v + h*f <<endl;

	// VectorXd x = (M + damping(1)*pow(h,2)*K).ldlt().solve(M*v + h*f);

	

	// Collision detection
	// Vector3d s0;
	// double xs, xb, ys, yb, zs, zb =0.0;

	// for(int i = 0; i< spheres.size(); i++){
	// 	Vector3d pos_sphere = spheres[i]->x;
	// 	if(pos_sphere(0)>xb){
	// 		xb = pos_sphere(0);
	// 	}
	// 	if(pos_sphere(0)<xs){
	// 		xs = pos_sphere(0);
	// 	}
	// 	if(pos_sphere(1)>yb){
	// 		yb = pos_sphere(1);
	// 	}
	// 	if(pos_sphere(1)<ys){
	// 		ys = pos_sphere(1);
	// 	}
	// 	if(pos_sphere(2)>zb){
	// 		zb = pos_sphere(2);
	// 	}
	// 	if(pos_sphere(2)<zs){
	// 		zs = pos_sphere(2);
	// 	}
	// }
	// s0 << (xb - xs)/2, (yb - ys)/2, (zb - zs)/2;


	// for(int i = 0; i< particles.size(); i++){
	// 	Vector3d pos_cloth = particles[i]->x;
	// 	Vector3d v_cloth = particles[i]->v;
	// 	Vector3d nor = pos_cloth - s0;
	// 	double norv = sqrt(pow(nor(0),2)+pow(nor(1),2)+ pow(nor(2),2));


	// 	if(norv < 0.1)
	// 	{	Vector3d proj = v_cloth.dot(nor)/pow(norv,2)*nor;

	// 		Vector3d pos_new = s0 + nor/norv;
	// 		particles[i]->x = pos_new;
			
	// 		particles[i]->v += proj;

			

	// 	}
		
	// }

 	// cout << x << endl;
 	// cout << endl;

	for(int i = 0; i < particles.size(); i++){
		
		if(particles[i]->i != -1){
			particles[i]->v = x.segment<3>(particles[i]->i);
		}
		// particles[i]->v = x.segment<3>(particles[i]->i);
	}

	for(int i = 0; i < particles.size(); i++){

		if(particles[i]->i != -1){
			particles[i]->x += particles[i]->v * h;

		}
	}

	// Collision Detection and Response (simple version)
	int n_col = 0;
	std::vector<int> index_cols;//used for the constraint matrix

	vector<T> C_; // for constraints
	//A_.push_back(T(idx+j, idx+j, mass));
	int row = 0;

	for (int i = 0; i < spheres.size(); i++){
	 	
	 	shared_ptr<Particle> p_s = spheres[i];
      	for (int j = 0; j < particles.size(); j++) {

         	shared_ptr<Particle> p = particles[j];
         	Vector3d dx = p->x - p_s->x;
         	Vector3d dx_n = dx.normalized();

         	// Collisions with the sphere

         	if (dx.norm() <= p->r + p_s->r) {
         		p->isCol = true;

            	// p->x = ( (p_s->r + p->r) * dx_n + p_s->x);
            
            	Vector3d v_normal = (p->v.dot(dx_n) * dx).normalized();

            	// save the normal of each particle that has a collision event

            	p->normal = v_normal;// normalize
            	// filling the constraint matrix
            	for(int t=0; t<3; t++){
            		C_.push_back(T(row, (p->i) + t, v_normal(t)));
            	}
            	row ++;

            	index_cols.push_back(p->i);

            	// p->v -= v_normal;
         	}
      	}
	}

	int num_collisions = index_cols.size();
	Eigen::SparseMatrix<double> C(num_collisions, n);
	C.setFromTriplets(C_.begin(), C_.end());
	VectorXd lc = VectorXd:: Zero(num_collisions); // lc m-long vector of 0 for linear inequality
	VectorXd uc = VectorXd:: Zero(num_collisions); // m-long +Inifity
	VectorXd lx = VectorXd:: Zero(n); 
	VectorXd ux = VectorXd:: Zero(n);
	VectorXd cc = -b;

	for(int i = 0; i< n; i++){
		lx(i) = - Infinity;
		ux(i) = + Infinity;
		uc(i) = + Infinity;
	}


	VectorXd results;
	igl::mosek::MosekData mosek_data;

	bool r = mosek_quadprog(A, cc, 0, C, lc, uc, lx, ux, mosek_data, results);
	// Collision Detection and Response (complex version) using mosek


	// int numvar = n;
	// int numcon = 1;
	// int i,j =0;
	// // Init
	// MSKenv_t env = NULL;
	// MSKtask_t task = NULL;
	// MSKrescodee r;

	// // Create the mosek environment
	// r = MSK_makeenv(&env,NULL);
	// if( r==MSK_RES_OK){
	// 	//Create the optimization task
	// 	r = MSK_maketask(env,numcon,numvar,&task);

	// 	if(r==MSK_RES_OK){
	// 		r=MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);

	// 	}

	// 	if(r== MSK_RES_OK){
	// 		r=MSK_appendcons(task,numcon);
	// 	}

	// 	if(r==MSK_RES_OK){
	// 		r=MSK_appendvars(task,numvar);
	// 	}

	// 	for(j=0; j<numvar && r==MSK_RES_OK; j++){
	// 		//Set the linear term b_j in the objective
	// 		if(r == MSK_RES_OK){
	// 			r = MSK_putcj(task,j,-b[j]);
	// 		}

	// 		//Set the bounds on var j
	// 		if(r == MSK_RES_OK){
	// 			r = MSK_putvarbound(task,
	// 								j,
	// 								MSK_BK_FR, // the constraint is free
	// 								-MSK_INFINITY,
	// 								+MSK_INFINITY);

	// 		}
	// 	}

	// 	vector<int> row_idx;
	// 	for(int i=0; i<index_cols.size(); i++){
	// 		row_idx.push_back(i);
	// 	}

	// 	//Input column j of matrix A
	// 	if(r == MSK_RES_OK){
	// 		for(int i = 0; i < index_cols.size(); i++){

	// 			int t = index_cols[i];
	// 			Vector3d normal = particles[(int)t/3]->normal;

	// 			// filling normal vector into A
	// 			for(int k=0; k< 3; k++){

	// 				r = MSK_putacol(task,
	// 								t+k,
	// 								1,
	// 								&row_idx[0],
	// 								&normal(k));
	// 			}
	// 			// r = MSK_putacol(task,
	// 			// 				t,//column index
	// 			// 				1,// number of non-zero in col j
	// 			// 				&row_idx,// pointer to row indexes of col j
	// 			// 				normal(0));

	// 		}													
	// 	}

	// 	//Set bounds of constraint
	// 	r = MSK_putconbound(task, 0, MSK_BK_LO, 0, +MSK_INFINITY);

	// 	//Q matrix
	// 	std::vector<MSKint32t> qsubi;
	// 	std::vector<MSKint32t> qsubj;
	// 	vector<double> qval;
	// 	double xx[numvar];

	// 	qsubi.resize(A_.size());
	// 	qsubj.resize(A_.size());
	// 	qval.resize(A_.size());

	// 	int count = 0;

	// 	for(int i = 0; i < A_.size(); i++){
	// 		if(A_[i].row() >= A_[i].col()){
	// 			qsubi[count] = (MSKint32t) A_[i].row();
	// 			qsubj[count] = (MSKint32t) A_[i].col();
	// 			qval[count] = A_[i].value();
	// 			++ count;
	// 		}
	// 	}

	// 	// 
	// 	r = MSK_putqobj(task, count, &qsubi[0], &qsubj[0], &qval[0]);
	// 	if (r == MSK_RES_OK )
	// 	{
	// 		MSKrescodee trmcode;
	// 		//run optimizer
	// 		r = MSK_optimizetrm(task, &trmcode);
	// 		MSKsolstae solsta;
	// 		if( r== MSK_RES_OK){
	// 			r = MSK_getsolsta(task,
	// 							  MSK_SOL_ITR,
	// 							  &solsta);
	// 			switch(solsta){
	// 				case MSK_SOL_STA_OPTIMAL:
	// 				case MSK_SOL_STA_NEAR_OPTIMAL:
	// 					MSK_getxx(task, MSK_SOL_ITR, xx);
	// 					break;
	// 				default:
	// 					printf("Other solution status.");
	// 					break;
	// 			}
	// 		}
	// 	}

	// 	MSK_deletetask(&task);
	// 	MSK_deleteenv(&env);
	// }

	// Update position and normal buffers
	updatePosNor();
}

void Cloth::init()
{
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	
	glGenBuffers(1, &norBufID);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size()*sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	
	glGenBuffers(1, &texBufID);
	glBindBuffer(GL_ARRAY_BUFFER, texBufID);
	glBufferData(GL_ARRAY_BUFFER, texBuf.size()*sizeof(float), &texBuf[0], GL_STATIC_DRAW);
	
	glGenBuffers(1, &eleBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size()*sizeof(unsigned int), &eleBuf[0], GL_STATIC_DRAW);
	
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	
	assert(glGetError() == GL_NO_ERROR);
}

void Cloth::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p) const
{
	// Draw mesh
	glUniform3fv(p->getUniform("kdFront"), 1, Vector3f(1.0, 0.0, 0.0).data());
	glUniform3fv(p->getUniform("kdBack"),  1, Vector3f(1.0, 1.0, 0.0).data());
	MV->pushMatrix();
	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	int h_pos = p->getAttribute("aPos");
	glEnableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	int h_nor = p->getAttribute("aNor");
	glEnableVertexAttribArray(h_nor);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size()*sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	for(int i = 0; i < rows; ++i) {
		glDrawElements(GL_TRIANGLE_STRIP, 2*cols, GL_UNSIGNED_INT, (const void *)(2*cols*i*sizeof(unsigned int)));
	}
	glDisableVertexAttribArray(h_nor);
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();
}
