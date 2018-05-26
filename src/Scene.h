#pragma once
#ifndef MUSCLEMASS_SRC_SCENE_H_
#define MUSCLEMASS_SRC_SCENE_H_

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Particle;
class MatrixStack;
class Program;
class Shape;
class Rigid;
class Solver;
class Joint;

class Scene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Scene();
	virtual ~Scene();
	
	void load(const std::string &RESOURCE_DIR);
	void init();
	void tare();
	void reset();
	void step();
	
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog) const;
	void computeEnergy();
	void saveData(int num_steps);
	double getTime() const { return t; }
	
private:
	double t;
	double h;
	int step_i;
	Eigen::Vector3d grav;
	int n_step;
	Eigen::VectorXd theta_list;
	double t_start;
	double t_stop;

	double V; // potential energy
	double K; // kinetic energy
	std::vector < double > Kvec;
	std::vector < double > Vvec;
	std::vector < double > Tvec;

	std::vector<double> y;
	std::vector<double> yp;

	std::shared_ptr<Shape> boxShape;
	std::vector< std::shared_ptr<Rigid> > boxes;
	std::shared_ptr<Solver> solver;
};

#endif // MUSCLEMASS_SRC_SCENE_H_