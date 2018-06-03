#pragma once
#ifndef MUSCLEMASS_SRC_SCENE_H_
#define MUSCLEMASS_SRC_SCENE_H_

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <json.hpp>

class Particle;
class MatrixStack;
class Program;
class Shape;
class Rigid;
class Solver;
class Joint;
class Vector;
class WrapSphere;
class WrapCylinder;
class WrapDoubleCylinder;

typedef Eigen::Matrix<double, 12, 1> Vector12d;

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
	
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;
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
	double V0;
	double K0;

	//double y_t;
	Vector12d phi; // twist
	std::vector < double > Kvec;
	std::vector < double > Vvec;
	std::vector < double > Tvec;
	//std::vector < double > y_t_vec;
	std::vector < Vector12d > Twist_vec;

	std::vector<double> y;
	std::vector<double> yp;

	std::shared_ptr<Shape> sphereShape;
	std::shared_ptr<Shape> boxShape;
	std::shared_ptr<Shape> cylinderShape;


	std::vector< std::shared_ptr<Rigid> > boxes;
	std::shared_ptr<Solver> solver;
	std::vector< std::shared_ptr<WrapDoubleCylinder> > wrap_doublecylinders;

	nlohmann::json js;
};

#endif // MUSCLEMASS_SRC_SCENE_H_