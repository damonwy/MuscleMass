#pragma once
#ifndef MUSCLEMASS_SRC_SCENE_H_
#define MUSCLEMASS_SRC_SCENE_H_

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <json.hpp>
#include "MLCommon.h"


class MatrixStack;
class Program;

class Rigid;
class Solver;
class Spring;
class Joint;
class Vector;
class MLWorld;

class MLScene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	explicit MLScene();
	virtual ~MLScene();

	void load(const std::string &RESOURCE_DIR);
	void init();
	void tare();
	void reset();
	void step();

	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;

	double getTime() const { return t; }

private:
	double t;
	double h;
	Eigen::Vector3d grav;

	nlohmann::json js;

	std::shared_ptr<MLWorld> m_world;
};

#endif // MUSCLEMASS_SRC_SCENE_H_