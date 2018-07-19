#pragma once
#ifndef MUSCLEMASS_SRC_MLCOMP_H_
#define MUSCLEMASS_SRC_MLCOMP_H_

#include "MLCommon.h"
#include "MLObject.h"

class MLBody;
class MLWorld;
class Shape;


class MLComp : public MLObject, public std::enable_shared_from_this<MLComp> {
public:
	explicit MLComp(std::shared_ptr<MLBody> parent);
	virtual ~MLComp();

	virtual MLError load(const nlohmann::json &elem, std::shared_ptr<MLWorld> world, const std::string &folder);
	virtual void save(std::ofstream &ofs);

	virtual void init();
	virtual void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P);

	std::shared_ptr<MLBody> getParent() const { return m_parent; }
	const Eigen::Matrix4d * getTransformLocal() const { return &m_cToP; }
	MLError getTransformWorld(Eigen::Matrix4d &sToWord, const Eigen::Matrix4d *pToWorld) const;
	// std::shared_ptr<MLPolygon> getPolygon() const { return m_polygon; }
	const std::vector<int> * getCollisions() const { return &m_collisions; }

protected:
	std::shared_ptr<MLBody> m_parent;
	std::shared_ptr<Shape> m_shape;
	Eigen::Matrix4d m_cToP; // transform wrt parent
	// MLPolygon m_polygon;
	std::vector<int> m_collisions;
};

#endif // MUSCLEMASS_SRC_MLCOMP_H_
