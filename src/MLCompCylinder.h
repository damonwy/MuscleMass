#pragma once
#ifndef MUSCLEMASS_SRC_MLCOMPSHPERE_H_
#define MUSCLEMASS_SRC_MLCOMPSPHERE_H_

#include "MLComp.h"


class MLCompCylinder : public MLComp
{
public:
	explicit MLCompCylinder(std::shared_ptr<MLBody> parent);
	virtual ~MLCompCylinder();

	virtual void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P) const;

	virtual MLError load(const nlohmann::json &elem, std::shared_ptr<MLWorld> world, const std::string &folder);
	virtual void save(std::ofstream &ofs);

	double getRadius() const { return m_radius; }
	double getHeight() const { return m_height; }

protected:
	double m_radius;
	double m_height;
};

#endif // MUSCLEMASS_SRC_MLCOMPSPHERE_H_