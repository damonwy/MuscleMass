#pragma once
#ifndef MUSCLEMASS_SRC_MLCOMPSHPERE_H_
#define MUSCLEMASS_SRC_MLCOMPSPHERE_H_

#include "MLComp.h"

class MLCompSphere : public MLComp
{
public:
	explicit MLCompSphere(std::shared_ptr<MLBody> parent);
	virtual ~MLCompSphere();

	virtual void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog) const;

	virtual MLError load(const nlohmann::json &elem, std::shared_ptr<MLWorld> world, const std::string &folder);
	virtual void save(std::ofstream &ofs);

	double getRadius() const { return m_radius; }

protected:
	double m_radius;

};

#endif // MUSCLEMASS_SRC_MLCOMPSPHERE_H_