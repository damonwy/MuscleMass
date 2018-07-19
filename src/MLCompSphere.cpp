#include "MLCompSphere.h"

#include <fstream>

#include "MatrixStack.h"
#include "Program.h"

MLCompSphere::MLCompSphere(std::shared_ptr<MLBody> parent) : MLComp(parent){

}

MLCompSphere::~MLCompSphere() {

}

void MLCompSphere::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog) const {

}

void MLCompSphere::save(std::ofstream &ofs) {

}