#include "MLObject.h"

#include <iostream>
#include <fstream>

#include "MLError.h"
#include "MLWorld.h"
#include "MatrixStack.h"
#include "Program.h"

using json = nlohmann::json;

int MLObject::UIDMAX = 1000;

MLObject::MLObject(): m_name(""), m_uid(0)
{
}

MLObject::~MLObject() {

}

MLError MLObject::load(const json &elem, std::shared_ptr<MLWorld> world, const std::string &/*folder*/) {
	try{
		m_name = elem.at("name").get<std::string>();
		m_uid = elem.at("uid").get<int>();
	}
	catch (json::exception& e) {
		// output exception information
		std::cout << "message: " << e.what() << '\n' << "exception id: " << e.id << std::endl;
		return MLError("error parsing object");
	}

	return MLError();
}

void MLObject::save(std::ofstream &ofs) const {
	ofs << " NAME " << getName() << std::endl;
	ofs << " UID " << getUID() << std::endl;
}

int MLObject::genNextUID() {
	++MLObject::UIDMAX;
	return MLObject::UIDMAX;
}
