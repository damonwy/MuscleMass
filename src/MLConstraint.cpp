#include "MLConstraint.h"

#include <sstream>

#include "MLContact.h"

MLConstraint::MLConstraint() : m_indexBilateral(0), m_indexContact(0) {
	m_mu = 0.5;
}

MLConstraint::~MLConstraint() {
	while (!m_contacts.empty()){
		delete m_contacts.back();
		m_contacts.pop_back();
	}
}

void MLConstraint::init() {
	for (int i = 0; i < (int)m_contacts.size(); ++i) {
		//m_contacts[i]->init();
	}
}

void MLConstraint::draw() {
	// can draw contacts if needed
	// now ignore
}

unsigned int MLConstraint::getNumRowsContact() const {
	return m_contacts.size();
}

unsigned int MLConstraint::getNumRowsFriction() const {
	// Number of tangents (both inactive and active)
	unsigned int k = 0;
	for (int i = 0; i < (int)m_contacts.size(); ++i) {
		//k += m_contacts[i]->getNumRowsFriction();
	}
	return k;
}

void MLConstraint::setIndexBilateral(unsigned int *indexBilateral) {
	m_indexBilateral = *indexBilateral;
	*indexBilateral += getNumRowsBilateral();
}

void MLConstraint::setIndexContact(unsigned int *indexContact) {
	m_indexContact = *indexContact;
	// indices for submatrices for each contact
	for (int i = 0; i < (int)m_contacts.size(); ++i) {
		//m_contacts[i]->setIndexContact(indexContact, i);
	}
}

void MLConstraint::setIndexFriction(unsigned int *indexFriction) {
	m_indexFriction = *indexFriction;
	for (int i = 0; i < (int)m_contacts.size(); ++i) {
		//m_contacts[i]->setIndexFriction(indexFriction);
	}
}

int MLConstraint::updateCache(const Eigen::VectorXd *q) {
	int status = 1;
	for (int i = 0; i < (int)m_contacts.size(); ++i) {
		//status *= m_contacts[i]->updateCache(q);
	}
	return status;
}


std::string MLConstraint::getColLabel() const
{
	std::stringstream label;
	label << 'V' << getUID();
	return label.str();
}

std::string MLConstraint::getRowLabel(char prefix) const
{
	std::stringstream label;
	label << prefix << getUID();
	return label.str();
}