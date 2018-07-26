#include "MLCommon.h"

glm::mat3 eigen_to_glm(const Eigen::Matrix3d &m) {
	return glm::make_mat3x3((const double *)m.data());
}

glm::mat4 eigen_to_glm(const Eigen::Matrix4d &m) {
	return glm::make_mat4x4((const double *)m.data());
}

glm::vec3 eigen_to_glm(const Eigen::Vector3d &v) {
	return glm::vec3(v[0], v[1], v[2]);
}

Eigen::Vector3d glm_to_eigen(const glm::vec3 &v) {
	return Eigen::Vector3d(v[0], v[1], v[2]);
}

Eigen::Matrix3d glm_to_eigen(const glm::mat3 &m) {
	const float *_m = glm::value_ptr(m);
	float m_cp[3 * 3];
	memcpy(m_cp, _m, sizeof(m_cp));
	Eigen::Map<Eigen::Matrix3f> result(m_cp);
	return result.cast<double>();
}

Eigen::Matrix4d glm_to_eigen(const glm::mat4 &m) {
	const float *_m = glm::value_ptr(m);
	float m_cp[4 * 4];
	memcpy(m_cp, _m, sizeof(m_cp));
	Eigen::Map<Eigen::Matrix4f> result(m_cp);
	return result.cast<double>();
}