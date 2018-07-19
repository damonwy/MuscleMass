#pragma once
#ifndef MUSCLEMASS_SRC_MLCOMMON_H_
#define MUSCLEMASS_SRC_MLCOMMON_H_

#include <json.hpp>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 12, 1> Vector12d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 3, 6> Matrix3x6d;
typedef Eigen::Matrix<double, 6, 3> Matrix6x3d;
typedef Eigen::Matrix<double, 5, 6> Matrix5x6d;
typedef Eigen::Matrix<double, 3, 12> Matrix3x12d;
typedef Eigen::Matrix<double, 12, 12> Matrix12d;

enum Integrator { RKF45, SYMPLECTIC };

// Eigen types to/from GLM types
glm::mat3 eigen_to_glm(const Eigen::Matrix3d &m);
glm::mat4 eigen_to_glm(const Eigen::Matrix4d &m);
glm::vec3 eigen_to_glm(const Eigen::Vector3d &v);
Eigen::Vector3d glm_to_eigen(const glm::vec3 &v);
Eigen::Matrix3d glm_to_eigen(const glm::mat3 &m);
Eigen::Matrix4d glm_to_eigen(const glm::mat4 &m);

#endif // MUSCLEMASS_SRC_MLCOMMON_H_