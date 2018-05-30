#include "WrapDoubleCylinder.h"

void WrapDoubleCylinder::compute()
{
	// compute Matrix U and V
	Eigen::Vector3d OP = this->point_P - this->point_U;
	OP = OP / OP.norm();
	Eigen::Vector3d vec_Z_U = vec_z_U / vec_z_U.norm();
	Eigen::Vector3d vec_X_U = vec_Z_U.cross(OP);
	vec_X_U = vec_X_U / vec_X_U.norm();
	Eigen::Vector3d vec_Y_U = vec_Z_U.cross(vec_X_U);
	vec_Y_U = vec_Y_U / vec_Y_U.norm();

	Eigen::Vector3d OS = this->point_S - this->point_V;
	OS = OS / OS.norm();
	Eigen::Vector3d vec_Z_V = vec_z_V / vec_z_V.norm();
	Eigen::Vector3d vec_X_V = vec_Z_V.cross(OS);
	vec_X_V = vec_X_V / vec_X_V.norm();
	Eigen::Vector3d vec_Y_V = vec_Z_V.cross(vec_X_V);
	vec_Y_V = vec_Y_V / vec_Y_V.norm();

	this->M_U.resize(3, 3);
	this->M_U.row(0) = vec_X_U.transpose();
	this->M_U.row(1) = vec_Y_U.transpose();
	this->M_U.row(2) = vec_Z_U.transpose();

	this->M_V.resize(3, 3);
	this->M_V.row(0) = vec_X_V.transpose();
	this->M_V.row(1) = vec_Y_V.transpose();
	this->M_V.row(2) = vec_Z_V.transpose();

	// step 1: compute H and T
	Eigen::Vector3d pv = this->M_V * (this->point_P - this->point_V);
	Eigen::Vector3d sv = this->M_V * (this->point_S - this->point_V);

	double denom_h = pv(0)*pv(0) + pv(1)*pv(1);
	double denom_t = sv(0)*sv(0) + sv(1)*sv(1);
	double Rv = this->radius_V;

	double root_h = sqrt(denom_h - Rv*Rv);
	double root_t = sqrt(denom_t - Rv*Rv);

	Eigen::Vector3d h(0.0, 0.0, 0.0);
	Eigen::Vector3d t(0.0, 0.0, 0.0);
	h(0) = (pv(0) * Rv*Rv + Rv * pv(1) * root_h) / denom_h;
	h(1) = (pv(1) * Rv*Rv - Rv * pv(0) * root_h) / denom_h;
	t(0) = (sv(0) * Rv*Rv - Rv * sv(1) * root_t) / denom_t;
	t(1) = (sv(1) * Rv*Rv + Rv * sv(0) * root_t) / denom_t;

	this->status = wrap;

	std::complex<double> ht_i = 1.0 - 0.5 *
		((h(0) - t(0)) * (h(0) - t(0))
			+ (h(1) - t(1)) * (h(1) - t(1))) / (Rv*Rv);
	std::complex<double> ph_i = 1.0 - 0.5 *
		((pv(0) - h(0)) * (pv(0) - h(0))
			+ (pv(1) - h(1)) * (pv(1) - h(1))) / (Rv*Rv);
	std::complex<double> ts_i = 1.0 - 0.5 *
		((t(0) - sv(0)) * (t(0) - sv(0))
			+ (t(1) - sv(1)) * (t(1) - sv(1))) / (Rv*Rv);

	double ht_xy = abs(Rv * acos(ht_i));
	double ph_xy = abs(Rv * acos(ph_i));
	double ts_xy = abs(Rv * acos(ts_i));

	h(2) = pv(2) + (sv(2) - pv(2)) * ph_xy / (ph_xy + ht_xy + ts_xy);
	t(2) = sv(2) - (sv(2) - pv(2)) * ts_xy / (ph_xy + ht_xy + ts_xy);

	Eigen::Vector3d H = this->M_V.transpose() * h + this->point_V;
	Eigen::Vector3d T = this->M_V.transpose() * t + this->point_V;
	Eigen::Vector3d H0 = H;

	Eigen::Vector3d q(0.0, 0.0, 0.0);
	Eigen::Vector3d g(0.0, 0.0, 0.0);
	Eigen::Vector3d Q, G;

	double len = 0.0;

	for (int i = 0; i < 30; i++)
	{
		len = 0.0;

		// step 2: compute Q and G
		Eigen::Vector3d pu = this->M_U * (this->point_P - this->point_U);
		Eigen::Vector3d hu = this->M_U * (H - this->point_U);

		double denom_q = pu(0)*pu(0) + pu(1)*pu(1);
		double denom_g = hu(0)*hu(0) + hu(1)*hu(1);
		double Ru = -this->radius_U;

		double root_q = sqrt(denom_q - Ru*Ru);
		double root_g = sqrt(denom_g - Ru*Ru);

		q(0) = (pu(0) * Ru*Ru + Ru * pu(1) * root_q) / denom_q;
		q(1) = (pu(1) * Ru*Ru - Ru * pu(0) * root_q) / denom_q;
		g(0) = (hu(0) * Ru*Ru - Ru * hu(1) * root_g) / denom_g;
		g(1) = (hu(1) * Ru*Ru + Ru * hu(0) * root_g) / denom_g;

		std::complex<double> qg_i = 1.0 - 0.5 *
			((q(0) - g(0)) * (q(0) - g(0))
				+ (q(1) - g(1)) * (q(1) - g(1))) / (Ru*Ru);
		std::complex<double> pq_i = 1.0 - 0.5 *
			((pu(0) - q(0)) * (pu(0) - q(0))
				+ (pu(1) - q(1)) * (pu(1) - q(1))) / (Ru*Ru);
		std::complex<double> gh_i = 1.0 - 0.5 *
			((g(0) - hu(0)) * (g(0) - hu(0))
				+ (g(1) - hu(1)) * (g(1) - hu(1))) / (Ru*Ru);

		double qg_xy = abs(Rv * acos(qg_i));
		double pq_xy = abs(Rv * acos(pq_i));
		double gh_xy = abs(Rv * acos(gh_i));
		len += qg_xy;

		q(2) = pu(2) + (hu(2) - pu(2)) * pq_xy / (pq_xy + qg_xy + gh_xy);
		g(2) = hu(2) - (hu(2) - pu(2)) * gh_xy / (pq_xy + qg_xy + gh_xy);

		Q = this->M_U.transpose() * q + this->point_U;
		G = this->M_U.transpose() * g + this->point_U;

		// step 3: compute H based on G and T
		Eigen::Vector3d gv = this->M_V * (G - this->point_V);

		double denom_h = gv(0)*gv(0) + gv(1)*gv(1);
		double root_h = sqrt(denom_h - Rv*Rv);

		h = Eigen::Vector3d(0.0, 0.0, 0.0);
		h(0) = (gv(0) * Rv*Rv + Rv * gv(1) * root_h) / denom_h;
		h(1) = (gv(1) * Rv*Rv - Rv * gv(0) * root_h) / denom_h;

		std::complex<double> ht_i = 1.0 - 0.5 *
			((h(0) - t(0)) * (h(0) - t(0))
				+ (h(1) - t(1)) * (h(1) - t(1))) / (Rv*Rv);
		gh_i = 1.0 - 0.5 *
			((gv(0) - h(0)) * (gv(0) - h(0))
				+ (gv(1) - h(1)) * (gv(1) - h(1))) / (Rv*Rv);
		std::complex<double> ts_i = 1.0 - 0.5 *
			((t(0) - sv(0)) * (t(0) - sv(0))
				+ (t(1) - sv(1)) * (t(1) - sv(1))) / (Rv*Rv);

		double ht_xy = abs(Rv * acos(ht_i));
		gh_xy = abs(Rv * acos(gh_i));
		double ts_xy = abs(Rv * acos(ts_i));
		len += ht_xy;

		h(2) = gv(2) + (sv(2) - gv(2)) * gh_xy / (gh_xy + ht_xy + ts_xy);
		t(2) = sv(2) - (sv(2) - gv(2)) * ts_xy / (gh_xy + ht_xy + ts_xy);

		H = this->M_V.transpose() * h + this->point_V;
		T = this->M_V.transpose() * t + this->point_V;

		len += (G - H).norm();

		double dist = (H - H0).norm();
		if (dist == 0) break;

		H0 = H;
	}

	this->path_length = len;
	this->point_q = q;
	this->point_g = g;
	this->point_h = h;
	this->point_t = t;
	/*
	std::cout << Q.transpose() << std::endl << G.transpose() << std::endl
	<< H.transpose() << std::endl << T.transpose() << std::endl;
	*/
}


Eigen::MatrixXd WrapDoubleCylinder::getPoints(int num_points)
{
	double theta_q = atan(this->point_q(1) / this->point_q(0));
	if (this->point_q(0) < 0.0)
		theta_q += PI;

	double theta_g = atan(this->point_g(1) / this->point_g(0));
	if (this->point_g(0) < 0.0)
		theta_g += PI;

	double theta_h = atan(this->point_h(1) / this->point_h(0));
	if (this->point_h(0) < 0.0)
		theta_h += PI;

	double theta_t = atan(this->point_t(1) / this->point_t(0));
	if (this->point_t(0) < 0.0)
		theta_t += PI;

	Eigen::MatrixXd points(3, 3 * num_points + 1);

	double theta_s, theta_e, z_s, z_e;

	// q to g
	if (theta_q < theta_g)
	{
		theta_s = theta_q; theta_e = theta_g;
		z_s = this->point_q(2); z_e = this->point_g(2);
	}
	else
	{
		theta_s = theta_g; theta_e = theta_q;
		z_s = this->point_g(2); z_e = this->point_q(2);
	}

	if (theta_e - theta_s > theta_s + 2 * PI - theta_e)
	{
		double tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2 * PI;
		tmp = z_s; z_s = z_e; z_e = tmp;
	}

	int col = 0;
	double z_i = z_s, dz = (z_e - z_s) / num_points;
	for (double i = theta_s; i <= theta_e + 0.001;
		i += (theta_e - theta_s) / num_points)
	{	
		if (col == num_points + 1) {
			break;
		}
		Eigen::Vector3d point = this->M_U.transpose() *
			Eigen::Vector3d(this->radius_U * cos(i),
				this->radius_U * sin(i), z_i) +
			this->point_U;
		z_i += dz;
		points.col(col++) = point;
	}

	// g to h
	Eigen::Vector3d G = this->M_U.transpose() * this->point_g + this->point_U;
	Eigen::Vector3d H = this->M_V.transpose() * this->point_h + this->point_V;
	Eigen::Vector3d diff = H - G;

	for (int i = 1; i < num_points; i++)
	{
		Eigen::Vector3d point = G + diff / num_points * i;
		points.col(col++) = point;
	}

	// h to t
	if (theta_h < theta_t)
	{
		theta_s = theta_h; theta_e = theta_t;
		z_s = this->point_h(2); z_e = this->point_t(2);
	}
	else
	{
		theta_s = theta_t; theta_e = theta_h;
		z_s = this->point_t(2); z_e = this->point_h(2);
	}

	if (theta_e - theta_s > theta_s + 2 * PI - theta_e)
	{
		double tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2 * PI;
		tmp = z_s; z_s = z_e; z_e = tmp;
	}

	z_i = z_s;
	dz = (z_e - z_s) / num_points;
	for (double i = theta_s; i <= theta_e + 0.001;
		i += (theta_e - theta_s) / num_points)
	{
		if (col == 3 * num_points + 1) {
			break;
		}
		Eigen::Vector3d point = this->M_V.transpose() *
			Eigen::Vector3d(this->radius_V * cos(i),
				this->radius_V * sin(i), z_i) +
			this->point_V;
		z_i += dz;
		points.col(col++) = point;
	}

	return points;
}