#include "RKF45Integrator.h"

#include "Rigid.h"
#include "Joint.h"
#include "MatlabDebug.h"
#include "Spring.h"
#include "Particle.h"

#include <iostream>
#include <iomanip>

# include <cstdlib>

# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;
using namespace Eigen;
#include <unsupported/Eigen/MatrixFunctions> // TODO: avoid using this later, write a func instead

RKF45Integrator::RKF45Integrator(vector< shared_ptr<Rigid> > _boxes, vector< shared_ptr<Spring> > _springs, bool _isReduced) :
	num_joints(_boxes.size() - 1),
	boxes(_boxes),
	springs(_springs),
	isReduced(_isReduced),
	epsilon(1e-8),
	grav(0.0, -9.8, 0.0)
{

	if (isReduced) {
		m = 6 * (int)boxes.size();
		n = 6 + 1 * ((int)boxes.size() - 1);
		M.resize(m, m);
		J.resize(m, n);
		f.resize(m);
		M.setZero();
		J.setZero();
		f.setZero();

	}
	else {
		n = 6 * (int)boxes.size() + 6 + 5 * ((int)boxes.size() - 1);
	}

	A.resize(n, n);
	x.resize(n);
	b.resize(n);

	A.setZero();
	x.setZero();
	b.setZero();
}