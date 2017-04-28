#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	// initialize
	VectorXd rmse = VectorXd::Zero(4);
	VectorXd residual = VectorXd(4);
	int n_estimations = estimations.size();

	for (int i=0; i<n_estimations; i++) {
		// calculate residual
		residual = estimations[i] - ground_truth[i];

		// square residual
		residual = residual.array() * residual.array();

		// sum residuals
		rmse += residual;
	}

	// find mean of squared residuals
	rmse /= n_estimations;

	// take square root of mean
	rmse = rmse.array().sqrt();

	// return rmse
	return rmse;
}
