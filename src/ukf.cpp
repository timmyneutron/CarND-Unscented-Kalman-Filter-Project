#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // size of state vector
  n_x_ = 5;

  // size of augmented state vector
  n_aug_ = n_x_ + 2;

  // number of sigma points
  n_sig_ = 2 * n_aug_ + 1;

  // lambda coefficient for sigma point generation
  lambda_ = 3 - n_aug_;

  // sigma point weights
  weights_ = VectorXd(n_sig_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < n_sig_; ++i)
  {
  	weights_(i) = 0.5 / (lambda_ + n_aug_);
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

	// if first measurement, initialize state vector
	if (!is_initialized_)
	{
		if (meas_package.sensor_type_ == MeasurementPackage::LASER)
		{
			double px = meas_package.raw_measurements_(0);
			double py = meas_package.raw_measurements_(1);

			x_ << px, py, 0, 0, 0;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
		{
			double rho = meas_package.raw_measurements_(0);
			double phi = meas_package.raw_measurements_(1);

			double px = rho * cos(phi);
			double py = rho * sin(phi);

			x_ << px, py, 0, 0, 0;
		}

		// initialize time stamp
		time_us_ = meas_package.timestamp_;

		is_initialized_ = true;
		return;
	}

	// define delta t as time since last measurement
	double delta_t = (meas_package.timestamp_ - time_us_) / 1e6;

	// reset time_us_ as time of last measurement
	time_us_ = meas_package.timestamp_;

	// run prediction step
	Prediction(delta_t);

	// run measurement update
	if (meas_package.sensor_type_ == MeasurementPackage::LASER)
	{
		UpdateLidar(meas_package);
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
		UpdateRadar(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

	// define augmented P matrix and square root (A) matrix
	MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(5, 5) = pow(std_a_, 2);
	P_aug(6, 6) = pow(std_yawdd_, 2);
	MatrixXd A = P_aug.llt().matrixL();

	// generate sigma points
	MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
	VectorXd x_aug = VectorXd(n_aug_);

	x_aug << x_, 0, 0;
	Xsig_aug.col(0) = x_aug;

	for (int i = 0; i < n_aug_; ++i)
	{
		Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
		Xsig_aug.col(i+n_aug_+1) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
	}

	// initialize matrix for predicted sigma points
	Xsig_pred_ = MatrixXd(n_x_, n_sig_);

	// initialize state transition vector
	VectorXd F = VectorXd(n_x_);

	// initialize noise vector
	VectorXd nu = VectorXd(n_x_);

	// fill F and nu vectors
	for (int i = 0; i < n_sig_; ++i)
	{
		double px = Xsig_aug(0, i);
		double py = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double psi = Xsig_aug(3, i);
		double psi_dot = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_psi_ddot = Xsig_aug(6, i);

		// if psi_dot = 0, use constant velocity model,
		// otherwise use CTRV model
		if (fabs(psi_dot) < 0.001)
		{
			F(0) = v * cos(psi) * delta_t;
			F(1) = v * sin(psi) * delta_t;
		}
		else
		{
			F(0) = v / psi_dot * (sin(psi + psi_dot * delta_t) - sin(psi));
			F(1) = v / psi_dot * (-cos(psi + psi_dot * delta_t) + cos(psi));
		}

		F(2) = 0;
		F(3) = psi_dot * delta_t;
		F(4) = 0;

		double delta_t2 = pow(delta_t, 2);

		nu(0) = 0.5 * delta_t2 * cos(psi) * nu_a;
		nu(1) = 0.5 * delta_t2 * sin(psi) * nu_a;
		nu(2) = delta_t * nu_a;
		nu(3) = 0.5 * delta_t2 * nu_psi_ddot;
		nu(4) = delta_t * nu_psi_ddot;

		// calculate predicted sigma point based on F and nu vectors
		Xsig_pred_.col(i) = Xsig_aug.col(i).head(n_x_) + F + nu;
	}

	// initialize predicted state vector
	x_ << 0, 0, 0, 0, 0;

	// calculate predicted state vector as weighted average of
	// predicted sigma points
	for (int i = 0; i < n_sig_; ++i)
	{
		x_ += weights_(i) * Xsig_pred_.col(i);
	}

	// calculate predicted covariance matrix
	P_ = MatrixXd::Zero(n_x_, n_x_);

	for (int i = 0; i < n_sig_; ++i)
	{
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		P_ += weights_(i) * x_diff * x_diff.transpose();
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	MatrixXd Zsig = MatrixXd(2, n_sig_);

	// fill prediction sigma points matrix
	for (int i = 0; i < n_sig_; ++i)
	{
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);

		Zsig.col(i) << px, py;
	}

	// calculate prediction mean as weighted average of sigma points
	VectorXd z_pred = VectorXd::Zero(2);

	for (int i = 0; i < n_sig_; ++i)
	{
		z_pred += weights_(i) * Zsig.col(i);
	}

	// calculate Kalman gain
	MatrixXd S = MatrixXd::Zero(2, 2);
	MatrixXd R = MatrixXd(2, 2);
	R << pow(std_laspx_, 2), 0,
		 0, pow(std_laspy_, 2);

	for (int i = 0; i < n_sig_; ++i)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S += weights_(i) * z_diff * z_diff.transpose();
	}

	S += R;

	MatrixXd T = MatrixXd::Zero(5, 2);

	for (int i = 0; i < n_sig_; ++i)
	{
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		VectorXd z_diff = Zsig.col(i) - z_pred;
		T += weights_(i) * x_diff * z_diff.transpose();
	}

	MatrixXd K = T * S.inverse();

	// update state vector and state covariance matrix
	x_ += K * (meas_package.raw_measurements_ - z_pred);
	P_ -= K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
