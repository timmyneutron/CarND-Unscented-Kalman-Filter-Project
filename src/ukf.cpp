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
  std_a_ = 3.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.5;

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

  // augmented state matrix
  x_aug_ = VectorXd(n_aug_);

  // augmented covariance matrix
  P_aug_ = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug_(5, 5) = std_a_ * std_a_;
  P_aug_(6, 6) = std_yawdd_ * std_yawdd_;

  // augmented sigma points matrix
  Xsig_aug_ = MatrixXd(n_aug_, n_sig_);

  // sigma point weights
  weights_ = VectorXd::Constant(n_sig_, 0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // measurement noise covariance matrix for laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
  			  0, std_laspy_ * std_laspy_;

  // measurement noise covariance matrix for radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
  			  0, std_radphi_ * std_radphi_, 0,
  			  0, 0, std_radrd_ * std_radrd_;

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

	// define delta_t as time since last measurement (in seconds)
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
	P_aug_.topLeftCorner(n_x_, n_x_) = P_;
	MatrixXd A = P_aug_.llt().matrixL();

	// generate sigma points
	x_aug_ << x_, 0, 0;
	Xsig_aug_.col(0) = x_aug_;

	for (int i = 0; i < n_aug_; ++i)
	{
		Xsig_aug_.col(i+1) = x_aug_ + sqrt(lambda_ + n_aug_) * A.col(i);
		Xsig_aug_.col(i+n_aug_+1) = x_aug_ - sqrt(lambda_ + n_aug_) * A.col(i);
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
		double px = Xsig_aug_(0, i);
		double py = Xsig_aug_(1, i);
		double v = Xsig_aug_(2, i);
		double psi = Xsig_aug_(3, i);
		double psi_dot = Xsig_aug_(4, i);
		double nu_a = Xsig_aug_(5, i);
		double nu_psi_ddot = Xsig_aug_(6, i);

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
		Xsig_pred_.col(i) = Xsig_aug_.col(i).head(n_x_) + F + nu;
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
	
	// fill prediction sigma points matrix
	MatrixXd Zsig = MatrixXd(2, n_sig_);

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

	for (int i = 0; i < n_sig_; ++i)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S += weights_(i) * z_diff * z_diff.transpose();
	}

	S += R_laser_;

	MatrixXd T = MatrixXd::Zero(5, 2);

	for (int i = 0; i < n_sig_; ++i)
	{
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		VectorXd z_diff = Zsig.col(i) - z_pred;
		T += weights_(i) * x_diff * z_diff.transpose();
	}

	MatrixXd S_inverse = S.inverse();
	MatrixXd K = T * S_inverse;

	// update state vector and state covariance matrix
	VectorXd z_meas_diff = meas_package.raw_measurements_ - z_pred;
	x_ += K * z_meas_diff;
	P_ -= K * S * K.transpose();

	// calculate normalized innovation squared for lidar
	NIS_laser_ = z_meas_diff.transpose() * S_inverse * z_meas_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	
	// fill prediction points matrix
	MatrixXd Zsig = MatrixXd(3, n_sig_);

	for (int i = 0; i < n_sig_; ++i)
	{
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double psi = Xsig_pred_(3, i);

		double vx = v * cos(psi);
		double vy = v * sin(psi);

		double rho = sqrt(px * px + py * py);
		double phi = atan2(py, px);
		double rho_dot = 0.5 / sqrt(px * px + py * py) * (2 * px * vx + 2 * py * vy);
	
		Zsig.col(i) << rho, phi, rho_dot;
	}
	
	// calculate prediction mean as weighted average of sigma points
	VectorXd z_pred = VectorXd::Zero(3);

	for (int i = 0; i < n_sig_; ++i)
	{
		z_pred += weights_(i) * Zsig.col(i);
	}

	// calculate Kalman gain
	MatrixXd S = MatrixXd::Zero(3, 3);

	for (int i = 0; i < n_sig_; ++i)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S += weights_(i) * z_diff * z_diff.transpose();
	}

	S += R_radar_;

	MatrixXd T = MatrixXd::Zero(5, 3);

	for (int i = 0; i < n_sig_; ++i)
	{
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		VectorXd z_diff = Zsig.col(i) - z_pred;

		T += weights_(i) * x_diff * z_diff.transpose();
	}

	MatrixXd S_inverse = S.inverse();
	MatrixXd K = T * S_inverse;

	// update state vector and state covariance matrix
	VectorXd z_meas_diff = meas_package.raw_measurements_ - z_pred;
	x_ += K * z_meas_diff;
	P_ -= K * S * K.transpose();

	// calculate normalized innovation squared for radar
	NIS_radar_ = z_meas_diff.transpose() * S_inverse * z_meas_diff;
}
