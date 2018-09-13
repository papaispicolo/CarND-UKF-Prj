#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  ///* State dimension
  n_x_ = x_.size();

  ///* Augmented state dimension
  n_aug_ = n_x_ + 2; // with 2 noise parameters

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  ///* time when the state is true, in us
  // time_us_ = 0

  ///* Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_ / (n_aug_ + lambda_);
  for (int i = 1; i <= n_aug_; i++) {
    weights_(i) = 0.5 / (n_aug_ + lambda_);
    weights_(i + n_aug_) = weights_(i);
  }

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  ///* initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
      cout << "UKF: " << endl;
      x_ << 1, 1, 0, 0, 0;
      P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

      if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
          /**
           Convert radar from polar to cartesian coordinates and initialize state.
           */
          double rho = meas_package.raw_measurements_[0];
          double phi = meas_package.raw_measurements_[1];

          x_(0) = rho * cos(phi);
          x_(1) = rho * sin(phi);
      }
      else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
          /**
           Initialize state.
           */
          x_(0) = meas_package.raw_measurements_[0];
          x_(1) = meas_package.raw_measurements_[1];
      }
      time_us_ = meas_package.timestamp_;

      // done initializing, no need to predict or update
      is_initialized_ = true;
      return;
  }

  // Ignore radar measurement
  if (!use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) return;
  // Ignore laser measurement
  if (!use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) return;

  // prediction
  double dt;
  dt = (meas_package.timestamp_ - time_us_) / 1000000.0;    //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;
  Prediction(dt);

  // update
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      UpdateRadar(meas_package);
  } else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //////////////////////////////////////////////
  // Create Augmented sigma points
  //////////////////////////////////////////////

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_ + 1) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  //////////////////////////////////////////////
  // Predict Augmented sigma points
  //////////////////////////////////////////////
  for (int i = 0; i <= 2 * n_aug_; i++) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // Predicted state values
    double px_p, py_p;

    // Avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = px + v / yawd * (sin (yaw + yawd * delta_t) - sin(yaw));
        py_p = py + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
        px_p = px + v * delta_t * cos(yaw);
        py_p = py + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // Add noise
    px_p += 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p += 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p += nu_a * delta_t;

    yaw_p += 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p += nu_yawdd * delta_t;

    // Write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  // Predicted state
  x_ = Xsig_pred_ * weights_;

  // Compute the predicted state comvariance matrix
  P_.fill(0.0);
  for (int i = 0; i <= 2 * n_aug_; i++) {
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // Angle normalization
    while (x_diff(3) >= M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  // Create measurement matrix H
  MatrixXd H = MatrixXd(2, n_x_);
  H.fill(0.0);
  H(0, 0) = 1;
  H(1, 1) = 1;

  MatrixXd S = H * P_ * H.transpose();
  S(0, 0) += std_laspx_ * std_laspx_;
  S(1, 1) += std_laspy_ * std_laspy_;

  MatrixXd Sinv = S.inverse();
  MatrixXd K = P_ * H.transpose() * Sinv;

  // Update state mean and covariance matrix
  VectorXd z_diff = meas_package.raw_measurements_ - x_.head(2);
  x_ += K * z_diff;
  P_ -= K * H * P_;

  // Compute laser NIS
  double nis = z_diff.transpose() * Sinv * z_diff;
  cout << "Lidar (laser) NIS: " << nis << endl;
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

  // Transforme Sigma points to the measurement space
  //  and compute predicted measurement mean z_pred
  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);

  for (int i = 0; i <= 2 * n_aug_; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double vx = v * cos(yaw);
    double vy = v * sin(yaw);

    Zsig(0, i) = sqrt(px * px + py * py);
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = (px * vx + py * vy) / Zsig(0, i);
  }

  MatrixXd z_pred = Zsig * weights_;

  // Create innovation covariacne matrix S
  // and cross-corrleation matrix T
  MatrixXd S = MatrixXd(3, 3);
  MatrixXd T = MatrixXd(n_x_, 3);
  S.fill(0.0);
  T.fill(0.0);

  for (int i = 0; i <= 2 * n_aug_; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1) > M_PI)  z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    // state diff
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI)  x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
    T += weights_(i) * x_diff * z_diff.transpose();
  }

  // Add noise
  S(0, 0) += std_radrd_  * std_radrd_;
  S(1, 1) += std_radphi_ * std_radphi_;
  S(2, 2) += std_radrd_  * std_radrd_;

  MatrixXd Sinv = S.inverse();

  // Kalman gain K;
  MatrixXd K = T * Sinv;

  // Residual
  VectorXd z_res = meas_package.raw_measurements_ - z_pred;

  // Angle normalization
  while (z_res(1) > M_PI) z_res(1) -= 2. * M_PI;
  while (z_res(1) < -M_PI) z_res(1) += 2. * M_PI;

  // Update state mean and covariance matrix
  x_ += K * z_res;
  P_ -= K * S * K.transpose();

  // Compute radar NIS
  double nis = z_res.transpose() * Sinv * z_res;
  cout << "radar NIS: " << nis << endl;
}
