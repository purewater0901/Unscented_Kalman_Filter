#include "unscented_kalman_filter.h"
#include <iostream>

UKF::UKF(int n_x, int n_aug, bool use_laser, bool use_radar, bool use_nis)
{
    use_laser_ = use_laser;
    use_radar_ = use_radar;
    use_nis_ = use_nis;

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 3.0;
    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 1.6;
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

    is_initialized_ = false;

    n_x_ = n_x; //状態ベクトルの次元
    n_aug_ = n_aug; //ノイズも含めた次元
    n_aug_size_ = 2 * n_aug_ + 1;

    lambda_ = 3 - n_aug_;

    //initial state vector
    x_ = Eigen::VectorXd(n_x_);
    x_.fill(0.1);

    P_ = Eigen::MatrixXd(n_x_, n_aug_size_);

    predicted_sigma_pts_ = Eigen::MatrixXd(n_x_, n_aug_size_);

    //init count forNormalized Innovation Squared
    update_count_ = 0;
    nis_thresh_count_ = 0;
}

void UKF::ProcessMeasurement(MeasurementPackage& meas_package)
{
    //初期化
    if(!is_initialized_)
    {
        weights_ = Eigen::VectorXd(n_aug_size_);
        weights_(0) = (lambda_)/(lambda_ + n_aug_);
        for(int i=1; i<n_aug_size_; ++i)
            weights_(i) = 0.5 / (n_aug_ + lambda_);

        if(meas_package.sensor_type_ == MeasurementPackage::LASER)
            x_.head(2) = meas_package.raw_measurements_;
        else if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
        {
            double range = meas_package.raw_measurements_(0);
            double theta = meas_package.raw_measurements_(1);
            double vel = meas_package.raw_measurements_(2);

            x_(0) = range * cos(theta);
            x_(1) = range * sin(theta);
            x_(2) = vel;
        }

        is_initialized_ = true;
        return;
    }

    /* Prediction */
    double dt = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;
    previous_timestamp_ = meas_package.timestamp_;

    //Predict the objects position at the current time step
    UKF::Prediction(dt);

    //Update
    if(meas_package.sensor_type_ == MeasurementPackage::LASER)
        UKF::UpdateLidar(meas_package);
    else if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
        UKF::UpdateRadar(meas_package);
}


void UKF::SigmaPointPrediction(Eigen::MatrixXd &sigma_pts, Eigen::MatrixXd &pred_sigma_pts, double delta_t)
{
    pred_sigma_pts.fill(0.0);

    //predict sigma points
    for(int i=0; i<n_aug_size_; ++i)
    {
        double pos_x = sigma_pts(0, i);
        double pos_y = sigma_pts(1, i);
        double vel = sigma_pts(2, i);
        double yaw = sigma_pts(3,i);
        double yaw_dot = sigma_pts(4,i);
        double a_pos = sigma_pts(5,i); // acceleration process noise
        double a_yaw_dot = sigma_pts(6,i); // change in angle process noise

        // check for divide by zero
        if (fabs(yaw_dot) > 0.001) {
            pos_x += vel/yaw_dot * (sin(yaw + yaw_dot * delta_t) - sin(yaw));
            pos_y += vel/yaw_dot * (cos(yaw) - cos(yaw + yaw_dot * delta_t));
        } else {
            pos_x += vel * delta_t * cos(yaw);
            pos_y += vel * delta_t * sin(yaw);
        }

        //write predicted sigma point into right column
        pred_sigma_pts(0,i) = pos_x + 0.5 * a_pos * delta_t * delta_t * cos(yaw);
        pred_sigma_pts(1,i) = pos_y + 0.5 * a_pos * delta_t * delta_t * sin(yaw);
        pred_sigma_pts(2,i) = vel + a_pos * delta_t;
        pred_sigma_pts(3,i) = yaw + yaw_dot * delta_t + 0.5 * a_yaw_dot * delta_t * delta_t;
        pred_sigma_pts(4,i) = yaw_dot + a_yaw_dot * delta_t;
    }
}

void UKF::PredictMeanAndCovariance(Eigen::VectorXd &x, Eigen::MatrixXd &P, Eigen::MatrixXd &pred_sigma_pts,
                                   const int yaw_pos=-1)
{
    x.fill(0.0);
    for(int i=0; i<n_aug_size_; ++i)
        x += weights_[i] * pred_sigma_pts.col(i);

    //共分散行列
    P.fill(0.0);
    for(int i=0; i<n_aug_size_; i++) //sigmaポイント全てに対して差を取る
    {
        //状態の差
        Eigen::VectorXd x_diff = pred_sigma_pts.col(i) - x;
        if(yaw_pos >= 0)
        {
            //angle normalization
            while (x_diff(yaw_pos)> M_PI) x_diff(yaw_pos) -= 2.*M_PI;
            while (x_diff(yaw_pos)<-M_PI) x_diff(yaw_pos) += 2.*M_PI;
        }

        P+=weights_(i) * x_diff * x_diff.transpose();
    }
}

void UKF::Prediction(double delta_t)
{
    Eigen::VectorXd x_aug = Eigen::VectorXd(n_aug_, n_aug_);
    x_aug << x_, 0, 0;

    //Calculate the augmentation matrix
    Eigen::MatrixXd P_aug = Eigen::MatrixXd(n_aug_, n_aug_);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(5,5) = std_a_ * std_a_;
    P_aug(6,6) = std_yawdd_ * std_yawdd_;

    Eigen::MatrixXd sigma_pts = Eigen::MatrixXd(n_aug_, n_aug_size_);
    sigma_pts.col(0) = x_aug;

    //create square root matrix
    Eigen::MatrixXd sqrt_P = P_aug.llt().matrixL();

    //Generate sigma points
    for(int i=0; i<n_aug_; ++i)
    {
        sigma_pts.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * sqrt_P.col(i);
        sigma_pts.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * sqrt_P.col(i);
    }

    /*
     * Predict the sigma point values for this the time step
     */
    UKF::SigmaPointPrediction(sigma_pts, predicted_sigma_pts_, delta_t);

    //calculate mean position and covariance
    UKF::PredictMeanAndCovariance(x_, P_, predicted_sigma_pts_, 3);

}

//Sは観測共分散行列
void UKF::UpdateState(Eigen::VectorXd &x, Eigen::MatrixXd &P, Eigen::MatrixXd &pred_sigma_pts, Eigen::MatrixXd &S,
                      Eigen::MatrixXd &Zsig, Eigen::VectorXd &z_pred, const int &n_z,
                      MeasurementPackage &meas_package)
{
    //観測行列の共分散
    Eigen::MatrixXd Tc = Eigen::MatrixXd(n_x_ , n_z);

    Tc.fill(0.0);
    for(int i=0; i<n_aug_size_; ++i)
    {
        //residual
        //Zsigはセンサーによって予測されたシグマポイントの点
        //z_predは現在いると思われている点
        Eigen::VectorXd z_diff = Zsig.col(i) - z_pred;
        Eigen::VectorXd x_diff = pred_sigma_pts.col(i) - x;
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K
    Eigen::MatrixXd K = Tc * S.inverse();
    //residual(観測値と予測値の差)
    Eigen::VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

    if(use_nis_)
        NISState(S, z_diff, n_z);

    x += K * z_diff;
    P = P - K*S*K.transpose();
}








