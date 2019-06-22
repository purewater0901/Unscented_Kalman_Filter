//
// Created by yutaka on 19/06/22.
//

#ifndef UNSCENTEDKALMANFILTER_UNSCENTED_KALMAN_FILTER_H
#define UNSCENTEDKALMANFILTER_UNSCENTED_KALMAN_FILTER_H

#include "measurement_package.h"
#include <Eigen/Eigen>

class UKF
{
public:
    //状態ベクトル(x, y, v, yaw, 角速度)
    Eigen::VectorXd x_;
    Eigen::MatrixXd P_;

    //状態の共分散行列
    UKF(int n_x, int n_aug, bool use_laser=true, bool use_radar=true, bool use_nis=false);

    virtual ~UKF() = default;

    void ProcessMeasurement(MeasurementPackage& meas_package);

private:

    bool is_initialized_;

    bool use_laser_;

    bool use_radar_;

    //nisを計算するかどうか
    bool use_nis_;

    Eigen::MatrixXd predicted_sigma_pts_;
    Eigen::VectorXd weights_;

    //加速度による縦方向の分散
    double std_a_;

    //yawの加速度の分散
    double std_yawdd_;

    //x方向の測定プロセスの分散
    double std_laspx_;
    //x方向の測定プロセスの分散
    double std_laspy_;

    //radarによる半径の測定ノイズのぶんさん
    double std_radr_;
    double std_radphi_;
    double std_radrd_;

    //状態の次元
    int n_x_;
    //拡張した状態の次元
    int n_aug_;

    //time when the state is true, in us
    long long previous_timestamp_;

    //sigma points size of the augmented matrix(全部で何個のシグマポイントを打つか)
    int n_aug_size_;

    //sigmaポイントを広げるパラメータ
    double lambda_;

    //updateを実行した回数
    long update_count_;
    //thresholdを超えているノイズの数
    long nis_thresh_count_;

    void Prediction(double delta_t);

    void UpdateLidar(MeasurementPackage& meas_package);

    void UpdateRadar(MeasurementPackage& meas_package);

    /*
     * Predicts the previous sigma points at the current time step
     */
    void SigmaPointPrediction(Eigen::MatrixXd& sigma_pts, Eigen::MatrixXd& pred_sigma_pts, double delta_t);

    //pred_sigma_ptsはシグマポイントを格納している
    void PredictMeanAndCovariance(Eigen::VectorXd& x, Eigen::MatrixXd& P, Eigen::MatrixXd& pred_sigma_pts, const int yaw_pos);

    //x 予測した平均ベクトル
    // P 予測した共分散行列
    //S 測定した共分散行列
    //Zsig センサーがよそしたシグマポイント
    //z_pred 予測した位置
    //n_z 何個のセンサーを使ったか
    //meas_package
    void UpdateState(Eigen::VectorXd& x, Eigen::MatrixXd& P, Eigen::MatrixXd& pred_sigma_pts,
                     Eigen::MatrixXd& S, Eigen::MatrixXd& Zsig, Eigen::VectorXd& z_pred, const int& n_z,
                     MeasurementPackage& meas_package);

    // S 測定した共分散行列
    //z_diff どれくらい予測値と測定値が離れているか
    void NISState(Eigen::MatrixXd& S, Eigen::VectorXd& z_diff, const int& n_z);
};


#endif //UNSCENTEDKALMANFILTER_UNSCENTED_KALMAN_FILTER_H
