#include <iostream>
#include <Eigen/Eigen>
#include "unscented_kalman_filter.h"

int main() {

    Eigen::MatrixXd P(3,3);
    P << 6, 0, 0, 0, 4, 0, 0, 0, 7;
    Eigen::MatrixXd L( P.llt().matrixL() );
    std::cout << L << std::endl;

    return 0;
}