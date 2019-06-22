#ifndef UNSCENTEDKALMANFILTER_MEASUREMENT_PACKAGE_H
#define UNSCENTEDKALMANFILTER_MEASUREMENT_PACKAGE_H

#include <Eigen/Eigen>

class MeasurementPackage
{
public:
    long timestamp_;

    enum SensorType
    {
        LASER,
        RADAR
    }sensor_type_;

    Eigen::VectorXd raw_measurements_;
};

#endif //UNSCENTEDKALMANFILTER_MEASUREMENT_PACKAGE_H
