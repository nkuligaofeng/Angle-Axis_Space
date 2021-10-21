//==============================================================================
/*
This header file is to provide the function to calculate the perpendicular curve
theory in rotation group SO(3).

\author    <gaofeng.li@iit.it>
\author    Gaofeng Li
\date      Sep 23, 2019
*/
//==============================================================================

#ifndef myRotPerpH
#define myRotPerpH

#include <math.h>
#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Dense>

#define PI 3.1415926

namespace myRotPerp{

    // @author
    //			Gaofeng Li: gaofeng.li@iit.it
    // @data
    //			23-Sep-2019
    // @description
    //			Return the rotation Matrix by given a vector in angle-axis space
    // @param
    //			Eigen::Vector3d& omega:  the vector in angle-axis space.
    // @return
    //			Eigen::Matrix3d: the corresponding rotation
    inline Eigen::Matrix3d ExpRotation(const Eigen::Vector3d& omega)
    {
        Eigen::Matrix3d Rot = Eigen::Matrix3d::Identity();
        double w_norm = omega.norm();
        if(w_norm != 0)
        {
            if(w_norm - PI > 0.005)
            {
                Logger::warning() << "Warning:!!!!!----Invalid input in ExpRotation(const Eigen::Vector3d& omega)--- the length of input vector is greater than PI!!!!!!!!!!!!" << Logger::endl();
            }
            Eigen::Matrix3d skewMat = (Eigen::Matrix3d() << 0.0, -omega.z()/w_norm, omega.y()/w_norm,
                    omega.z()/w_norm, 0.0, -omega.x()/w_norm,
                    -omega.y()/w_norm, omega.x()/w_norm, 0.0).finished();
            Rot = Eigen::Matrix3d::Identity() + sin(w_norm) * skewMat + (1 - cos(w_norm)) * skewMat * skewMat;
            return Rot;
        }
    }

    // @author
    //			Gaofeng Li: gaofeng.li@iit.it
    // @data
    //			09-May-2019
    // @description
    //			Return the corresponding vector in angle-axis space by given a rotation matrix
    // @param
    //			Eigen::Matrix3d& omega:  the given rotation matrix.
    // @return
    //			Eigen::Vector3d: the corresponding vector
    inline Eigen::Vector3d LogRotation(const Eigen::Matrix3d& Rot)
    {
        Eigen::Matrix3d temp1Rot;
        temp1Rot = Eigen::Matrix3d::Identity() - Rot.transpose() * Rot;//tempRot = I - R^T * R, to check whether the input is a rotation matrix
        if(temp1Rot.norm() > 0.005)// this threshold is to consider the calculation accuracy, so we can let it larger to gain the robustness
        {
            Logger::error() << "Warning!!!!!   ||R^T*R-I||=" << temp1Rot.norm() << ", Rot is not a rotation matrix, please check again!!!!" << Logger::endl();
            return Eigen::Vector3d(PI, PI, PI);
        }// check whether R is a rotation matrix
        else{
            // special case 1: the rotation matrix is Identity
            Eigen::Matrix3d temp2Rot;
            temp2Rot = Eigen::Matrix3d::Identity() - Rot;
            if(temp2Rot.norm() < 0.005)// special case 1: the rotation matrix is Identity. this threshold is to consider the calculation accuracy, so we can let it larger to gain the robustness
            {
                return Eigen::Vector3d(0.0, 0.0, 0.0);
            }// check whether R is identity
            else{
                Eigen::Matrix3d skewMat;
                skewMat = 0.5 * (Rot - Rot.transpose());
                Eigen::Vector3d omega_candi(-skewMat(1,2), skewMat(0,2), -skewMat(0,1));// candidate omega
                if(omega_candi.norm() < 0.00001)// special case 2, the angle is PI; choose the one in the positive half sphere; this threshold is to avoid "divide 0", so it can smaller
                {
                    double wx = sqrt( (Rot(0,0) + 1) * PI * PI/2 ); //  we always choose the point on the positive half sphere
                    double wy, wz;
                    if (wx < 0.00001)// this means that the vector is on the the circle that is the intersection of the y-z plane and the sphere
                    {
                        wx = 0;
                        wy = sqrt( (Rot(1, 1) + 1) * PI * PI / 2 );
                        if (wy < 0.00001)
                        {
                            wy = 0;
                            wz = sqrt((Rot(2, 2) + 1) * PI * PI / 2);
                        }
                        else
                        {
                            wz = sqrt((Rot(1, 2) + 1) * PI * PI / (2 * wy));
                        }
                    }
                    else// wx is not zero
                    {
                        wy = Rot(0, 1) * PI * PI / (2 * wx);
                        wz = Rot(0, 2) * PI * PI / (2 * wx);
                    }
                    return Eigen::Vector3d(wx, wy, wz);
                }//check whether Rot is the special case 2
                else{
                    double cosAngle = (Rot(0,0) + Rot(1,1) + Rot(2,2) - 1.0)/2.0;
                    // To avoid invalid input for acos
                    if(cosAngle < -1)
                    {
                        Logger::warning() << "Warning!!!!!---Invalid Input in LogRotation(const Eigen::Matrix3d& Rot):  the input for acos is less than -1---" << Logger::endl();
                        cosAngle = -1;
                    }
                    if(cosAngle > 1)
                    {
                        Logger::warning() << "Warning!!!!!---Invalid Input in LogRotation(const Eigen::Matrix3d& Rot):  the input for acos is larger than 1---" << Logger::endl();
                        cosAngle = 1;
                    }
                    double omegaNorm = acos( cosAngle );
                    Eigen::Vector3d omega;
                    omega = ( omegaNorm / omega_candi.norm() ) * omega_candi;
                    return omega;
                }// end else for checking whether Rot is the special case 2
            }//end else for checking whether R is identity
        } // end else for checking whether R is a rotation matrix
    }//LogRotation is finished

    // @author
    //			Gaofeng Li: gaofeng.li@iit.it
    // @data
    //			24-Sep-2019
    // @description
    //			Get the Foot of Perpendicular (FoP).
    // @param
    //			Eigen::Matrix3d& R_tar:  target rotation matrix
    //          Eigen::Matrix3d& R_cur:  current rotation matrix
    //          Eigen::Vector3d& kappa:  the given direction
    // @return
    //			Eigen::Matrix3d: the FoP
    inline Eigen::Matrix3d GetRPerp(const Eigen::Matrix3d& R_tar, const Eigen::Matrix3d& R_cur, const Eigen::Vector3d& kappa)
    {
        /**
         * Input check
         * ****/
        Eigen::Matrix3d tempRot;
        tempRot = Eigen::Matrix3d::Identity() - R_tar.transpose() * R_tar;//tempRot = I - R^T * R, to check whether the input is a rotation matrix
        if(tempRot.norm() > 0.005)// this threshold is to consider the calculation accuracy, so we can let it larger to gain the robustness
        {
            Logger::error() << "Error in GetRPerp(const Eigen::Matrix3d& R_tar, const Eigen::Matrix3d& R_cur, const Eigen::Vector3d& kappa):   " << tempRot.norm();
            Logger::error() << "||R^T*R-I||=" << tempRot.norm() << ", The input R_tar is not a rotation matrix, please check again!!!!" << Logger::endl();
            return Eigen::Matrix3d::Zero();
        }

        tempRot = Eigen::Matrix3d::Identity() - R_cur.transpose() * R_cur;//tempRot = I - R^T * R, to check whether the input is a rotation matrix
        if(tempRot.norm() > 0.005)// this threshold is to consider the calculation accuracy, so we can let it larger to gain the robustness
        {
            Logger::error() << "Error in GetRPerp(const Eigen::Matrix3d& R_tar, const Eigen::Matrix3d& R_cur, const Eigen::Vector3d& kappa):   " << tempRot.norm();
            Logger::error() << "||R^T*R-I||=" << tempRot.norm() << ", The input R_cur is not a rotation matrix, please check again!!!!" << Logger::endl();
            return Eigen::Matrix3d::Zero();
        }
        /**
         * Input check finished
         * ****/

        Eigen::Matrix3d R_perp;
        Eigen::Vector3d omega_1 = LogRotation( R_tar.transpose() * R_cur );
        Eigen::Vector3d omega_2 = 2 * ( kappa.dot(omega_1) ) * kappa - omega_1;


        Eigen::Vector3d omega_bar = LogRotation( R_cur.transpose() * R_tar * ExpRotation(omega_2) );//symmetry point

        if(omega_bar.norm() < 0.001)
        {
            R_perp = R_cur;
        }
        else{
            Eigen::Matrix3d R_middle = R_tar.transpose() * R_cur * ExpRotation(0.5 * omega_bar);
            if(kappa.cross(LogRotation(R_middle)).norm() < 0.001)//this threshold can not be too small. otherwise the R_perp would jump between the two candidates
            {
                R_perp = R_tar * R_middle;
            }
            else{
                R_perp = R_cur * ExpRotation( -( ( 2*PI - omega_bar.norm() )/( 2*omega_bar.norm() ) ) * omega_bar );
            }
        }
        return R_perp;
    }

    // @author
    //			Gaofeng Li: gaofeng.li@iit.it
    // @data
    //			24-Sep-2019
    // @description
    //			convert the quaternion numbers to a rotation matrix
    // @param
    //			double* qw, qx, qy, qz
    // @return
    //			Eigen::Matrix3d: the rotation matrix
    template <class Number>
    inline Eigen::Matrix3d QuaternionNum2Rot(const Number qw, const Number qx, const Number qy, const Number qz)
    {
        Number qw_, qx_, qy_, qz_, q_norm;
        qw_ = qw;
        qx_ = qx;
        qy_ = qy;
        qz_ = qz;
        q_norm = sqrt(qw_*qw_ + qx_*qx_ + qy_*qy_ + qz_*qz_);
        if (fabs(q_norm - 1) > 0.0001)//check whether it is a unit quaternion
        {
            std::cout << "Warning!!!!---Invalid input in QuaternionNum2Rot(const Number qw, const Number qx, const Number qy, const Number qz)---The input quaternion is not a unit quaternion." << std::endl;
            qw_ = qw_/q_norm;
            qx_ = qx_/q_norm;
            qy_ = qy_/q_norm;
            qz_ = qz_/q_norm;
            std::cout << "We have normalized the input quaternion" << std::endl;
        }
        if(qw_ < 0) // restrict the unit quaternion on the positive half sphere of S^3
        {
            //std::cout << "Debug: on the negtive half sphere" << std::endl;
            qw_ = -qw_;
            qx_ = -qx_;
            qy_ = -qy_;
            qz_ = -qz_;
        }
        return Eigen::Quaterniond(
                qw_,
                qx_,
                qy_,
                qz_).toRotationMatrix();
    }

    // @author
    //			Gaofeng Li: gaofeng.li@iit.it
    // @data
    //			24-Sep-2019
    // @description
    //			convert the quaternion numbers to a rotation matrix
    // @param
    //			double* q: qx, qy, qz, qw
    // @return
    //			Eigen::Matrix3d: the rotation matrix
    template <class Number>
    inline Eigen::Matrix3d QuaternionArray2Rot(const Number* q)
    {
        Number qw_, qx_, qy_, qz_, q_norm;
        qw_ = q[3];
        qx_ = q[0];
        qy_ = q[1];
        qz_ = q[2];
        q_norm = sqrt(qw_*qw_ + qx_*qx_ + qy_*qy_ + qz_*qz_);
        if (fabs(q_norm - 1) > 0.0001)//check whether it is a unit quaternion
        {
            std::cout << "Warning!!!!---Invalid input in QuaternionNum2Rot(const Number qw, const Number qx, const Number qy, const Number qz)---The input quaternion is not a unit quaternion." << std::endl;
            qw_ = qw_/q_norm;
            qx_ = qx_/q_norm;
            qy_ = qy_/q_norm;
            qz_ = qz_/q_norm;
            std::cout << "We have normalized the input quaternion" << std::endl;
        }
        if(qw_ < 0) // restrict the unit quaternion on the positive half sphere of S^3
        {
            //std::cout << "Debug: on the negtive half sphere" << std::endl;
            qw_ = -qw_;
            qx_ = -qx_;
            qy_ = -qy_;
            qz_ = -qz_;
        }
        return Eigen::Quaterniond(
                qw_,
                qx_,
                qy_,
                qz_).toRotationMatrix();
    }

    // @author
    //			Gaofeng Li: gaofeng.li@iit.it
    // @data
    //			18-July-2019
    // @description
    //			convert the quaternion to angle-axis space
    // @param
    //			double* quater: the quaternion need to be converted
    //                         x, y, z, w
    // @return
    //			Eigen::Vector3d: the axis-angle space
    template <class Number>
    inline Eigen::Vector3d Quaternion2Omega(const Number* quater)
    {
        Number q1 = quater[3];//w
        Number q2 = quater[0];//x
        Number q3 = quater[1];//y
        Number q4 = quater[2];//z
        Number wx, wy, wz;

        //normalise
        Number length = sqrt(q1*q1 + q2*q2 + q3*q3 + q4*q4);
        if (fabs(length - 1.0) > 0.0001) {
            q1 = q1 / length;
            q2 = q2 / length;
            q3 = q3 / length;
            q4 = q4 / length;
        }

        // For a quaternion, rotation around k by theta is equal to rotation around -k by -theta. namely, q and -q represent the same rotation
        // we restrict the value: q1 = cos(theta/2) to be positive, such that the angle is restricted by: theta \in [0, PI]
        if(q1 < 0)
        {
            q1 = -q1;
            q2 = -q2;
            q3 = -q3;
            q4 = -q4;
        }

        // To avoid invalid input for acos
        if(q1 - 1 > 0.005)
        {
            std::cout << "Warning!!!!!---Invalid Input in Quaternion2Omega(const Number* quater):  the input for acos is larger than 1---" << std::endl;
            q1 = 1;
        }
        Number angle = 2 * acos(q1);//calculate the angle. Since we have restricted the quaternion to be on the positive half spehre of S^3, the q1 \in [0, 1] ==> acos(q1) \in [0, PI/2] ==> angle \in [0, PI]
        Number norm = sqrt(1 - q1*q1);// the norm for the direction vector is sin(theta/2) = sqrt(1 - cos(theta/2)^2) = sqrt(1 - q1^2),
        if (norm < 0.005)//if the angle == 0, norm = sin(angle/2) will be 0
        {
            wx = q2;
            wy = q3;
            wz = q4;//original: q2, q3, q4
        }
        else
        {
            if (fabs(angle - PI) < 0.0001) // The angle is PI. choose the direction in the positive direction. Here we should use fabs instead of abs. I don't know why. But abs(0.53XXX) return 0. It convert the float value to integer
            {
                if (q2 < 0)//to make sure q2 is positive
                {
                    q2 = -q2;
                    q3 = -q3;
                    q4 = -q4;
                }
                else {
                    if (q2 == 0) // if on the circle that is the intersection of the y-z plane and the sphere
                    {
                        if (q3 < 0)
                        {
                            q3 = -q3;
                            q4 = -q4;
                        }
                        else
                        {
                            if (q3 == 0)//[0 0 -PI] or [0 0 PI]
                            {
                                if (q4 < 0)
                                {
                                    q4 = -q4;
                                }//end of q4 <0;
                            }//end of q3 == 0
                        }//end of q3 < 0
                    }//end of q2 == 0
                }//end of q2 < 0
            }

            wx = angle * (q2 / norm);
            wy = angle * (q3 / norm);
            wz = angle * (q4 / norm);
        }//end of angle ! = 0

        return Eigen::Vector3d(static_cast<double>(wx), static_cast<double>(wy), static_cast<double>(wz));
    }


    // @author
    //			Gaofeng Li: gaofeng.li@iit.it
    // @data
    //			23-Spe-2019
    // @description
    //			convert the Angle-Axis Space to Quaternion
    // @param
    //			Eigen::Vector3d: the axis-angle space vector need to be converted
    // @return
    //			Eigen::Quaternion: the axis-angle space
    inline Eigen::Quaterniond Omega2Quaternion(const Eigen::Vector3d& omega)
    {
        double w, x, y, z;
        double angle_in_radian = omega.norm();
        if(angle_in_radian == 0)
        {
            return Eigen::Quaterniond::Identity();
        }
        else{
            w = cos(angle_in_radian / 2.0);

            x = sin(angle_in_radian / 2.0) * omega.x() / angle_in_radian;
            y = sin(angle_in_radian / 2.0) * omega.y() / angle_in_radian;
            z = sin(angle_in_radian / 2.0) * omega.z() / angle_in_radian;

            Eigen::Quaterniond q_result(w, x, y, z);
            return q_result;
        }

    }

    // @author
    //			Gaofeng Li: gaofeng.li@iit.it
    // @data
    //			23-Spe-2019
    // @description
    //			convert the Rotation matrix to Quaternion
    // @param
    //			Eigen::Matrix3d: the rotation matrix need to be converted
    // @return
    //			Eigen::Quaternion: the axis-angle space
    inline Eigen::Quaterniond Rot2Quaternion(const Eigen::Matrix3d& Rot)
    {
        Eigen::Vector3d omega = LogRotation(Rot);
        //std::cout << "Omega = " << omega << std::endl;
        Eigen::Quaterniond q_result = Omega2Quaternion(omega);
        //std::cout << "Quaternion = [" << q_result.x() << ", " << q_result.y() << ", "  << q_result.z() << ", "  << q_result.w() << "]. "  << std::endl;
        return q_result;
    }

    /// @author
    ///     Gaofeng Li: gaofeng.li@iit.it
    /// @date
    ///     06-Oct-2020
    /// @description
    ///     transform a wrench from the current frame to its target frame
    /// @param
    ///     Eigen::Affine3d& Transformation: the transformation matrix of the target frame relative to the current frame
    ///     Eigen::Vector3d& force_cur: force in current frame
    ///     Eigen::Vector3d& torque_cur: torque in current frame
    ///     Eigen::Vector3d& force_tar: force in target frame
    ///     Eigen::Vector3d& torque_tar: torque in target frame
    inline void WrenchTransformation(const Eigen::Affine3d& Transformation, const Eigen::Vector3d& force_cur, const Eigen::Vector3d& torque_cur, Eigen::Vector3d& force_tar, Eigen::Vector3d& torque_tar)
    {
        Eigen::Matrix3d R = Transformation.linear();
        Eigen::Vector3d p = Transformation.translation();

        Eigen::Matrix3d skew_p = (Eigen::Matrix3d() << 0.0, -p.z(), p.y(),
                p.z(), 0.0, -p.x(),
                -p.y(), p.x(), 0.0).finished();

        force_tar = R.transpose() * force_cur;
        torque_tar = -R.transpose() * skew_p * force_cur + R.transpose() * torque_cur;
    }
}

#endif //myRotPerpH
