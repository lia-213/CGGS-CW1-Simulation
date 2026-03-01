#ifndef CONSTRAINTS_HEADER_FILE
#define CONSTRAINTS_HEADER_FILE

using namespace Eigen;
using namespace std;

typedef enum ConstraintType{DISTANCE, COLLISION} ConstraintType;
typedef enum ConstraintEqualityType{EQUALITY, INEQUALITY} ConstraintEqualityType;

class Constraint{
public:
  
  int m1, m2;
  int v1, v2;
  double invMass1, invMass2;
  double refValue;
  bool isUpper;
  RowVector3d refVector;
  double CRCoeff;
  ConstraintType constraintType;
  ConstraintEqualityType constraintEqualityType;
  
  Constraint(const ConstraintType _constraintType, const ConstraintEqualityType _constraintEqualityType, const bool _isUpper, const int& _m1, const int& _v1, const int& _m2, const int& _v2, const double& _invMass1, const double& _invMass2, const RowVector3d& _refVector, const double& _refValue, const double& _CRCoeff):constraintType(_constraintType), constraintEqualityType(_constraintEqualityType), isUpper(_isUpper), m1(_m1), v1(_v1), m2(_m2), v2(_v2), invMass1(_invMass1), invMass2(_invMass2),  refValue(_refValue), CRCoeff(_CRCoeff){
    refVector=_refVector;
  }
  
  ~Constraint(){}
  
  bool resolve_velocity_constraint(const MatrixXd& currCOMPositions, const MatrixXd& currVertexPositions, const MatrixXd& currCOMVelocities, const MatrixXd& currAngVelocities, const Matrix3d& invInertiaTensor1, const Matrix3d& invInertiaTensor2, MatrixXd& correctedCOMVelocities, MatrixXd& correctedAngVelocities, double tolerance){
    
    RowVector3d com1 = currCOMPositions.row(0);
    RowVector3d com2 = currCOMPositions.row(1);
    RowVector3d p1 = currVertexPositions.row(0);
    RowVector3d p2 = currVertexPositions.row(1);
    
    RowVector3d v1 = currCOMVelocities.row(0);
    RowVector3d w1 = currAngVelocities.row(0);
    RowVector3d v2 = currCOMVelocities.row(1);
    RowVector3d w2 = currAngVelocities.row(1);
    
    RowVector3d r1 = p1 - com1;
    RowVector3d r2 = p2 - com2;
    
    // n = unit vector from p1 to p2 (direction of increasing distance)
    RowVector3d diff = p2 - p1;
    double currentDist = diff.norm();
    if (currentDist < 1e-10) {
      correctedCOMVelocities = currCOMVelocities;
      correctedAngVelocities = currAngVelocities;
      return true;
    }
    RowVector3d n = diff / currentDist;
    
    // C = |p2 - p1| - refValue
    // Cdot = (vel2 - vel1).dot(n): rate of distance increase (positive = separating)
    RowVector3d vel1 = v1 + w1.cross(r1);
    RowVector3d vel2 = v2 + w2.cross(r2);
    double Cdot = (vel2 - vel1).dot(n);
    
    // The velocity constraint always targets Cdot = 0 (zero relative velocity along n).
    // Return valid only if Cdot is already zero within tolerance.
    if (std::abs(Cdot) <= tolerance) {
      correctedCOMVelocities = currCOMVelocities;
      correctedAngVelocities = currAngVelocities;
      return true;
    }
    
    // Effective mass denominator: J * M_inv * J^T
    RowVector3d r1CrossN = r1.cross(n);
    RowVector3d r2CrossN = r2.cross(n);
    
    double denominator = invMass1 + invMass2 +
                         r1CrossN.dot((invInertiaTensor1 * r1CrossN.transpose()).transpose()) +
                         r2CrossN.dot((invInertiaTensor2 * r2CrossN.transpose()).transpose());
    
    // lambda = -Cdot / denominator
    double lambda = -Cdot / denominator;
    RowVector3d impulse = lambda * n;
    
    correctedCOMVelocities.resize(2, 3);
    correctedAngVelocities.resize(2, 3);
    
    // Body 1: J1 = -n, so dv1 = -invMass1 * impulse
    correctedCOMVelocities.row(0) = v1 - invMass1 * impulse;
    // Body 2: J2 = +n, so dv2 = +invMass2 * impulse
    correctedCOMVelocities.row(1) = v2 + invMass2 * impulse;
    
    // dw1 = -I1_inv * (r1 x impulse)
    correctedAngVelocities.row(0) = w1 - (invInertiaTensor1 * r1.cross(impulse).transpose()).transpose();
    // dw2 = +I2_inv * (r2 x impulse)
    correctedAngVelocities.row(1) = w2 + (invInertiaTensor2 * r2.cross(impulse).transpose()).transpose();
    
    return false;
  }
  
  bool resolve_position_constraint(const MatrixXd& currCOMPositions, const MatrixXd& currConstPositions, MatrixXd& correctedCOMPositions, double tolerance){
    
    RowVector3d p1 = currConstPositions.row(0);
    RowVector3d p2 = currConstPositions.row(1);
    
    RowVector3d diff = p2 - p1;
    double currentDist = diff.norm();
    
    // C = |p2 - p1| - refValue
    double C = currentDist - refValue;
    
    // Check constraint validity
    if (constraintEqualityType == EQUALITY) {
      if (std::abs(C) <= tolerance) {
        correctedCOMPositions = currCOMPositions;
        return true;
      }
    } else {
      if (isUpper) {
        // Cu: distance <= refValue  =>  C <= 0 is valid
        if (C <= tolerance) {
          correctedCOMPositions = currCOMPositions;
          return true;
        }
      } else {
        // Cl: distance >= refValue  =>  C >= 0 is valid
        if (C >= -tolerance) {
          correctedCOMPositions = currCOMPositions;
          return true;
        }
      }
    }
    
    // n = unit vector from p1 to p2
    RowVector3d n = diff / currentDist;
    
    double totalInvMass = invMass1 + invMass2;
    double wt1 = invMass1 / totalInvMass;
    double wt2 = invMass2 / totalInvMass;
    
    // Project COMs by C along n:
    // COM1 += wt1 * C * n  (moves p1 toward p2 when C > 0)
    // COM2 -= wt2 * C * n  (moves p2 toward p1 when C > 0)
    correctedCOMPositions.resize(2, 3);
    correctedCOMPositions.row(0) = currCOMPositions.row(0) + wt1 * C * n;
    correctedCOMPositions.row(1) = currCOMPositions.row(1) - wt2 * C * n;
    
    return false;
  }
  
};

#endif /* constraints_h */
