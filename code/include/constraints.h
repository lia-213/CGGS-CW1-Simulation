#ifndef CONSTRAINTS_HEADER_FILE
#define CONSTRAINTS_HEADER_FILE

using namespace Eigen;
using namespace std;

typedef enum ConstraintType{DISTANCE, COLLISION} ConstraintType;   //You can expand it for more constraints
typedef enum ConstraintEqualityType{EQUALITY, INEQUALITY} ConstraintEqualityType;

//there is such constraints per two variables that are equal. That is, for every attached vertex there are three such constraints for (x,y,z);
class Constraint{
public:
  
  int m1, m2;                     //Two participating meshes (can be the same)  - auxiliary data for users (constraint class shouldn't use that)
  int v1, v2;                     //Two vertices from the respective meshes - auxiliary data for users (constraint class shouldn't use that)
  double invMass1, invMass2;       //inverse masses of two bodies
  double refValue;                //Reference values to use in the constraint, when needed (like distance)
  bool isUpper;                   //in case this is an inequality constraints, whether it's an upper or a lower bound
  RowVector3d refVector;             //Reference vector when needed (like vector)
  double CRCoeff;                 //velocity bias
  ConstraintType constraintType;  //The type of the constraint, and will affect the value and the gradient. This SHOULD NOT change after initialization!
  ConstraintEqualityType constraintEqualityType;  //whether the constraint is an equality or an inequality
  
  Constraint(const ConstraintType _constraintType, const ConstraintEqualityType _constraintEqualityType, const bool _isUpper, const int& _m1, const int& _v1, const int& _m2, const int& _v2, const double& _invMass1, const double& _invMass2, const RowVector3d& _refVector, const double& _refValue, const double& _CRCoeff):constraintType(_constraintType), constraintEqualityType(_constraintEqualityType), isUpper(_isUpper), m1(_m1), v1(_v1), m2(_m2), v2(_v2), invMass1(_invMass1), invMass2(_invMass2),  refValue(_refValue), CRCoeff(_CRCoeff){
    refVector=_refVector;
  }
  
  ~Constraint(){}
  
  
  
  //computes the impulse needed for all particles to resolve the velocity constraint, and corrects the velocities accordingly.
  //The velocities are a vector (vCOM1, w1, vCOM2, w2) in both input and output.
  //returns true if constraint was already valid with "currVelocities", and false otherwise (false means there was a correction done)
  bool resolve_velocity_constraint(const MatrixXd& currCOMPositions, const MatrixXd& currVertexPositions, const MatrixXd& currCOMVelocities, const MatrixXd& currAngVelocities, const Matrix3d& invInertiaTensor1, const Matrix3d& invInertiaTensor2, MatrixXd& correctedCOMVelocities, MatrixXd& correctedAngVelocities, double tolerance){
    
    
    /***************************TODO: implement this function**********************/
    
    // Extract positions
    RowVector3d com1 = currCOMPositions.row(0);
    RowVector3d com2 = currCOMPositions.row(1);
    RowVector3d p1 = currVertexPositions.row(0);
    RowVector3d p2 = currVertexPositions.row(1);
    
    // Extract velocities
    RowVector3d v1 = currCOMVelocities.row(0);
    RowVector3d w1 = currAngVelocities.row(0);
    RowVector3d v2 = currCOMVelocities.row(1);
    RowVector3d w2 = currAngVelocities.row(1);
    
    // Compute arms from COM to constrained points
    RowVector3d r1 = p1 - com1;
    RowVector3d r2 = p2 - com2;
    
    // Compute constraint direction (gradient of distance constraint)
    RowVector3d diff = p2 - p1;
    double currentDist = diff.norm();
    if (currentDist < 1e-10) {
      // Degenerate case
      correctedCOMVelocities = currCOMVelocities;
      correctedAngVelocities = currAngVelocities;
      return true;
    }
    RowVector3d n = diff / currentDist;
    
    // Compute point velocities at constrained vertices
    RowVector3d vel1 = v1 + w1.cross(r1);
    RowVector3d vel2 = v2 + w2.cross(r2);
    
    // Relative velocity along constraint direction
    double vRel = (vel2 - vel1).dot(n);
    
    // Check if constraint is already satisfied
    if (constraintEqualityType == EQUALITY) {
      if (std::abs(vRel) <= tolerance) {
        correctedCOMVelocities = currCOMVelocities;
        correctedAngVelocities = currAngVelocities;
        return true;
      }
    } else {  // INEQUALITY
      if (isUpper) {
        // Upper bound constraint: velocity should not increase distance beyond bound
        if (vRel <= tolerance) {
          correctedCOMVelocities = currCOMVelocities;
          correctedAngVelocities = currAngVelocities;
          return true;
        }
      } else {
        // Lower bound constraint: velocity should not decrease distance below bound
        if (vRel >= -tolerance) {
          correctedCOMVelocities = currCOMVelocities;
          correctedAngVelocities = currAngVelocities;
          return true;
        }
      }
    }
    
    // Constraint is violated, compute correction
    // Compute the effective mass along the constraint direction
    RowVector3d r1CrossN = r1.cross(n);
    RowVector3d r2CrossN = r2.cross(n);
    
    double denominator = invMass1 + invMass2 +
                        r1CrossN.dot(invInertiaTensor1 * r1CrossN.transpose()) +
                        r2CrossN.dot(invInertiaTensor2 * r2CrossN.transpose());
    
    // Compute Lagrange multiplier (impulse magnitude)
    double lambda = -vRel / denominator;
    
    // Apply velocity corrections using impulse = lambda * n
    RowVector3d impulse = lambda * n;
    
    correctedCOMVelocities.resize(2, 3);
    correctedAngVelocities.resize(2, 3);
    
    // Linear velocity corrections
    correctedCOMVelocities.row(0) = v1 - impulse * invMass1;
    correctedCOMVelocities.row(1) = v2 + impulse * invMass2;
    
    // Angular velocity corrections: delta_w = I^-1 * (r x impulse)
    RowVector3d torque1 = r1.cross(impulse);
    RowVector3d torque2 = r2.cross(-impulse);
    
    correctedAngVelocities.row(0) = w1 + (invInertiaTensor1 * torque1.transpose()).transpose();
    correctedAngVelocities.row(1) = w2 + (invInertiaTensor2 * torque2.transpose()).transpose();
    
    return false;
  }
  
  //projects the position unto the constraint
  //returns true if constraint was already good
  bool resolve_position_constraint(const MatrixXd& currCOMPositions, const MatrixXd& currConstPositions, MatrixXd& correctedCOMPositions, double tolerance){
    
    /***************************TODO: implement this function**********************/
    
    // Extract the two constrained vertex positions
    RowVector3d p1 = currConstPositions.row(0);
    RowVector3d p2 = currConstPositions.row(1);
    
    // Compute current distance between constrained points
    RowVector3d diff = p2 - p1;
    double currentDist = diff.norm();
    
    // Check constraint validity based on type
    double C;  // Constraint function value
    if (constraintEqualityType == EQUALITY) {
      C = currentDist - refValue;
    } else {  // INEQUALITY
      if (isUpper) {
        // Upper bound: C = |p2 - p1| - u*d <= tolerance
        C = currentDist - refValue;
        if (C <= tolerance) {
          correctedCOMPositions = currCOMPositions;
          return true;
        }
      } else {
        // Lower bound: C = |p2 - p1| - l*d >= -tolerance
        C = currentDist - refValue;
        if (C >= -tolerance) {
          correctedCOMPositions = currCOMPositions;
          return true;
        }
      }
    }
    
    // Constraint is violated, need to correct positions
    // Compute the constraint gradient (normalized direction)
    RowVector3d n = diff / currentDist;  // Unit vector from p1 to p2
    
    // Compute correction magnitude using inverse mass weights
    double totalInvMass = invMass1 + invMass2;
    double deltaDist = currentDist - refValue;
    
    // Compute position corrections for each COM
    double w1 = invMass1 / totalInvMass;
    double w2 = invMass2 / totalInvMass;
    
    RowVector3d deltaP1 = w1 * deltaDist * n;
    RowVector3d deltaP2 = -w2 * deltaDist * n;
    
    // Apply corrections to COMs
    correctedCOMPositions.resize(2, 3);
    correctedCOMPositions.row(0) = currCOMPositions.row(0) + deltaP1;
    correctedCOMPositions.row(1) = currCOMPositions.row(1) + deltaP2;
    
    return false;
  }
  
};



#endif /* constraints_h */
