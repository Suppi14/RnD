/*Copyright  (C)  2009, 2011

Version: 1.0
Author: Ruben Smits, Herman Bruyninckx, Azamat Shakhimardanov
Maintainer: Ruben Smits, Azamat Shakhimardanov
URL: http://www.orocos.org/kdl

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distribute  d in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "my_chainidsolver_vereshchagin.hpp"
#include "kdl/frames_io.hpp"
#include "kdl/utilities/svd_eigen_HH.hpp"
#include <Eigen/SVD>
#include <vector>


namespace my_KDL
{
using namespace Eigen;

ChainIdSolver_Vereshchagin::ChainIdSolver_Vereshchagin(const Chain& chain_, Twist root_acc, unsigned int _nc) :
    chain(chain_), nj(chain.getNrOfJoints()), ns(chain.getNrOfSegments()), nc(_nc),
    results(ns + 1, segment_info(nc))
{
    acc_root = root_acc;
    // std::cout << "No of constraints" << nc << std::endl;
    //Provide the necessary memory for computing the inverse of M0
    nu_sum.resize(nc);
    M_0_inverse.resize(nc, nc);
    Um = MatrixXd::Identity(nc, nc);
    Vm = MatrixXd::Identity(nc, nc);
    Sm = VectorXd::Ones(nc);
    tmpm = VectorXd::Ones(nc);
}

void ChainIdSolver_Vereshchagin::updateInternalDataStructures() {
    ns = chain.getNrOfSegments();
    results.resize(ns+1,segment_info(nc));
}

int ChainIdSolver_Vereshchagin::CartToJnt(const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Jacobian& alfa, const JntArray& beta, const Wrenches& f_ext, JntArray &torques)
{
    nj = chain.getNrOfJoints();
    if(ns != chain.getNrOfSegments())
        return (error = E_NOT_UP_TO_DATE);
    //Check sizes always
    if (q.rows() != nj || q_dot.rows() != nj || q_dotdot.rows() != nj || torques.rows() != nj || f_ext.size() != ns)
        return (error = E_SIZE_MISMATCH);
    if (alfa.columns() != nc || beta.rows() != nc)
        return (error = E_SIZE_MISMATCH);
    //do an upward recursion for position and velocities
    this->initial_upwards_sweep(q, q_dot, q_dotdot, f_ext);
    //do an inward recursion for inertia, forces and constraints
    this->downwards_sweep(alfa, torques);
    //Solve for the constraint forces
    this->constraint_calculation(beta,alfa);
    //do an upward recursion to propagate the result
    this->final_upwards_sweep(q_dotdot, torques);
    return (error = E_NOERROR);
}

void ChainIdSolver_Vereshchagin::initial_upwards_sweep(const JntArray &q, const JntArray &qdot, const JntArray &qdotdot, const Wrenches& f_ext)
{
    //if (q.rows() != nj || qdot.rows() != nj || qdotdot.rows() != nj || f_ext.size() != ns)
    //        return -1;
    // std::cout<< "Vereschagin here" <<std::endl;
    unsigned int j = 0;
    F_total = Frame::Identity();
    for (unsigned int i = 0; i < ns; i++)
    {
        //Express everything in the segments reference frame (body coordinates)
        //which is at the segments tip, i.e. where the next joint is attached.

        //Calculate segment properties: X,S,vj,cj
        const Segment& segment = chain.getSegment(i);
        segment_info& s = results[i + 1];
        //The pose between the joint root and the segment tip (tip expressed in joint root coordinates)
        s.F = segment.pose(q(j)); //X pose of each link in link coord system

        F_total = F_total * s.F; //X pose of the each link in root coord system
        s.F_base = F_total; //X pose of the each link in root coord system for getter functions

        //The velocity due to the joint motion of the segment expressed in the segments reference frame (tip)
        Twist vj = s.F.M.Inverse(segment.twist(q(j), qdot(j))); //XDot of each link
        //Twist aj = s.F.M.Inverse(segment.twist(q(j), qdotdot(j))); //XDotDot of each link
        // std::cout<<"Twist" << vj<<std::endl;
        //The unit velocity due to the joint motion of the segment expressed in the segments reference frame (tip)
        s.Z = s.F.M.Inverse(segment.twist(q(j), 1.0));
        //Put Z in the joint root reference frame:
        s.Z = s.F * s.Z;

        //The total velocity of the segment expressed in the segments reference frame (tip)
        if (i != 0)
        {
            s.v = s.F.Inverse(results[i].v) + vj; // recursive velocity of each link in segment frame
            //s.A=s.F.Inverse(results[i].A)+aj;
            s.A = s.F.M.Inverse(results[i].A);
        }
        else
        {
            s.v = vj;
            s.A = s.F.M.Inverse(acc_root);
        }
        //c[i] = cj + v[i]xvj (remark: cj=0, since our S is not time dependent in local coordinates)
        //The velocity product acceleration
        //std::cout << i << " Initial upward" << s.v << std::endl;
        s.C = s.v*vj; //This is a cross product: cartesian space BIAS acceleration in local link coord.
        //Put C in the joint root reference frame
        s.C = s.F * s.C; //+F_total.M.Inverse(acc_root));
        //The rigid body inertia of the segment, expressed in the segments reference frame (tip)
        s.H = segment.getInertia();

        //wrench of the rigid body bias forces and the external forces on the segment (in body coordinates, tip)
        //external forces are taken into account through s.U.
        Wrench FextLocal = F_total.M.Inverse() * f_ext[i];
        s.U = s.v * (s.H * s.v) - FextLocal; //f_ext[i];
        if (segment.getJoint().getType() != Joint::None)
            j++;
    }

}

void ChainIdSolver_Vereshchagin::downwards_sweep(const Jacobian& alfa, const JntArray &torques)
{
    int j = nj - 1;
    for (int i = ns; i >= 0; i--)
    {
        //Get a handle for the segment we are working on.
        segment_info& s = results[i];
        //For segment N,
        //tilde is in the segment refframe (tip, not joint root)
        //without tilde is at the joint root (the childs tip!!!)
        //P_tilde is the articulated body inertia
        //R_tilde is the sum of external and coriolis/centrifugal forces
        //M is the (unit) acceleration energy already generated at link i
        //G is the (unit) magnitude of the constraint forces at link i
        //E are the (unit) constraint forces due to the constraints
        if (i == (int)ns)
        {
            s.P_tilde = s.H;
            s.R_tilde = s.U;
            s.Projection_tilde = Matrix6d::Identity();
            s.M.setZero();
            s.G.setZero();
            //changeBase(alfa_N,F_total.M.Inverse(),alfa_N2);
            for (unsigned int r = 0; r < 3; r++)
                for (unsigned int c = 0; c < nc; c++)
                {
                    //copy alfa constrain force matrix in E~
                    s.E_tilde(r, c) = alfa(r+3, c);
                    s.E_tilde(r + 3, c) = alfa(r, c);
                }
              // std::cout << "Alpha" << s.E_tilde << std::endl;
            //Change the reference frame of alfa to the segmentN tip frame
            //F_Total holds end effector frame, if done per segment bases then constraints could be extended to all segments
            Rotation base_to_end = F_total.M.Inverse();
            for (unsigned int c = 0; c < nc; c++)
            {
                Wrench col(Vector(s.E_tilde(3, c), s.E_tilde(4, c), s.E_tilde(5, c)),
                           Vector(s.E_tilde(0, c), s.E_tilde(1, c), s.E_tilde(2, c)));
                col = base_to_end*col;
                s.E_tilde.col(c) << Vector3d::Map(col.torque.data), Vector3d::Map(col.force.data);
            }
        }

        else
        {
            //For all others:
            //Everything should expressed in the body coordinates of segment i
            segment_info& child = results[i + 1];
            //Copy PZ into a vector so we can do matrix manipulations, put torques above forces
            Vector6d vPZ;
            vPZ << Vector3d::Map(child.PZ.torque.data), Vector3d::Map(child.PZ.force.data);
            Matrix6d PZDPZt;
            PZDPZt.noalias() = vPZ * vPZ.transpose();
            PZDPZt /= child.D;
            /**/
            Vector6d vZ;//
            vZ << Vector3d::Map(child.Z.rot.data), Vector3d::Map(child.Z.vel.data);
            s.Projection_tilde = Matrix6d::Identity() - (vPZ * vZ.transpose() / child.D);
            if (i != (int)(ns)-1)
            {
                s.Projection_tilde *= child.Projection;
            }
            // std::cout << "****Projection_tilde matrix**** " << std::endl << s.Projection_tilde << std::endl;
            Eigen::BDCSVD<Matrix6d> svd(s.Projection_tilde);
            // std::cout << "****Rank of Projection_tilde**** " << svd.rank() << std::endl;

            //equation a) (see Vereshchagin89) PZDPZt=[I,H;H',M]
            //Azamat:articulated body inertia as in Featherstone (7.19)
            s.P_tilde = s.H + child.P - ArticulatedBodyInertia(PZDPZt.bottomRightCorner<3,3>(), PZDPZt.topRightCorner<3,3>(), PZDPZt.topLeftCorner<3,3>());
        //    s.P_tilde = s.H + child.P - ArticulatedBodyInertia(PZDPZt.bottomRightCorner<3,3>(), PZDPZt.topRightCorner<3,3>(), PZDPZt.topLeftCorner<3,3>());
            //equation b) (see Vereshchagin89)
            //Azamat: bias force as in Featherstone (7.20)
            s.R_tilde = s.U + child.R + child.PC + (child.PZ / child.D) * child.u;
            //equation c) (see Vereshchagin89)
            s.E_tilde = child.E;

            //Azamat: equation (c) right side term
            s.E_tilde.noalias() -= (vPZ * child.EZ.transpose()) / child.D;


            //equation d) (see Vereshchagin89)
            s.M = child.M;
            //Azamat: equation (d) right side term
            s.M.noalias() -= (child.EZ * child.EZ.transpose()) / child.D;
            // std::cout<<"EZ matrix" <<child.D << std::endl;

            // std::cout << "****M matrix in downwards_sweep**** " <<std::endl << s.M<< std::endl;
            //equation e) (see Vereshchagin89)
            s.G = child.G;
            Twist CiZDu = child.C + (child.Z / child.D) * child.u;
            Vector6d vCiZDu;
            vCiZDu << Vector3d::Map(CiZDu.rot.data), Vector3d::Map(CiZDu.vel.data);
            s.G.noalias() += child.E.transpose() * vCiZDu;
        }
        if (i != 0)
        {
            //Transform all results to joint root coordinates of segment i (== body coordinates segment i-1)
            //equation a)
            s.P = s.F * s.P_tilde;
            //Projection transformation
            Matrix6d force_transform ;
            Matrix3d rx;
            force_transform.setZero();
            force_transform.block(0,0,3,3) = Map<Eigen::Matrix3d>(s.F.M.data,3,3);
            force_transform.block(3,3,3,3) = Map<Eigen::Matrix3d>(s.F.M.data,3,3);
            rx << 0, -s.F.p.z(), s.F.p.y(),
                  s.F.p.z(), 0, -s.F.p.x(),
                - s.F.p.y(), s.F.p.x(), 0 ;
            // std::cout << "rx" << rx<< std::endl;
            force_transform.block(0,3,3,3) = rx * Map<Eigen::Matrix3d>(s.F.M.data,3,3);
            s.Projection = force_transform * s.Projection_tilde;
            Eigen::BDCSVD<Matrix6d> svd(s.Projection);
            // std::cout << "****Projection****" << std::endl << s.Projection << std::endl;
            // std::cout << "****Rank of Projection****  " << svd.rank() << std::endl;
            //equation b)
            s.R = s.F * s.R_tilde;
            //equation c), in matrix: torques above forces, so switch and switch back
            for (unsigned int c = 0; c < nc; c++)
            {
                Wrench col(Vector(s.E_tilde(3, c), s.E_tilde(4, c), s.E_tilde(5, c)),
                           Vector(s.E_tilde(0, c), s.E_tilde(1, c), s.E_tilde(2, c)));
                col = s.F*col;
                s.E.col(c) << Vector3d::Map(col.torque.data), Vector3d::Map(col.force.data);
            }

            //needed for next recursion
            s.PZ = s.P * s.Z;
            s.D = dot(s.Z, s.PZ);
            s.PC = s.P * s.C;

            //u=(Q-Z(R+PC)=sum of external forces along the joint axes,
            //R are the forces coming from the children,
            //Q is taken zero (do we need to take the previous calculated torques?

            //projection of coriolis and centrepital forces into joint subspace (0 0 Z)
            s.totalBias = -dot(s.Z, s.R + s.PC);
            s.u = torques(j) + s.totalBias;

            //Matrix form of Z, put rotations above translations
            Vector6d vZ;
            vZ << Vector3d::Map(s.Z.rot.data), Vector3d::Map(s.Z.vel.data);
            s.EZ.noalias() = s.E.transpose() * vZ;

            if (chain.getSegment(i - 1).getJoint().getType() != Joint::None)
                j--;
        }
    }
}

void ChainIdSolver_Vereshchagin::constraint_calculation(const JntArray& beta,const Jacobian& alfa)
{
    //equation f) nu = M_0_inverse*(beta_N - E0_tilde`*acc0 - G0)
    //M_0_inverse, always nc*nc symmetric matrix
    //std::cout<<"M0: in constraint calc"<<results[0].M<<std::endl;
    // std::cout<<"M0 transpose"<< results[0].M.transpose() << std::endl;//To verify its symmetricity
    //results[0].M-=MatrixXd::Identity(nc,nc);
    //std::cout<<"augmented M0: "<<results[0].M<<std::endl;

    //ToDo: Need to check ill conditions
    int nz = 0;
    for(int r=0 ; r<alfa.rows(); r++){
      for(int c=0; c<alfa.columns(); c++){
        if(alfa.operator()(r,c)!=0){
          nz = nz+1;
        }
      }
    }
    std::cout<<"***No of non zero rows i.e Rank of matrix alpha is*** " << nz << std::endl;
    // M_0_inverse=results[0].M.inverse();
    // std::cout<<"M0_inverse"<<results[0].M.inverse()<<std::endl;

    svd_eigen_HH(results[0].M, Um, Sm, Vm, tmpm);
    // std::cout<<"Um"<<Um<<std::endl;
     std::cout<<"Sm"<<Sm<<std::endl;
    // std::cout<<"Vm"<<Vm<<std::endl;
    //check singularity for task specification
    int nz_sigma=0;
    for(int r=0 ; r<Sm.rows(); r++){
      for(int c=0; c<Sm.cols(); c++){
        if(fabs(Sm.operator()(r,c))>= 1e-4){
          nz_sigma = nz_sigma+1;
        }
      }
    }
    std::cout<<" ***No of non zero rows i.e Rank of matrix sigma is*** " << nz_sigma << std::endl;
    if(nz_sigma<nz){
      std::cout << "****The kinematic chain cannot satisfy this constraint for the current task specification****" << std:: endl;
    }

    int counter=0;
    // int undentifiable = 0;
  //  epsilon = std::numeric_limits<double>::epsilon();
     epsilon = 1e-4; //value close to zero
    // std::cout << "epsilon value" << epsilon << std::endl;
    // int counter =0;
    for(int r=0;r<Sm.rows();r++){
    //  for(int c=0;c<Sm.cols(); c++){
        if(fabs(Sm(r)) >= epsilon){
        counter+=1;
        }
//      }
    }

    std::string direction[6] = {"linear_x","linear_y","linear_z","angular_x","angular_y","angular_z"};
    std::vector<int> non_zero_values;
    Vm_block = Vm.block(0,0,Vm.rows(),counter);
    std::cout<< " ****blocking Vm****   "<< std::endl<< Vm_block << std::endl;
    for(int r=0; r<Vm_block.rows(); r++){
      int nZV = 0;
      for(int c=0; c<Vm_block.cols(); c++){
          if(fabs(Vm_block.operator()(r,c)) >= 1e-4){
              nZV++;
          }

      }
      non_zero_values.push_back(nZV);
    }
    // if(non_zero_values !=0){
      for (int i=0; i< non_zero_values.size(); i++){
        if(non_zero_values[i] == 0){
          std::cout<< " ***Task not achievable in " << "direction " << direction[i] << std::endl;
        }
        else{
          std::cout<< " ***Task achievable in " << "direction " << direction[i] << " with nZV*** " << non_zero_values[i]<< std::endl;
        }

      }


     // std::cout<<"Number of non zero column elements ******" << non_zero_values << std::endl;
    // std::cout<< "unidentifiable" << undentifiable<< std::endl;
    // std::cout<< "Sm values counter" << counter << std::endl;
    // Sm >= epsilon;
//     for(int r=0 ; r< 6; r++){
//       for(int c=0; c< 6; c++){
//         bool var = Sm  >= epsilon;
//   }
// }
    // nullspace_calc(results[0].M);
    // std::cout<< "sm epsilon values" << var << std::endl;
  //   for(int r=0;r<6;r++){
  //     for(int c=0; c<6; c++){
  //       if(Sm.operator() >= std::numeric_limits<double>::epsilon())
  //       // std::cout<<Sm.operator()(r,c);
  //   }
  // }

     // std::cout<<"Vm transpose" << Vm.transpose()<<std::endl;
    //truncated svd, what would sdls, dls physically mean?
    for(unsigned int i = 0; i < nc; i++)
        if (Sm(i) < 1e-14)
            Sm(i) = 0.0;
        else
            Sm(i) = 1 / Sm(i);
    // std::cout<<"Sm"<<Sm<<std::endl;
    //To check for singularity
    results[0].M.noalias() = Vm * Sm.asDiagonal();
    M_0_inverse.noalias() = results[0].M * Um.transpose();
    //results[0].M.ldlt().solve(MatrixXd::Identity(nc,nc),&M_0_inverse);
    //results[0].M.computeInverse(&M_0_inverse);
    Vector6d acc;
    acc << Vector3d::Map(acc_root.rot.data), Vector3d::Map(acc_root.vel.data);
    nu_sum.noalias() = -(results[0].E_tilde.transpose() * acc);
    //nu_sum.setZero();
    nu_sum += beta.data;
    nu_sum -= results[0].G;

    //equation f) nu = M_0_inverse*(beta_N - E0_tilde`*acc0 - G0)
    nu.noalias() = M_0_inverse * nu_sum;
    // std::cout << "****M0**** " << std::endl<< results[0].M <<std::endl;
}

void ChainIdSolver_Vereshchagin::final_upwards_sweep(JntArray &q_dotdot, JntArray &torques)
{
    unsigned int j = 0;

    for (unsigned int i = 1; i <= ns; i++)
    {
        segment_info& s = results[i];
        //Calculation of joint and segment accelerations
        //equation g) qdotdot[i] = D^-1*(Q - Z'(R + P(C + acc[i-1]) + E*nu))
        // = D^-1(u - Z'(P*acc[i-1] + E*nu)
        Twist a_g;
        Twist a_p;
        if (i == 1)
        {
            a_p = acc_root;
        }
        else
        {
            a_p = results[i - 1].acc;
        }

        //The contribution of the constraint forces at segment i
        Vector6d tmp = s.E*nu;
        Wrench constraint_force = Wrench(Vector(tmp(3), tmp(4), tmp(5)),
                                         Vector(tmp(0), tmp(1), tmp(2)));

        //acceleration components are also computed
        //Contribution of the acceleration of the parent (i-1)
        Wrench parent_force = s.P*a_p;
        double parent_forceProjection = -dot(s.Z, parent_force);
        double parentAccComp = parent_forceProjection / s.D;

        //The constraint force and acceleration force projected on the joint axes -> axis torque/force
        double constraint_torque = -dot(s.Z, constraint_force);
        //The result should be the torque at this joint

        torques(j) = constraint_torque;
        //s.constAccComp = torques(j) / s.D;
        s.constAccComp = constraint_torque / s.D;
        s.nullspaceAccComp = s.u / s.D;
        // external_forces = s.u;
        std::cout<<"external force torque" << s.u << std::endl;
        //total joint space acceleration resulting from accelerations of parent joints, constraint forces and
        // nullspace forces.
        q_dotdot(j) = (s.nullspaceAccComp + parentAccComp + s.constAccComp);
        s.acc = s.F.Inverse(a_p + s.Z * q_dotdot(j) + s.C);//returns acceleration in link distal tip coordinates. For use needs to be transformed
        // std::cout << "***Acceleration*** " << s.acc << std::endl;
        if (chain.getSegment(i - 1).getJoint().getType() != Joint::None)
            j++;

    }

}

void ChainIdSolver_Vereshchagin::get_transformed_link_acceleration(Twists& xDotDot)
{
    //Assersions need to be replaced with run-time errors!
    //For example errors specified in SolverI class
    //Because KDL compiles in RELEASE mode!
    assert(xDotDot.size() == ns + 1);
    xDotDot[0] = acc_root;
    for (int i = 1; i < ns + 1; i++) {
        xDotDot[i] = results[i].F_base.M * results[i].acc;
    }
}


/*
void ChainIdSolver_Vereshchagin::getLinkCartesianPose(Frames& x_base)
{
    for (int i = 0; i < ns; i++)
    {
        x_base[i] = results[i + 1].F_base;
    }
    return;
}

void ChainIdSolver_Vereshchagin::getLinkCartesianVelocity(Twists& xDot_base)
{

    for (int i = 0; i < ns; i++)
    {
        xDot_base[i] = results[i + 1].F_base.M * results[i + 1].v;
    }
    return;
}

void ChainIdSolver_Vereshchagin::getLinkCartesianAcceleration(Twists& xDotDot_base)
{

    for (int i = 0; i < ns; i++)
    {
        xDotDot_base[i] = results[i + 1].F_base.M * results[i + 1].acc;
        //std::cout << "XDotDot_base[i] " << xDotDot_base[i] << std::endl;
    }
    return;
}

void ChainIdSolver_Vereshchagin::getLinkPose(Frames& x_local)
{
    for (int i = 0; i < ns; i++)
    {
        x_local[i] = results[i + 1].F;
    }
    return;
}

void ChainIdSolver_Vereshchagin::getLinkVelocity(Twists& xDot_local)
{
    for (int i = 0; i < ns; i++)
    {
        xDot_local[i] = results[i + 1].v;
    }
    return;

}

void ChainIdSolver_Vereshchagin::getLinkAcceleration(Twists& xDotdot_local)
{
     for (int i = 0; i < ns; i++)
    {
        xDotdot_local[i] = results[i + 1].acc;
    }
    return;

}

void ChainIdSolver_Vereshchagin::getJointBiasAcceleration(JntArray& bias_q_dotdot)
{
    for (int i = 0; i < ns; i++)
    {
        //this is only force
        double tmp = results[i + 1].totalBias;
        //this is acceleration
        bias_q_dotdot(i) = tmp / results[i + 1].D;

        //s.totalBias = - dot(s.Z, s.R + s.PC);
        //std::cout << "totalBiasAccComponent" << i << ": " << bias_q_dotdot(i) << std::endl;
        //bias_q_dotdot(i) = s.totalBias/s.D

    }
    return;

}

void ChainIdSolver_Vereshchagin::getJointConstraintAcceleration(JntArray& constraint_q_dotdot)
{
    for (int i = 0; i < ns; i++)
    {
        constraint_q_dotdot(i) = results[i + 1].constAccComp;
        //double tmp = results[i + 1].u;
        //s.u = torques(j) + s.totalBias;
        // std::cout << "s.constraintAccComp" << i << ": " << results[i+1].constAccComp << std::endl;
        //nullspace_q_dotdot(i) = s.u/s.D

    }
    return;


}

//Check the name it does not seem to be appropriate

void ChainIdSolver_Vereshchagin::getJointNullSpaceAcceleration(JntArray& nullspace_q_dotdot)
{
    for (int i = 0; i < ns; i++)
    {
        nullspace_q_dotdot(i) = results[i + 1].nullspaceAccComp;
        //double tmp = results[i + 1].u;
        //s.u = torques(j) + s.totalBias;
        //std::cout << "s.nullSpaceAccComp" << i << ": " << results[i+1].nullspaceAccComp << std::endl;
        //nullspace_q_dotdot(i) = s.u/s.D

    }
    return;


}

//This is not only a bias force energy but also includes generalized forces
//change type of parameter G
//this method should return array of G's

void ChainIdSolver_Vereshchagin::getLinkBiasForceAcceleratoinEnergy(Eigen::VectorXd& G)
{
    for (int i = 0; i < ns; i++)
    {
        G = results[i + 1].G;
        //double tmp = results[i + 1].u;
        //s.u = torques(j) + s.totalBias;
        //std::cout << "s.G " << i << ":  " << results[i+1].G << std::endl;
        //nullspace_q_dotdot(i) = s.u/s.D

    }
    return;

}

//this method should return array of R's

void ChainIdSolver_Vereshchagin::getLinkBiasForceMatrix(Wrenches& R_tilde)
{
    for (int i = 0; i < ns; i++)
    {
        R_tilde[i] = results[i + 1].R_tilde;
        //Azamat: bias force as in Featherstone (7.20)
        //s.R_tilde = s.U + child.R + child.PC + child.PZ / child.D * child.u;
        std::cout << "s.R_tilde " << i << ":  " << results[i + 1].R_tilde << std::endl;
    }
    return;
}

*/
// void nullspace_calc(Eigen::MatrixXd& M){
//
// }

}//namespace
