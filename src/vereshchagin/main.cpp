/*
Author(s): Supriya Vadiraj, Sven Schneider
Copyright (c) [2018]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#include <iostream>
#include <kdl/chain.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/frames_io.hpp>
#include <kdl/chainidsolver_recursive_newton_euler.hpp>
#include <kdl/kinfam_io.hpp>
#include <my_chainidsolver_vereshchagin.hpp>
#include <kdl/chainjnttojacsolver.hpp>
#include <kdl/jacobian.hpp>
#include <Eigen/SVD>

const int NUMBER_OF_JOINTS = 7;
const int NUMBER_OF_SEGMENTS = 7 ;
const int NUMBER_OF_CONSTRAINTS = 6;

class my_robot
{
  public:
      KDL::Chain my_bot;
      void solve_fk(KDL::Chain &r);
      void solve_jac(KDL::Chain &r);
};

class motion_specification
{
    public:
        motion_specification(int nj,int ns, int nc) : q(nj), q_dot(nj), q_dot_dot(nj), feedforward_torque(nj),
                                                      alpha(nc),
                                                      beta(nc),
                                                      external_force(ns),jac(nj){
        }
        KDL::JntArray q;
        KDL::JntArray q_dot;
        KDL::JntArray q_dot_dot;
        KDL::JntArray feedforward_torque;
        KDL::Jacobian alpha;
        KDL::JntArray beta;
        KDL::Wrenches external_force;
        KDL::Jacobian jac;
};

//Jacobain

// void my_robot::solve_jac(KDL::Chain &r){
//   KDL::ChainJntToJacSolver jacsolver(r);
//   KDL::JntArray q(r.getNrOfJoints());
//   q(0) = M_PI / 2.0;
//   q(1) = M_PI / 2.0;
//   q(2) = M_PI / 1.0;
//   q(3) = -M_PI / 2.0;
//   KDL::Jacobian jac(r.getNrOfJoints());
//   jacsolver.JntToJac(q,jac);
//   std::cout << "Traditional singularity "<< jac << std::endl;
//   Eigen::BDCSVD<Eigen::Matrix<double,6,Eigen::Dynamic> > svd(jac.data);
//   std::cout << "Rank of Jacobain " << svd.rank() << std::endl;
// }

// recursive forward kinematics solver
void my_robot::solve_fk(KDL::Chain &r){
  KDL::ChainFkSolverPos_recursive fksolver(r);
  KDL::JntArray q(r.getNrOfJoints());
   q(0) = 0;
   q(1) = 0;
   q(2) = 0;
   q(3) = 0;
   KDL::Frame end_effector_pose;
   fksolver.JntToCart(q, end_effector_pose);
   std::cout<<" FKSolver: "<< end_effector_pose << std::endl;
}


void create_motion_specification(motion_specification &m)
{
  int number_of_joints = 7;
  int number_of_segments = m.external_force.size();

  //specifying q, qd, qdd
  m.q(0) = 0.0;
  m.q(1) = 0.0;
  m.q(2) = 0.0;
  m.q(3) = 0.0;
  m.q(4) = 0.0;
  m.q(5) = 0.0;
  m.q(6) = 0.0;


  for(int i=0;i<number_of_joints;i++){
    // m.q(i) = 0.0;
    m.q_dot(i) = 0.0;
    m.q_dot_dot(i) = 0.0;
  }

  m.jac.resize(number_of_joints);

//specifying torques
  for(int i=0;i < number_of_joints; i++){
    m.feedforward_torque(i) = 0.0;
  }

  m.external_force[6] = KDL::Wrench(
                            KDL::Vector(0.0, 0.0, 0.0), //Force
                            KDL::Vector(0.0, 0.0, 1.0));
  // for (int i = 0; i < number_of_segments - 1; i++) {
  //         KDL::Wrench externalForce(
  //             KDL::Vector(15.0, 1.0, 7.0), //Force
  //             KDL::Vector(0.0, 0.0, 0.0)); //Torque
  //         m.external_force[i] = externalForce;
  // }

      KDL::Twist alfa_x(
              KDL::Vector(1.0, 0.0, 0.0),     // linear
              KDL::Vector(0.0, 0.0, 0.0));    // angular
      m.alpha.setColumn(0, alfa_x);
      m.beta(0) = 2.0;

      KDL::Twist alfa_y(
              KDL::Vector(0.0, 1.0, 0.0),     // linear
              KDL::Vector(0.0, 0.0, 0.0));    // angular
      m.alpha.setColumn(1, alfa_y);
      m.beta(1) = 2.0;

      KDL::Twist alfa_z(
              KDL::Vector(0.0, 0.0, 1.0),     // linear
              KDL::Vector(0.0, 0.0, 0.0));    // angular
      m.alpha.setColumn(2, alfa_z);
      m.beta(2) = 2.0;
      //
      KDL::Twist alfa_xr(
              KDL::Vector(0.0, 0.0, 0.0),     // linear
              KDL::Vector(1.0, 0.0, 0.0));    // angular
      m.alpha.setColumn(3, alfa_xr);
      m.beta(3) = 2.0;

      KDL::Twist alfa_yr(
              KDL::Vector(0.0, 0.0, 0.0),     // linear
              KDL::Vector(0.0, 1.0, 0.0));    // angular
      m.alpha.setColumn(4, alfa_yr);
      m.beta(4) = 2.0;

      KDL::Twist alfa_zr(
              KDL::Vector(0.0, 0.0, 0.0),     // linear
              KDL::Vector(0.0, 0.0, 1.0));    // angular
      m.alpha.setColumn(5, alfa_zr);
      m.beta(5) = 2.0;
  }


//main function
int main(int argc, char **argv)
{
  my_robot m;
  motion_specification my_motion(NUMBER_OF_JOINTS, NUMBER_OF_SEGMENTS, NUMBER_OF_CONSTRAINTS);

  //	joint 0
  	// m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::None),
  	// 			  KDL::Frame::Frame::DH_Craig1989(0.0, 0.0, 0.31, 0.0)
  	// 			  ));
    /////////////////////////////////////////////////////////////////////////
//  	joint 1
  	m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
  				  KDL::Frame::Frame::DH_Craig1989(0.0, 1.5707963, 0.0, 0.0),
  				  KDL::Frame::Frame::DH_Craig1989(0.0, 1.5707963, 0.0, 0.0).Inverse()*KDL::RigidBodyInertia(2,
  												 KDL::Vector::Zero(),
  												 KDL::RotationalInertia(0.0,0.0,0.0115343,0.0,0.0,0.0))));

  	//joint 2
  	m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
  				  KDL::Frame::Frame::DH_Craig1989(0.0, -1.5707963, 0.4, 0.0),
  				  KDL::Frame::Frame::DH_Craig1989(0.0, -1.5707963, 0.4, 0.0).Inverse()*KDL::RigidBodyInertia(2,
  												   KDL::Vector(0.0,-0.3120511,-0.0038871),
  												   KDL::RotationalInertia(-0.5471572,-0.0000302,-0.5423253,0.0,0.0,0.0018828))));

  	//joint 3
  	m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
  				  KDL::Frame::Frame::DH_Craig1989(0.0, -1.5707963, 0.0, 0.0),
  				  KDL::Frame::Frame::DH_Craig1989(0.0, -1.5707963, 0.0, 0.0).Inverse()*KDL::RigidBodyInertia(2,
  												   KDL::Vector(0.0,-0.0015515,0.0),
  												   KDL::RotationalInertia(0.0063507,0.0,0.0107804,0.0,0.0,-0.0005147))));

  	//joint 4
  	m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
  				  KDL::Frame::Frame::DH_Craig1989(0.0, 1.5707963, 0.39, 0.0),
  				  KDL::Frame::Frame::DH_Craig1989(0.0, 1.5707963, 0.39, 0.0).Inverse()*KDL::RigidBodyInertia(2,
  												   KDL::Vector(0.0,0.5216809,0.0),
  												   KDL::RotationalInertia(-1.0436952,0.0,-1.0392780,0.0,0.0,0.0005324))));

  	//joint 5
  	m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
  				  KDL::Frame::Frame::DH_Craig1989(0.0, 1.5707963, 0.0, 0.0),
  				  KDL::Frame::Frame::DH_Craig1989(0.0, 1.5707963, 0.0, 0.0).Inverse()*KDL::RigidBodyInertia(2,
  												   KDL::Vector(0.0,0.0119891,0.0),
  												   KDL::RotationalInertia(0.0036654,0.0,0.0060429,0.0,0.0,0.0004226))));

  	//joint 6
  	m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
  				  KDL::Frame::Frame::DH_Craig1989(0.0, -1.5707963, 0.0, 0.0),
  				  KDL::Frame::Frame::DH_Craig1989(0.0, -1.5707963, 0.0, 0.0).Inverse()*KDL::RigidBodyInertia(2,
  												   KDL::Vector(0.0,0.0080787,0.0),
  												   KDL::RotationalInertia(0.0010431,0.0,0.0036376,0.0,0.0,0.0000101))));
  	//joint 7
  	m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
  				   KDL::Frame::Frame::Identity(),
  				   KDL::RigidBodyInertia(2,
  												   KDL::Vector::Zero(),
  												   KDL::RotationalInertia(0.000001,0.0,0.0001203,0.0,0.0,0.0))));

///////////////////////////////////////////////////////
  // m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
  // 				  KDL::Frame::DH_Craig1989(0.0, 1.5707963, 0.0, 0.0),
  // 				  KDL::Frame::DH_Craig1989(0.0, 1.5707963, 0.0, 0.0).Inverse()*KDL::RigidBodyInertia(2,
  // 										 KDL::Vector::Zero(),
  // 										 KDL::RotationalInertia(0.0,0.0,0.0115343,0.0,0.0,0.0))));
  //
  // 	//joint 2
  // 	m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
  // 				  KDL::Frame::DH_Craig1989(0.0, -1.5707963, 0.4, 0.0),
  // 				  KDL::Frame::DH_Craig1989(0.0, -1.5707963, 0.4, 0.0).Inverse()*KDL::RigidBodyInertia(2,
  // 										   KDL::Vector(0.0,-0.3120511,-0.0038871),
  // 										   KDL::RotationalInertia(-0.5471572,-0.0000302,-0.5423253,0.0,0.0,0.0018828))));
  //
  // 	//joint 3
  // 	m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
  // 				  KDL::Frame::DH_Craig1989(0.0, -1.5707963, 0.0, 0.0),
  // 				  KDL::Frame::DH_Craig1989(0.0, -1.5707963, 0.0, 0.0).Inverse()*KDL::RigidBodyInertia(2,
  // 										   KDL::Vector(0.0,-0.0015515,0.0),
  // 										   KDL::RotationalInertia(0.0063507,0.0,0.0107804,0.0,0.0,-0.0005147))));
  //
  // 	//joint 4
  // 	m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
  // 				  KDL::Frame::DH_Craig1989(0.0, 1.5707963, 0.39, 0.0),
  // 				  KDL::Frame::DH_Craig1989(0.0, 1.5707963, 0.39, 0.0).Inverse()*KDL::RigidBodyInertia(2,
  // 										   KDL::Vector(0.0,0.5216809,0.0),
  // 										   KDL::RotationalInertia(-1.0436952,0.0,-1.0392780,0.0,0.0,0.0005324))));
  //
  // 	//joint 5
  // 	m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
  // 				  KDL::Frame::DH_Craig1989(0.0, 1.5707963, 0.0, 0.0),
  // 				  KDL::Frame::DH_Craig1989(0.0, 1.5707963, 0.0, 0.0).Inverse()*KDL::RigidBodyInertia(2,
  // 										   KDL::Vector(0.0,0.0119891,0.0),
  // 										   KDL::RotationalInertia(0.0036654,0.0,0.0060429,0.0,0.0,0.0004226))));
  //
  // 	//joint 6
  // 	m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
  // 				  KDL::Frame::DH_Craig1989(0.0, -1.5707963, 0.0, 0.0),
  // 				  KDL::Frame::DH_Craig1989(0.0, -1.5707963, 0.0, 0.0).Inverse()*KDL::RigidBodyInertia(2,
  // 										   KDL::Vector(0.0,0.0080787,0.0),
  // 										   KDL::RotationalInertia(0.0010431,0.0,0.0036376,0.0,0.0,0.0000101))));
  //
  //     //Frame 8 - end-effector (link 8) frame - at same pose as joint 7, frame == indentity!
  // 	//joint 7
  // 	m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
  // 				   KDL::Frame::Identity(),
  // 				   KDL::RigidBodyInertia(2,
  //     									   KDL::Vector::Zero(),
  //     									   KDL::RotationalInertia(0.000001,0.0,0.0001203,0.0,0.0,0.0))));
      //In total 8 joints (counting fixed - 0), 8 segments (counting base link 0) and 9 frames
  //In total 7 joints (NOT counting fixed - 0), 7 segments (NOT counting base link 0) and 8 frames
 //  m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
 //           KDL::Frame(KDL::Rotation::RotY(M_PI_2),KDL::Vector(0.0,0.0,1.0)), //local assignments of frames
 //           KDL::RigidBodyInertia(1.0, KDL::Vector(2.0, 0.0, 0.0),
 //           KDL::RotationalInertia(1.0,1.0,1.0)))); //centre of gravity
 // //
 //  m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
 //           KDL::Frame(KDL::Rotation::Identity(),KDL::Vector(-1.0,0.0,0.0)),
 //           KDL::RigidBodyInertia(1.0, KDL::Vector(2.0,2.0,0.0), //centre of gravity
 //           KDL::RotationalInertia(1.0,1.0,1.0))));
 //  //
 //  m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotY),
 //           KDL::Frame(KDL::Rotation::Identity(),KDL::Vector(0.0,0.0,1.0)), //local assignments of frames
 //           KDL::RigidBodyInertia(1.0, KDL::Vector(2.0, 0.0, 0.0),
 //           KDL::RotationalInertia(1.0,1.0,1.0)))); //centre of gravity
 //
 //  m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
 //           KDL::Frame(KDL::Rotation::Identity(),KDL::Vector(0.0,0.0,1.0)), //local assignments of frames
 //           KDL::RigidBodyInertia(1.0, KDL::Vector(2.0, 0.0, 0.0),
 //           KDL::RotationalInertia(1.0,1.0,1.0)))); //centre of gravity

//PUMA 560
     // m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
     //                            KDL::Frame::DH(0.0,M_PI_2,0.0,0.0),
     //                            KDL::RigidBodyInertia(0,KDL::Vector::Zero(),KDL::RotationalInertia(0,0.35,0,0,0,0))));
     // m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
     //                            KDL::Frame::DH(0.4318,0.0,0.0,0.0),
     //                            KDL::RigidBodyInertia(17.4,KDL::Vector(-.3638,.006,.2275),KDL::RotationalInertia(0.13,0.524,0.539,0,0,0))));
     // m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
     //                            KDL::Frame::DH(0.0,-M_PI_2,0.0,0.0),
     //                            KDL::RigidBodyInertia(4.8,KDL::Vector(-.0203,-.0141,.070),KDL::RotationalInertia(0.066,0.086,0.0125,0,0,0))));
     // m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
     //                            KDL::Frame::DH(0.0,M_PI_2,0.4318,0.0),
     //                            KDL::RigidBodyInertia(0.82,KDL::Vector(0,.019,0),KDL::RotationalInertia(1.8e-3,1.3e-3,1.8e-3,0,0,0))));
     // m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
     //                            KDL::Frame::DH(0.0,-M_PI_2,0.0,0.0),
     //                            KDL::RigidBodyInertia(0.34,KDL::Vector::Zero(),KDL::RotationalInertia(.3e-3,.4e-3,.3e-3,0,0,0))));
     // m.my_bot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
     //                            KDL::Frame::DH(0.0,0.0,0.0,0.0),
     //                            KDL::RigidBodyInertia(0.09,KDL::Vector(0,0,.032),KDL::RotationalInertia(.15e-3,0.15e-3,.04e-3,0,0,0))));


  // m.solve_fk(m.my_bot);
  // m.solve_jac(m.my_bot);
  create_motion_specification(my_motion);
  //traditional singularity
  KDL::ChainJntToJacSolver jacsolver(m.my_bot);
  jacsolver.JntToJac(my_motion.q,my_motion.jac);
  std::cout << "Traditional singularity "<< my_motion.jac << std::endl;
  Eigen::BDCSVD<Eigen::Matrix<double,6,Eigen::Dynamic> > svd(my_motion.jac.data);
  std::cout << "Rank of Jacobain " << svd.rank() << std::endl;

  //recursive newton Euler
  // KDL::ChainIdSolver_RNE rne(m.my_bot,KDL::Vector(0.0,0.0,-9.0)); //vector of gravity and all depends on attaching the robot to the base, directions - right hand rule
  // rne.CartToJnt(my_motion.q,my_motion.q_dot,my_motion.q_dot_dot, my_motion.external_force, my_motion.feedforward_torque); //Assumption is centre of mass and centre of gravity is the same- caanot be applied to satellities
  // std::cout<< "RNE_Solver : "<< my_motion.feedforward_torque << std::endl;

//vereshchagin solver
  my_KDL::ChainIdSolver_Vereshchagin vereshchagin(m.my_bot,
                                 KDL::Twist(KDL::Vector(0.0,0.0,-9.8),KDL::Vector(0.0,0.0,0.0)),
                                 6);
 //vereshchagin.CartToJnt(q, q_dot, q_dot_dot, alpha, beta, f_ext, torques);
 int result = vereshchagin.CartToJnt(my_motion.q,my_motion.q_dot,my_motion.q_dot_dot, my_motion.alpha, my_motion.beta, my_motion.external_force, my_motion.feedforward_torque);
 // std::cout<< "Result :"  << result << std::endl;
 std::cout<< " Vereshchagin_solver q_dot_dot : " << my_motion.q_dot_dot << " Torques: " << my_motion.feedforward_torque << std::endl;
 return 0;
}
  //segment1
//   KDL::Chain my_robot;
//   KDL::Joint joint1("joint1", KDL::Joint::RotZ);
//   KDL::Frame tip1(
//   		KDL::Rotation::RotY(M_PI_2),//local assignments of frames
//   		KDL::Vector(0.0,0.0,1.0)
//   	);
//   KDL::RigidBodyInertia inertia_of_body1(1.0,
//   			KDL::Vector(2.0, 0.0, 0.0),
//   			KDL::RotationalInertia(1.0,1.0,1.0)
//   	);
//
//   KDL::Segment segment1("link1",joint1, tip1, inertia_of_body1);
//   my_robot.addSegment(segment1);
//
// //segment2
//   // my_robot.addSegment(
//   // 	KDL::Segment("link2",
//   // 		KDL::Joint("joint2", KDL::Joint::RotZ),
//   // 		KDL::Frame (
// 	//   		KDL::Rotation::Identity(),
// 	//   		KDL::Vector(-1.0,0.0,0.0)),
// 	//   	KDL::RigidBodyInertia(
// 	//   		1.0,
//   // 			KDL::Vector(2.0,2.0,0.0), //centre of gravity
//   // 			KDL::RotationalInertia(1.0,1.0,1.0)) //Inertia of the body
//   // 		));
//   // std::cout<<" Joints : " << my_robot.getNrOfJoints() << std::endl;
//   // std::cout<<" Links : "  << my_robot.getNrOfSegments() << std::endl;
//   my_robot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ),
//            KDL::Frame(KDL::Rotation::Identity(),KDL::Vector(-1.0,0.0,0.0)),
//            KDL::RigidBodyInertia(1.0, KDL::Vector(2.0,2.0,0.0), //centre of gravity
//            KDL::RotationalInertia(1.0,1.0,1.0))));//Inertia of the body
//   // std::cout<<" Joints : " << my_robot.getNrOfJoints() << std::endl;
//   // std::cout<<" Links : "  << my_robot.getNrOfSegments() << std::endl;
//
//   myrobot.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ, 1, 0, c.joint_inertia[4]),
//           KDL::Frame::DH_Craig1989(0.0, 1.5707963, 0.0, 0.0),
//           KDL::Frame::DH_Craig1989(0.0, 1.5707963, 0.0, 0.0).Inverse()*KDL::RigidBodyInertia(2,
//                        KDL::Vector(0.0,0.0119891,0.0),
//                        KDL::RotationalInertia(0.0036654,0.0,0.0060429,0.0,0.0,0.0004226))));
  //3rd segment
  // my_robot.addSegment(
  // 	KDL::Segment("link2",
  // 		KDL::Joint("joint2", KDL::Joint::RotZ),
  // 		KDL::Frame (
	//   		KDL::Rotation::Identity(),
	//   		KDL::Vector(-1.0,0.0,0.0)),
	//   	KDL::RigidBodyInertia(
	//   		1.0,
  // 			KDL::Vector(2.0,2.0,0.0), //centre of gravity
  // 			KDL::RotationalInertia(1.0,1.0,1.0)) //Inertia of the body
  // 		));
  // std::cout<<" Joints : " << my_robot.getNrOfJoints() << std::endl;
  // std::cout<<" Links : "  << my_robot.getNrOfSegments() << std::endl;

  //  KDL::ChainFkSolverPos_recursive fksolver(m.my_bot);
  //  KDL::JntArray q(m.my_bot.getNrOfJoints());
  //   q(0) = 0;
  //   q(1) = 0;
  //   KDL::Frame end_effector_pose;
  //   fksolver.JntToCart(q, end_effector_pose);
  //   std::cout<<" FKSolver: "<< end_effector_pose << std::endl;
  //
//   KDL::JntArray q_dot(my_robot.getNrOfJoints());
//   q_dot(0) = 0.0; //(radians) //effect of centrifugal force (eg :spin a ball with thread)- velocity int the first joint produces a torque in the second joint
//   q_dot(1) = 0.0;
//
//   KDL::JntArray q_dot_dot(my_robot.getNrOfJoints());
//   q_dot_dot(0) = 0.0;
//   q_dot_dot(1) = 0.0;
//
//
//   KDL::Wrenches f_ext(my_robot.getNrOfSegments());
//   f_ext[0] = KDL::Wrench(
//     KDL::Vector(0.0,0.0,0.0),
//     KDL::Vector(0.0,0.0,0.0)
//   );
//   KDL::JntArray torques(my_robot.getNrOfJoints());
//
//   KDL::ChainIdSolver_RNE rne(my_robot,KDL::Vector(0.0,0.0,-9.0)); //vector of gravity and all depends on attaching the robot to the base, directions - right hand rule
//   rne.CartToJnt(q,q_dot,q_dot_dot, f_ext, torques); //Assumption is centre of mass and centre of gravity is the same- caanot be applied to satellities
//   std::cout<< "RNE_Solver : "<< torques << std::endl;
//
//
//   torques(0) = 0.0;
//   torques(1) = 0.0;
//   KDL::JntArray beta(6);
//   KDL::Jacobian alpha(6);
//   alpha.setColumn(0,
//     KDL::Twist(KDL::Vector(1.0,0.0,0.0),
//     KDL::Vector(0.0,0.0,0.0))
//   );
//   alpha.setColumn(1,
//     KDL::Twist(KDL::Vector(0.0,1.0,0.0),
//     KDL::Vector(0.0,0.0,0.0))
//   );
//   alpha.setColumn(2,
//     KDL::Twist(KDL::Vector(0.0,0.0,1.0),
//     KDL::Vector(0.0,0.0,0.0))
//   );
//   alpha.setColumn(3,
//     KDL::Twist(KDL::Vector(0.0,0.0,0.0),
//     KDL::Vector(1.0,0.0,0.0))
//   );
//   alpha.setColumn(4,
//     KDL::Twist(KDL::Vector(0.0,0.0,0.0),
//     KDL::Vector(0.0,1.0,0.0))
//   );
//   alpha.setColumn(5,
//     KDL::Twist(KDL::Vector(0.0,0.0,0.0),
//     KDL::Vector(0.0,0.0,1.0))
//   );
//   beta(0) = 2.0;
//   beta(1) = 2.0;
//   beta(3) = 2.0;
//   beta(2) = 2.0;
//   beta(4) = 2.0;
//   beta(5) = 2.0;
//
//   my_KDL::ChainIdSolver_Vereshchagin vereshchagin(my_robot,
//                                   KDL::Twist(KDL::Vector(0.0,0.0,0.0),KDL::Vector(0.0,0.0,0.0)),
//                                   6);
//   //vereshchagin.CartToJnt(q, q_dot, q_dot_dot, alpha, beta, f_ext, torques);
//   int result = vereshchagin.CartToJnt(q, q_dot, q_dot_dot, alpha, beta, f_ext, torques);
//   std::vector<KDL::Twist> xDotDot;
//   xDotDot.resize((my_robot.getNrOfSegments()+1));
//   vereshchagin.get_transformed_link_acceleration(xDotDot);
//   // std::cout << "Frame ACC" << '\n';
//     for (size_t i = 0; i < my_robot.getNrOfSegments()+1; i++) {
//         std::cout << xDotDot[i] << '\n';
// }
//   std::cout<<"beta"<<beta<<std::endl;
//  // std::cout<< "Alpha - Unit_constraint_forces matrix" << alpha << '\n';
//  // std::cout<< "Rows in unit constrained matrix" << alpha.rows()<< std::endl;
//  // for(unsigned int i=0; i<alpha.rows() ;i++){
//  //   for(unsigned int j=0; j<alpha.columns();j++){
//  //     std::cout<<alpha[i][j]<< std::endl;
//  //   }
//  // }
//  // std::cout<< "Frame accelaration"<<xDotDot[2]<< '\n';
//  // std::cout<<"Beta:"<< beta << std::endl;
// //       for(int i=0; i<6; i++){
// //            if(xDotDot[i]==beta[j]){
// //              std::cout<< "same values";
// //          }
// //       }
// //    }
//   std::cout<< "Result : " << result << std::endl;
//   std::cout<< " Vereshchagin_solver q_dot_dot : " << q_dot_dot << " Torques: " << torques << std::endl;
//   return 0;



//try extending kinematic chain, with 6DOF and constraints -also depends on configurations - provide joint level torques and ext forces - what do they controlhow
