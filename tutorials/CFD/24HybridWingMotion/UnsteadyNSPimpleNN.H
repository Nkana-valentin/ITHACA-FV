#ifndef UnsteadyNSPimpleNN_H
#define UnsteadyNSPimpleNN_H

#include <torch/script.h>
#include "torch2Eigen.H"
//#include "forces.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "fsiBasic.H"
#include "IOmanip.H"
#include <Eigen/SparseLU>
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include <pointPatchField.H>
#include <points0MotionSolver.H>

using namespace ITHACAtorch::torch2Eigen;
using namespace torch::indexing;

class UnsteadyNSPimpleNN : public fsiBasic
{

    public:
	    /// Constructors
	    /// Construct Null
		UnsteadyNSPimpleNN();
		/// Construct with argc and argv
		UnsteadyNSPimpleNN(int argc, char* argv[]);
		~UnsteadyNSPimpleNN();

		/// Members of the class
		Eigen::MatrixXd bias_inp;
		Eigen::MatrixXd scale_inp;
		Eigen::MatrixXd bias_out;
		Eigen::MatrixXd scale_out;
		Eigen::MatrixXd coeffL2U, coeffL2Nut;
		//Eigen::MatrixXd zmomentforces, yliftforces;
		torch::Tensor coeffL2U_tensor;
		torch::Tensor coeffL2Nut_tensor;
		/// List of POD coefficients
		List<Eigen::MatrixXd> CoeffU, CoeffP,  CoeffsNut; 
		torch::jit::script::Module netTorchscript, lstmTorchscript;
		 /// List of snapshots for the solution for eddy viscosity
		PtrList<volScalarField> nutFields;
		// PtrList<surfaceScalarField> phiFields;
		/// List of POD modes for eddy viscosity
		volScalarModes nutModes;
		/// Eddy viscosity field
		autoPtr<volScalarField> _nut;
		List<scalar> Rx, Ry, Rz, Rt;


		void loadLstmNet(word filename);
		void loadNet(word filename);
		void getTurbNN();

		Eigen::MatrixXd evalNet(Eigen::MatrixXd a);
		Eigen::MatrixXd evalLstm(Eigen::MatrixXd a);


		void truthSolve(List<scalar> mu_now, fileName folder);

};

#endif