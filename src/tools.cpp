#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse <<0,0,0,0;	
	float sum;
	
	int i;

	if (estimations.size()!=ground_truth.size() || estimations.size()==0){
		cout << "Estimation or Ground Data is invalid, please check your data file."<< endl;
		return rmse;
	}

	sum = 0.0;
	for (i=0;i<estimations.size();++i){
		VectorXd diff = estimations[i] -  ground_truth[i];
		diff = diff.array() * diff.array();
		rmse += diff;
	}

	rmse = rmse/estimations.size();
	rmse = rmse.array().sqrt();

	//cout << " rmse = " << rmse << endl;
	return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

	MatrixXd Hj(3,4);
	//state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);
	
	float c1 = px*px + py*py;
	float c2 = sqrt(c1);
	float c3 = sqrt(c1*c1*c1);
	float tmp1 = vx*py-vy*px;
	float tmp2 = vy*px-vx*py;

	//cout << " c1 = " << c1 << endl;
	//cout << " c2 = " << c2 << endl;
	//cout << " c3 = " << c3 << endl;

	if (fabs(c1) < 0.0001 || fabs(c2) < 0.0001 || fabs(c3) < 0.0001){
		cout << "Error: divided by 0.0:" << endl;
		return Hj;
	}

	// Jacobian matrix
	
	Hj << px/c2, py/c2, 0, 0,
	      -py/c1,px/c1, 0, 0,
	      (py*tmp1)/c3,(px*tmp2)/c3,px/c2,py/c2;

	return Hj;		
}
