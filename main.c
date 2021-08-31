#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <nlopt.h>

#define pi 3.14159265
///* Body numbers from _i file */
#define GROUND	(-1)
#define BODY1	0
#define BODY2	1
#define BODY3	2

/* State variables.  The state variable corresponding to
 * a particular joint axis can be obtained from the
 * tp_pend_i file or by using the sdindx() routine.
 */
#define NQ 3
#define NU 3
#define NSTATE (NQ+NU)

//Define Parameters for Manipulator
const double l1 = 1.0;
const double l2 = 1.0;
const double l3 = 0.5;
const double x_ref = 1.0;
const double y_ref = 1.0;

//Initialize position in x and y
double x, y;
double cost;

//Initialize angles
double theta1, theta2, theta3;

//Assign x, y, & cost along with the respective optimization variables
void opt_variables(double opt_var[3]){

  //Joint angles
  opt_var[0] = theta1;
  opt_var[1] = theta2;
  opt_var[2] = theta3;
      
  //Global end effector position
  x = l1*cos(opt_var[0]) + l2*cos(opt_var[0] + opt_var[1]) + l3*cos(opt_var[0] + opt_var[1] + opt_var[2]);
  y = l1*sin(opt_var[0]) + l2*sin(opt_var[0] + opt_var[1]) + l3*sin(opt_var[0] + opt_var[1] + opt_var[2]);
      
  //The addition of all the squared errors of the constraints: x and y end effector position along with the end effector -90 angle
  cost = (x-x_ref)*(x-x_ref) + (y-y_ref)*(y-y_ref) + (opt_var[2] - (-pi/2.0 - opt_var[0] -opt_var[1]))*(opt_var[2] - (-pi/2.0 - opt_var[0] -opt_var[1]));
}

// Ojective Function

int count; //Keeps count of the number of iterations for convergence

double objective(unsigned n, double *opt_var, double *grad, void *objective_data)
{
  count++;
  
  //Gradient with respect to joint angles, but objective function is set to zero
  if (grad){
      grad[0] = 0.0;
      grad[1] = 0.0;
      grad[2] = 0.0;
  }
  return 0.0;
}

// Constraints

//Error for x end effector
double myconstraint0(unsigned n, const double *opt_var, double *grad, void *c0_data)
{
  if (grad) {
      grad[0] = -l1*sin(opt_var[0]) - l2*sin(opt_var[0] + opt_var[1]) - l3*sin(opt_var[0] + opt_var[1] + opt_var[2]);
      grad[1] = - l2*sin(opt_var[0] + opt_var[1]) - l3*sin(opt_var[0] + opt_var[1] + opt_var[2]);
      grad[2] = - l3*sin(opt_var[0] + opt_var[1] + opt_var[2]);
  }
  return (l1*cos(opt_var[0]) + l2*cos(opt_var[0] + opt_var[1]) + l3*cos(opt_var[0] + opt_var[1] + opt_var[2]) - x_ref);
  }

//Error for y end effector
double myconstraint1(unsigned n, const double *opt_var, double *grad, void *c1_data)
{
  if (grad) {
      grad[0] = l1*cos(opt_var[0]) + l2*cos(opt_var[0] + opt_var[1]) + l3*cos(opt_var[0] + opt_var[1] + opt_var[2]);
      grad[1] = l2*cos(opt_var[0] + opt_var[1]) + l3*cos(opt_var[0] + opt_var[1] + opt_var[2]);
      grad[2] = l3*cos(opt_var[0] + opt_var[1] + opt_var[2]);
  }
  return (l1*sin(opt_var[0]) + l2*sin(opt_var[0] + opt_var[1]) + l3*sin(opt_var[0] + opt_var[1] + opt_var[2]) - y_ref);
}

//-90 degree angle constraint for joint 3
double myconstraint2(unsigned n, const double *opt_var, double *grad, void *c2_data)
{
  if (grad) {
      grad[0] = 1.0;
      grad[1] = 1.0;
      grad[2] = 1.0;
  }
  return (opt_var[2] - (-opt_var[0] - opt_var[1] - pi/2.0));
}

int main()
{
  //Create a data file
  FILE * fp;
  fp = fopen("data.txt", "w+");

  double t=0,opt_var[NQ]={0},u[NU]={0};

  /* initialize sd/fast model */
   sdinit(); sdprinterr(stderr);

  /*set initial position opt_var and velocities u as needed*/
  opt_var[0] = 0;
  opt_var[1] = 0;
  opt_var[2] = -pi/2;
  fprintf(fp,  "%f\t %f\t %f\n", opt_var[0],opt_var[1],opt_var[2]);
  //fprintf(fp, "starting joint angles: %f\t %f\t %f\n",opt_var[0],opt_var[1],opt_var[2]);

  /*pass the state */
  sdstate(t,opt_var,u); sdprinterr(stderr);

  /*forward kinematics*/
  //location of end-effector on body3 wrt to center of mass in zero reference frame
  double local_end_effector[]={0.25,0,0}; //end-effector is at 0.25 wrt com of body3
  double global_end_effector[3]={0};
  sdpos(BODY3,local_end_effector,global_end_effector);

  fprintf(fp,"%f\t %f\n",global_end_effector[0],global_end_effector[1]); //x and y position
  //fprintf(fp, "end_effector start position: %f\t %f\n",global_end_effector[0],global_end_effector[1]); //x and y position

  // Setup and Run nlopt solver
  nlopt_opt opt; //Initialize opt
  opt = nlopt_create(NLOPT_LN_COBYLA, 3); /* algorithm and dimensionality where LD_MMA is for gradient based and LN_COBYLA is for non-gradient based; the # represents the number of optimization variables  */
  nlopt_set_min_objective(opt, objective, NULL); //Set to minimize objective function with previously specified algorithm

  nlopt_add_equality_constraint(opt, myconstraint0, NULL, 1e-8); //Initialize constraint 0 with zero inputs (NULL) and tolerance
  nlopt_add_equality_constraint(opt, myconstraint1, NULL, 1e-8); //Initialize constraint 1 with zero inputs (NULL) and tolerance
  nlopt_add_equality_constraint(opt, myconstraint2, NULL, 1e-8); //Initialize constraint 2 with zero inputs (NULL) and tolerance

  nlopt_set_xtol_rel(opt, 1e-4); //Set the tolerance for convergence criteria

  /*double opt_var[3] = { -pi/2, -pi/2, -pi/2 }; `*`some` `initial` `guess`*` */
  double minf; /* `*`the` `minimum` `objective` `value,` `upon` `return`*` */

  if (nlopt_optimize(opt, opt_var, &minf) < 0) {
        printf("nlopt failed!\n");
  }
  else {
      printf("found minimum after %d evaluations\n", count);
      /*printf("found minimum %0.10g\n for joint angles: theta1 = %g, theta2 = %g, theta3 = %g \n", minf, opt_var[0], opt_var[1], opt_var[2]);*/
      }
  nlopt_destroy(opt); //Necessary to close off the program

  //fprintf(fp, "final joint angles: %f\t %f\t %f\n",opt_var[0],opt_var[1],opt_var[2]);
  fprintf(fp,"%f\t %f\t %f\n",opt_var[0],opt_var[1],opt_var[2]);
    
  /*pass the state */
  sdstate(t,opt_var,u); sdprinterr(stderr);

  /*forward kinematics*/
  sdpos(BODY3,local_end_effector,global_end_effector);
    
  fprintf(fp,"%f\t %f\n",global_end_effector[0],global_end_effector[1]); //x and y position
  //fprintf(fp, "end_effector final position: %f\t %f\n",global_end_effector[0],global_end_effector[1]); //x and y position

  fclose(fp);

  return 0;
}


void sduforce(double t, double *q, double *u)
{

}
