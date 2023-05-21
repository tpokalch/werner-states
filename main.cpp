
#include </usr/local/Cellar/gsl/2.7.1/include/gsl/gsl_multimin.h>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <Eigen/Dense>
#include <math.h>

using namespace std;
const std::complex<double> im(0.0f, 1.0f); 

//typedef Eigen::Matrix<std::complex<double>, 2, 2> ComplexMatrix2d;
//typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> Eigen::MatrixXcd;

//	MyMatrix is COMPLEX despite it saying Xd (must be Xcd)



class MyMatrixXd : public Eigen::MatrixXcd
{
	public:
		MyMatrixXd(int row, int col) : Eigen::MatrixXcd(row, col)
		{
//			std::cout << "i hope the base class initializes the members" << std::endl;
		}
		MyMatrixXd()
		{
//			std::cout<<"\n Constructor executed";
		}
		void	Conjugate()
		{
			this->transposeInPlace();

			int rowb = this->rows();
			int colb = this->cols();

			printf("rows cols %d,%d\n", rowb, colb);

			for (int i = 0; i < colb; i++)
			{
				for (int j = 0; j < rowb; j++)
				{
					(*this)(i,j) = conj((*this)(i,j));
				}
			}			
	
		}
		void	Kron(MyMatrixXd &M2, MyMatrixXd &result)
		{
			MyMatrixXd *M1 = this;
			int rowa = M1->rows();
			int rowb = M2.rows();

			int cola = M1->cols();
			int colb = M2.cols();

			 for (int i = 0; i < rowa; i++) {
				for (int k = 0; k < rowb; k++) {
					for (int j = 0; j < cola; j++) {
						for (int l = 0; l < colb; l++) {
 
                    // Each element of matrix A is
                    // multiplied by whole Matrix B
                    // resp and stored as Matrix C
                    result(k + i * rowb, l + j * colb) = (*M1)(i, j) * M2(k, l);
                   // std::cout << result(i + l, j + k) << " ";
								}
							}
//							std::cout << std::endl;
						}
					}
		//	return result;
		}
		MyMatrixXd operator= (const Eigen::MatrixXcd& M2)
		{
			int rowb = M2.rows();
			int colb = M2.cols();
			MyMatrixXd result(rowb,colb);
//			std::cout << "assigned matrix is" << std::endl << M2;
//			printf("this doesn'w work, don't use!\n");
			for (int i = 0; i < colb; i++)
			{

				for (int j = 0; j < rowb; j++)
				{
					std::cout << "assign " << M2(j, i) << "to " << i << j << std::endl;

					result(j,i) = M2(j,i);
				}
//				printf("return result\n");
			}
			return (result);	
		}


/*		MyMatrixXd operator= (const MyMatrixXd& M2)
		{
			int rowb = M2.rows();
			int colb = M2.cols();
//			MyMatrixXd result(rowb,colb);

			for (int i = 0; i < colb; i++)
			{
				for (int j = 0; j < rowb; j++)
				{
					printf("rows of M2 %d\n", M2.rows());
					printf("rows of this %d\n", this->rows());
					printf("cols of M2 %d\n", M2.cols());
					printf("cols of this %d\n", this->cols());
					printf("assigning %d,%d\n", i, j);

					(*this)(j,i) = M2(j,i);
				}
				printf("return result\n");
			}
			return (*this);	
		}
*/

		MyMatrixXd operator+ (const MyMatrixXd& M2)
		{
			int rowb = M2.rows();
			int colb = M2.cols();
			MyMatrixXd result(rowb,colb);

			for (int i = 0; i < colb; i++)
			{
				for (int j = 0; j < rowb; j++)
				{
/*					printf("rows of M2 %d\n", M2.rows());
					printf("rows of this %d\n", this->rows());
					printf("cols of M2 %d\n", M2.cols());
					printf("cols of this %d\n", this->cols());
					printf("assigning %d,%d\n", i, j);
*/
					result(j,i) = (*this)(j, i) + M2(j,i);
				}
//				printf("return result\n");
			}
			return (result);	
		}


		virtual ~MyMatrixXd()
		{
//			std::cout<<"\n Destructor executed";
		}
};


typedef struct s_param
{
	double t;
	int d;

	MyMatrixXd Usw;
//	MyMatrixXd C;
//	MyMatrixXd S;
	Eigen::MatrixXcd W[4][4];


	double *p;
	int d2;

	int degrees_in_vector;

	int number_of_terms;
	int number_of_weights;

	int degrees_of_freedom;
	int line_length;
	int lines;

	MyMatrixXd Mx;
	MyMatrixXd My;
	MyMatrixXd MxMy;
	Eigen::MatrixXcd M;

	MyMatrixXd x;
	MyMatrixXd y;
	MyMatrixXd xt;
	MyMatrixXd yt;

	Eigen::MatrixXcd Zero;
	MyMatrixXd Identity;

	gsl_vector *v;
}	t_param;

//inits i'th vector;
void		init_with_known_solution(t_param *param, const gsl_vector *x, int i)
{
	double t = param->t;

	int j = i / param->d;
	int k = i % param->d;

	if (param->d == 2)
	{
	double a = 0.5 * (sqrt((1 + 3 * t) * 0.5) + sqrt(0.5 * (1 - t)));
	double b = 0.5 * (sqrt((1 + 3 * t) * 0.5) - sqrt(0.5 * (1 - t)));
	double r = 0.5 * (sqrt((1 + a) * (1 + a) - b * b) + sqrt((1 - a) * (1 - a) - b * b));
	double s = 0.5 * (sqrt((1 + a) * (1 + a) - b * b) - sqrt((1 - a) * (1 - a) - b * b));
	param->x(0,0) = r;
	param->x(1,0) = sqrt(1 - r * r) * polar(1.0, M_PI / 4.0);

	param->y(0,0) = s * polar(1.0, 0.0);
	param->y(1,0) = sqrt(1 - s * s) * polar(1.0, M_PI / 4.0);

//	now account for the fact that i'th vector is requested

//	printf("j, k is %d,%d\n", j, k);
//	printf("here\n");
	
//	std::cout << "x was " << std::endl << param->x << std::endl;

//	assigning form Eigen::Matrix to MyMatrix doesn'twork!!! don't do!
/*
	param->x = param->W[j][k] * param->x;
	param->y = param->W[j][k] * param->y;

*/


	//stupid mistake!

/*
	complex<double>x0 = (param->W[j][k] * param->x)(0,0);
	complex<double>x1 = (param->W[j][k] * param->x)(1,0);
	param->x(0,0) = x0;
	param->x(1,0) = x1;

	complex<double>y0 = (param->W[j][k] * param->y)(0,0);
	complex<double>y1 = (param->W[j][k] * param->y)(1,0);
	param->y(1,0) = y1;
	param->y(0,0) = y0;
*/



//	std::cout << "x is " << std::endl << param->x << std::endl;


	}

	else if (param->d == 3)
	{
		double t = param->t;


		if (t == 1.0/(param->d + 1))
		{		
			param->x(0,0) = sqrt(2/3.0);
			param->x(1, 0) = -1 / sqrt(6);
			param->x(2, 0) = -1/sqrt(6);

			param->y(0,0) = sqrt(2/3.0);
			param->y(1, 0) = -1 / sqrt(6);
			param->y(2, 0) = -1/sqrt(6);
		
//		printf("1: %f,%f,%f\n%f,%f,%f\n", std::abs(param->x(0, 0)), std::abs(param->x(1, 0)), std::abs(param->x(2, 0)), std::abs(param->y(0, 0)), std::abs(param->y(1, 0)), std::abs(param->y(2, 0)));
		}
/*
		complex<double>x0 = (param->W[j][k] * param->x)(0,0);
		complex<double>x1 = (param->W[j][k] * param->x)(1,0);
		complex<double>x2 = (param->W[j][k] * param->x)(2,0);


		param->x(0,0) = x0;
		param->x(1,0) = x1;
		param->x(2,0) = x2;

		complex<double>y0 = (param->W[j][k] * param->y)(0,0);
		complex<double>y1 = (param->W[j][k] * param->y)(1,0);
		complex<double>y2 = (param->W[j][k] * param->y)(2,0);


		param->y(1,0) = y1;
		param->y(0,0) = y0;
		param->y(2,0) = y2;

*/
		else{

		double a1 = 1 + 7 * t - 8 * t * t;
		double a = sqrt((5 + 4 * t + 4 * sqrt(a1)) / 81.0);
		double p = (1 + 2 * t - sqrt(a1)) / (1 - 4 * t);
		double l = (-1 + sqrt( 1 + 4 * p)) / 2.0;
		double B = -a * (l + 1);

//		printf("t, a, a1, p, l B: %f, %f, %f, %f, %f, %f\n", t, a, a1, p, l, B);

		Eigen::MatrixXd u(3,1);
		Eigen::MatrixXd v(3,1);

		u(0, 0) = 1;
		u(1, 0) = l;
		u(2, 0) = l;

		v(0, 0) = a;
		v(1, 0) = B;
		v(2, 0) = B;





		param->x(0, 0) = (sqrt(3)  * sqrt((v.adjoint() * v)(0, 0)) * u)(0, 0);
		param->x(1, 0) = (sqrt(3)  * sqrt((v.adjoint() * v)(0, 0)) * u)(1, 0);
		param->x(2, 0) = (sqrt(3)  * sqrt((v.adjoint() * v)(0, 0)) * u)(2, 0);


//		printf("|v| = %f\n", (v.adjoint() * v )(0, 0));
//		printf("|v| = %f\n", sqrt((v.adjoint() * v )(0, 0)));


		param->y(0, 0)= (v / sqrt((v.adjoint() * v )(0, 0)))(0, 0);
		param->y(1, 0)= (v / sqrt((v.adjoint() * v )(0, 0)))  (1, 0);
		param->y(2, 0)= (v / sqrt((v.adjoint() * v )(0, 0)))  (2, 0);


//		printf("2: %f,%f,%f\n%f,%f,%f\n", std::abs(param->x(0, 0)), std::abs(param->x(1, 0)), std::abs(param->x(2, 0)), std::abs(param->y(0, 0)), std::abs(param->y(1, 0)), std::abs(param->y(2, 0)));


		}
	}
	else if (param->d == 4)
	{
/*
		double r0 = (1 - 1.0 / sqrt(5)) / (2 * sqrt(2 - sqrt(2)));
		double r1 = (sqrt(2) - 1) * r0;
		double r_p = 0.5 * sqrt(1.0 + 1.0 / sqrt(5) + sqrt(1.0 / 5.0 + 1.0 / sqrt(5)));
//		double r_p = 0.5 * (r0 + sqrt(2 - 3 * r0 * r0));
//		double r_m = 0.5 * (r0 - sqrt(2 - 3 * r0 * r0));


		double r_m = 0.5 * sqrt(1.0 + 1.0 / sqrt(5) - sqrt(1.0 / 5.0 + 1.0 / sqrt(5)));

		printf("r0, r1, r_p, r_m %f,%f,%f,%f\n", r0, r1, r_p, r_m);

		printf("r0^2 + r1^2 + r_p^2 + r_m^2 = %f\n", r0 * r0 + r1 * r1 + r_p * r_p + r_m * r_m);


		double R = sqrt(r0 * r0 + r1 * r1 + r_p * r_p + r_m * r_m);
		printf("r0 / R,  %f\n", r0 / R);



		double a = acos(2 / sqrt(5 + sqrt(5)));
		double b = asin(2 / sqrt(5));
		double Th_p = a / 2 + b / 4 + M_PI / 4.0;
		double Th1 = M_PI / 2.0;
		double Th_m = -a / 2 + b / 4 + M_PI / 4.0; 

		param->x(0,0) = r0 / R ;
		param->x(1, 0) = r_p / R* polar(1.0, Th_p);
		param->x(2, 0) = r1 / R* polar(1.0, Th1);
		param->x(3, 0) = r_m / R* polar(1.0, Th_m);
*/

		param->x(0, 0) = 0.485712214091264;
		param->x(1, 0) = 0.600433696-0.4498963669 * im;
		param->x(2, 0) = -0.201188586 *im;
		param->x(3, 0) = -0.3992451-0.0358158471 * im;


		param->y(0,0) = param->x(0, 0);
		param->y(1, 0) = param->x(1, 0);
		param->y(2, 0) = param->x(2, 0);
		param->y(3, 0) = param->x(3, 0);

//////////////////////////////////
/*
		complex<double>x0 = (param->W[j][k] * param->x)(0,0);
		complex<double>x1 = (param->W[j][k] * param->x)(1,0);
		complex<double>x2 = (param->W[j][k] * param->x)(2,0);
		complex<double>x3 = (param->W[j][k] * param->x)(3,0);


		param->x(0,0) = x0;
		param->x(1,0) = x1;
		param->x(2,0) = x2;
		param->x(3,0) = x3;



		complex<double>y0 = (param->W[j][k] * param->y)(0,0);
		complex<double>y1 = (param->W[j][k] * param->y)(1,0);
		complex<double>y2 = (param->W[j][k] * param->y)(2,0);
		complex<double>y3 = (param->W[j][k] * param->y)(3,0);



		param->y(0,0) = y0;
		param->y(1,0) = y1;
		param->y(2,0) = y2;
		param->y(3,0) = y3;

*/


//		printf("not implemented!\n");
		if(t == 0)
		{
//			NOT FINISHED !!
			int d = param->d;
//			Eigen::MatrixXcd(d, 1) a[d];
//			Eigen::MatrixXcd(d, 1) b[d];
			Eigen::MatrixXcd a[d];
			Eigen::MatrixXcd b[d];


			for(int i = 0; i < d; i++)
			{
/*
				for(int j = 0; j < d; i++)
				{

					a[i](i, 0) = 1;
					b[i](i, 0) = 1;
				}*/
				a[i] = Eigen::MatrixXcd(d, 1);
				a[i] = param->Identity.col(i);
				std::cout << "a[i] is "<< a[i] << endl;

				b[i] = Eigen::MatrixXcd(d, 1);
				b[i] = param->Identity.col(i);
	
			}
			for(int k = 0; k < d; k++)
			{

				for(int l = 0; l < d; l++)
				{
					x[k][l] = a[k];
					y[k][l] = b[l];
				}
			}

			param->x(0, 0) = 0.485712214091264;
			param->x(1, 0) = 0.600433696-0.4498963669 * im;
			param->x(2, 0) = -0.201188586 *im;
			param->x(3, 0) = -0.3992451-0.0358158471 * im;


		}
	}

	complex<double>X[param->d];
	complex<double>Y[param->d];

	for (int l = 0; l < param->d; l++)
		X[l] = (param->W[j][k] * param->x)(l, 0);
	for (int l = 0; l < param->d; l++)
		param->x(l, 0) = X[l];
	for (int l = 0; l < param->d; l++)
		Y[l] = (param->W[j][k] * param->y)(l, 0);
	for (int l = 0; l < param->d; l++)
		param->y(l, 0) = Y[l];

	//param->x(0, 0) must be real!

//	printf("change phase\n");

	double phasex = std::arg(param->x(0, 0));
	for (int l = 0; l < param->d; l++)
		param->x(l, 0) = param->x(l, 0) * polar(1.0, -phasex);

//	printf("param->x(0, 0).imag() = %f\n", param->x(0, 0).imag()); 
	double phasey = std::arg(param->y(0, 0));
	for (int l = 0; l < param->d; l++)
		param->y(l, 0) = param->y(l, 0) * polar(1.0, -phasey);
//	printf("known solution :\n");

//	std::cout << "x["<<i<<"] " << "is " << std::endl << param->x << std::endl;
//	std::cout << "y["<<i<<"] " <<" is " << std::endl << param->y << std::endl;

	param->xt = param->x;
	param->yt = param->y;
	param->xt.adjointInPlace() ;
	param->yt.adjointInPlace();



}

 

//	inits i'th pair of vectors, writes to param->x and param-y
void		init_xy(t_param *param, const gsl_vector *x, int i)
{


//	if uncomment set minimizint iterations to 2

/*
		init_with_known_solution(param, x, i);
		return ;
*/

//#if 0
		if (param->d == 2)
		{
//		PARAMETRISATION FOR d = 2

		int DELETE_ONE_TEST;

		double B, A, Y, X, w;

		Y = gsl_vector_get(x, 0 + param->line_length * i);
		X = gsl_vector_get(x, 1 + param->line_length * i);
		B = gsl_vector_get(x, 2 + param->line_length * i);
		A = gsl_vector_get(x, 3 + param->line_length * i);


//		std::complex<double> x[2];
//		std::complex<double> y[2];

//		MyMatrixXd x(d,1);
//		MyMatrixXd y(d,1);

//		for future transposition
//		MyMatrixXd xt(d,1);
//		MyMatrixXd yt(d,1);


//		init |x> and |y> from parameters A, B, X, Y, alway explicitly
//		passed to my_f() so that <x|x> = <y|y> = 1

		param->x(0, 0) = std::abs(cos(Y)); //					cos(Y)
//					magnitude (>0 always)
//		x(1, 0) = sgn(sin(Y)) * polar(abs(sin(Y)), X);//	sin(Y) * e^iX
		param->x(1, 0) = sin(Y) * polar(1.0, X);//	sin(Y) * e^iX
		param->y(0, 0) = std::abs(cos(B)); //					cos(B)
//		y(1, 0) = sgn(sin(B)) * polar(abs(sin(B)), A);//	sin(B) * e^iA
		param->y(1, 0) = sin(B) * polar(1.0, A);//	sin(B) * e^iA

//		std::cout << "vector is " << param->x << std::endl;

////		init x[i > 0], y[i > 0] with weyl-heisenberg



	if (i > 0)
	{
#if 0
		//init fiducial vector
		double B1, A1, Y1, X1, w1;

		Y1 = gsl_vector_get(x, 0);
		X1 = gsl_vector_get(x, 1);
		B1 = gsl_vector_get(x, 2);
		A1 = gsl_vector_get(x, 3);


//		std::complex<double> x[2];
//		std::complex<double> y[2];

//		MyMatrixXd x(d,1);
//		MyMatrixXd y(d,1);

//		for future transposition
//		MyMatrixXd xt(d,1);
//		MyMatrixXd yt(d,1);


//		init |x> and |y> from parameters A, B, X, Y, alway explicitly
//		passed to my_f() so that <x|x> = <y|y> = 1

		param->x(0, 0) = std::abs(cos(Y1)); //					cos(Y)
//					magnitude (>0 always)
//		x(1, 0) = sgn(sin(Y)) * polar(abs(sin(Y)), X);//	sin(Y) * e^iX
		param->x(1, 0) = sin(Y1) * polar(1.0, X1);//	sin(Y) * e^iX
		param->y(0, 0) = std::abs(cos(B1)); //					cos(B)
//		y(1, 0) = sgn(sin(B)) * polar(abs(sin(B)), A);//	sin(B) * e^iA
		param->y(1, 0) = sin(B1) * polar(1.0, A1);//	sin(B) * e^iA


	int j = i / param->d;
	int k = i % param->d;

/*
	complex<double>x0 = (param->W[j][k] * param->x)(0,0);
	complex<double>x1 = (param->W[j][k] * param->x)(1,0);
	param->x(0,0) = x0;
	param->x(1,0) = x1;

	complex<double>y0 = (param->W[j][k] * param->y)(0,0);
	complex<double>y1 = (param->W[j][k] * param->y)(1,0);
	param->y(1,0) = y1;
	param->y(0,0) = y0;

*/
	complex<double>X[param->d];
	complex<double>Y[param->d];

	for (int l = 0; l < param->d; l++)
		X[l] = (param->W[j][k] * param->x)(l, 0);
	for (int l = 0; l < param->d; l++)
		param->x(l, 0) = X[l];
	for (int l = 0; l < param->d; l++)
		Y[l] = (param->W[j][k] * param->y)(l, 0);
	for (int l = 0; l < param->d; l++)
		param->y(l, 0) = Y[l];

double phasex = std::arg(param->x(0, 0));
	for (int l = 0; l < param->d; l++)
		param->x(l, 0) = param->x(l, 0) * polar(1.0, -phasex);

	double phasey = std::arg(param->y(0, 0));
	for (int l = 0; l < param->d; l++)
		param->y(l, 0) = param->y(l, 0) * polar(1.0, -phasey);

	//	std::cout << "x" << i << " : " << param->x << endl;

#endif

	}


		if (i == 0)
		{
//			|x0> = |0>
//			param->x(0, 0) = 1.0 + 0.0 * im;
//			param->x(1, 0) = 0.0;
//////////
/*
			param->x(0, 0) = 1.0 / sqrt(2);
			param->x(1, 0) = 1 / sqrt(2) * polar(1.0, 1.0);
*/
			// y дійсний
//			param->y(1, 0) = std::abs(param->y(1, 0));
		}

		else if (i == 1)
		{
			// |x1> дійсний
//			param->x(1, 0) = std::abs(param->x(1, 0));//misterious free param!

//			param->x(0, 0) = 0.0 + 0.0 * im;
//			param->x(1, 0) = 0.0 + 1.0 * im;

/*
			param->x(0, 0) = 0.22831118948877446306 + 0.00000000000000000000 * im;
			param->x(1, 0) =  0.71111332912290403652 + -0.66497506261360017632 * im;
*/


		}
		
		}
//#endif

//#if 0
		else if (param->d == 3)
		{

//		PARAMETRISATION FOR d = 3

//		printf("i is %d\n", i);

		double xal, xbe, xA, xB, xC;

		double yal, ybe, yA, yB, yC;


		int DELETE_ONE_TEST;

		xal = gsl_vector_get(x, 0 + param->line_length * i);
		xbe =  gsl_vector_get(x, 1 + param->line_length * i);
		xA =  gsl_vector_get(x, 2 + param->line_length * i);
		xB =  gsl_vector_get(x, 3 + param->line_length * i);

		yal = gsl_vector_get(x, 4 + param->line_length * i);
		ybe = gsl_vector_get(x, 5 + param->line_length * i);
		yA = gsl_vector_get(x, 6 + param->line_length * i);
		yB = gsl_vector_get(x, 7 + param->line_length * i);

		param->x(0, 0) = cos(xal) * cos(xbe);
		param->x(1, 0) = polar(1.0, xA) * cos(xal) * sin(xbe);
		param->x(2, 0) = polar(1.0, xB) * sin(xal);

		param->y(0, 0) = cos(yal) * cos(ybe);
		param->y(1, 0) = polar(1.0, yA) * cos(yal) * sin(ybe);
		param->y(2, 0) = polar(1.0, yB) * sin(yal);

		if (i == 0)
		{
/*			param->x(0, 0) = cos(xal) * cos(xbe);
			param->x(1, 0) = cos(xal) * sin(xbe);
			param->x(2, 0) = sin(xal);
*/
			param->x(0, 0) = param->x(0, 0).real();
			param->x(1, 0) = param->x(1, 0).real();
			param->x(2, 0) = param->x(2, 0).real();

		}

		}
//#endif 

//#if 0

		else if (param->d == 4)
		{
//		PARAMETRISATION FOR d = 4
		int DELETE_ONE_TEST;



		double xal, xbe, xga, xA, xB, xC;

		double yal, ybe, yga, yA, yB, yC;

//		поміняти починаючи з 1;

		xal = gsl_vector_get(x, 0 + param->line_length * i);
		xbe =  gsl_vector_get(x, 1 + param->line_length * i);
		xga =  gsl_vector_get(x, 2 + param->line_length * i);
		xA =  gsl_vector_get(x, 3 + param->line_length * i);
		xB =  gsl_vector_get(x, 4 + param->line_length * i);
		xC =  gsl_vector_get(x, 5 + param->line_length * i);

		yal = gsl_vector_get(x, 6 + param->line_length * i);
		ybe = gsl_vector_get(x, 7 + param->line_length * i);
		yga = gsl_vector_get(x, 8 + param->line_length * i);
		yA = gsl_vector_get(x, 9 + param->line_length * i);
		yB = gsl_vector_get(x, 10 + param->line_length * i);
		yC = gsl_vector_get(x, 11 + param->line_length * i);




		param->x(0, 0) = cos(xal) * cos(xbe) * cos(xga);
		param->x(1, 0) = polar(1.0, xA) * cos(xal) * sin(xbe) * cos(xga);
		param->x(2, 0) = polar(1.0, xB) * sin(xal) * cos(xga);
		param->x(3, 0) = polar(1.0, xC) * sin(xga);


		param->y(0, 0) = cos(yal) * cos(ybe) * cos(yga);
		param->y(1, 0) = polar(1.0, yA) * cos(yal) * sin(ybe) * cos(yga);
		param->y(2, 0) = polar(1.0, yB) * sin(yal) * cos(yga);
		param->y(3, 0) = polar(1.0, yC) * sin(yga);
		}

//# endif

		param->xt = param->x;
		param->yt = param->y;
		param->xt.adjointInPlace() ;
		param->yt.adjointInPlace();

//		std::cout << sqrt((t * (param->d2 - 1) + 1) / param->d) << std::endl;
	}

double*	get_angles_from_vector(t_param *p, Eigen::MatrixXcd &x, Eigen::MatrixXcd &y)
{
	double *ret = NULL;
	if (p->d == 2)
	{

		ret = (double *)malloc(sizeof(double) * 2 * 2);
/*		param->x(0, 0) = cos(Y); //					cos(Y)
//					magnitude (>0 always)
//		x(1, 0) = sgn(sin(Y)) * polar(abs(sin(Y)), X);//	sin(Y) * e^iX
		param->x(1, 0) = sin(Y) * polar(1.0, X);//	sin(Y) * e^iX
		param->y(0, 0) = cos(B); //					cos(B)
//		y(1, 0) = sgn(sin(B)) * polar(abs(sin(B)), A);//	sin(B) * e^iA
		param->y(1, 0) = sin(B) * polar(1.0, A);//	sin(B) * e^iA
*/
		double cosY = std::abs(x(0, 0));
		double Y = acos(cosY);
		complex<double> XsinY = x(1, 0);
		double X = std::arg(XsinY);
//		double sinY = std::abs(XsinY);

		double cosB = std::abs(y(0, 0));
		double B = acos(cosB);
		complex<double> AsinB = y(1, 0);
		double A = std::arg(AsinB);
//		double sinB = std::abs(AsinB);


		ret[0] = Y;
		ret[1] = X;
		ret[2] = B;
		ret[3] = A;

		printf("angles are %f,%f,%f,%f\n", Y, X, B, A);
/*
		Y = gsl_vector_get(x, 0 + param->line_length * i);
		X = gsl_vector_get(x, 1 + param->line_length * i);
		B = gsl_vector_get(x, 2 + param->line_length * i);
		A = gsl_vector_get(x, 3 + param->line_length * i);
*/

	}
	if(p->d == 3)
	{
		ret = (double *)malloc(sizeof(double) * 4 * 2);
//		printf("NOT IMPLEMENTED!!!\n");

/////////////////////////////////////////////////////
	//////////////////////////////////////////////

		double sinxal = std::abs(x(2, 0));
		double sinxal2 = std::norm(x(2, 0));
		double xB = std::arg(x(2, 0));
		double cosxal = sqrt(1 - sinxal2);
		double xal = asin(sinxal);
		complex<double> Xasinxbe = x(1, 0) / cosxal;
		double xA = std::arg(Xasinxbe);
		double sinxbe = std::abs(Xasinxbe);
		double sinxbe2 = std::norm(Xasinxbe);
		double xbe = asin(sinxbe);
//		double cosxbe = sqrt(1 - sinxbe2);/		double cosxal = x(0, 0) / cosxbe;
///////	for y

		double sinyal = std::abs(y(2, 0));
		double sinyal2 = std::norm(y(2, 0));
		double yB = std::arg(x(2, 0));
		double cosyal = sqrt(1 - sinyal2);
		double yal = asin(sinyal);
		complex<double> Yasinybe = y(1, 0) / cosyal;
		double yA = std::arg(Yasinybe);
		double sinybe = std::abs(Yasinybe);
		double sinybe2 = std::norm(Yasinybe);
		double ybe = asin(sinybe);
//		double cosybe = sqrt(1 - sinybe2);
//		double cosyal = y(0, 0) / cosybe;

		ret[0] = xal;
		ret[1] = xbe;
		ret[2] = xA;
		ret[3] = xB;

		ret[4] = yal;
		ret[5] = ybe;
		ret[6] = yA;
		ret[7] = yB;
	}
	if (p->d == 4)
	{
		printf("here4\n");
		ret = (double *)malloc(sizeof(double) * 6 * 2);
	double sinxga = std::abs(x(3, 0));
	double sinxga2 = std::norm(x(3, 0));
	double xga = asin(sinxga);
	double xC = std::arg(x(3, 0));

	double cosxga = sqrt(1 - sinxga2);
	complex<double> Xbsinxal = x(2, 0) / cosxga;
	
	double sinxal = std::abs(Xbsinxal);
	double sinxal2 = std::norm(Xbsinxal);
	double xal = asin(sinxal);
	double xB = std::arg(x(2, 0));

	double cosxal = sqrt(1 - sinxal2);
	complex<double> Xasinxbe = x(1, 0) / (cosxal * cosxga);
	double sinxbe = std::abs(Xasinxbe);
	double sinxbe2 = std::norm(Xasinxbe);
	double xbe = asin(sinxbe);
	double xA = std::arg(x(1, 0)); // Xa here

///	now same thing for y

	double sinyga = std::abs(y(3, 0));
	double sinyga2 = std::norm(y(3, 0));
	double yga = asin(sinyga);
	double yC = std::arg(y(3, 0));

	double cosyga = sqrt(1 - sinyga2);
	complex<double> Ybsinyal = y(2, 0) / cosyga;
	
	double sinyal = std::abs(Ybsinyal);
	double sinyal2 = std::norm(Ybsinyal);
	double yal = asin(sinyal);
	double yB = std::arg(y(2, 0));

	double cosyal = sqrt(1 - sinyal2);
	complex<double> Yasinybe = y(1, 0) / (cosyal * cosyga);
	double sinybe = std::abs(Yasinybe);
	double sinybe2 = std::norm(Yasinybe);
	double ybe = asin(sinybe);
	double yA = std::arg(y(1, 0));

/*
		xal = gsl_vector_get(x, 0 + param->line_length * i);
		xbe =  gsl_vector_get(x, 1 + param->line_length * i);
		xga =  gsl_vector_get(x, 2 + param->line_length * i);
		xA =  gsl_vector_get(x, 3 + param->line_length * i);
		xB =  gsl_vector_get(x, 4 + param->line_length * i);
		xC =  gsl_vector_get(x, 5 + param->line_length * i);

		yal = gsl_vector_get(x, 6 + param->line_length * i);
		ybe = gsl_vector_get(x, 7 + param->line_length * i);
		yga = gsl_vector_get(x, 8 + param->line_length * i);
		yA = gsl_vector_get(x, 9 + param->line_length * i);
		yB = gsl_vector_get(x, 10 + param->line_length * i);
		yC = gsl_vector_get(x, 11 + param->line_length * i);
*/

//	x vector params
	ret[0] = xal;
	ret[1] = xbe;
	ret[2] = xga;
	ret[3] = xA;
	ret[4] = xB;
	ret[5] = xC;

//	now y
	ret[6] = yal;
	ret[7] = ybe;
	ret[8] = yga;
	ret[9] = yA;
	ret[10] = yB;
	ret[11] = yC;

	printf("xal, xbe, xga, Xa, Xb, Xc: %f, %f, %f, %f, %f, %f\n", xal, xbe, xga, xA, xB, xC);
	printf("yal, ybe, yga, Ya, Yb, Yc: %f, %f, %f, %f, %f, %f\n", yal, ybe, yga, yA, yB, yC);


	}
	return (ret);

}


//		this function takes all pairs of vectors written to v and extracts angles.
//		write ith pair here 	//guess in vector form	//number of the pair of vectors
void	init_minimizing_guess(t_param *p, Eigen::MatrixXcd *v)
{
	int ll = p->line_length;
	double *vector_params;
	printf("here\n");
//	for (int j = 0; j < ll; j++)
//		vector_params[j] = get_angles_form_

//	line_length  = 2 * degrees
//	ith x vector
	int number_of_lines = p->number_of_terms;
	for (int i = 0; i <  number_of_lines; i++)
	{
	vector_params = get_angles_from_vector(p, v[2 * i], v[2 * i + 1]);
	printf("here5\n");
	for (int j = 0; j < ll; j++)
	{
		printf("i j %d, %d\n", i, j);
//								ith line
//		printf("setting param %f\n", vector_params[j]);
		gsl_vector_set(p->v, j + ll * i, vector_params[j]);
		printf("here7\n");
//		gsl_vector_set(p->x, j + n + param->line_length * i, vector_params[j + n]); // inits y params
	}
		printf("\n");
	}

/*
	gsl_vector_set (x, 0, 0.17278);		// param for x0
	gsl_vector_set (x, 1, M_PI / 4.0);	// param for x0
	gsl_vector_set (x, 2, 0.6936);		// param for y0
	gsl_vector_set (x, 3, M_PI / 4.0);	// param for y0
*/

//	for d = 2 structure of p->x
//	x1_param1, x1_param2, y1_param1, y1_param2,
//	x2_param1, x2_param2, y2_param1, y2_param2,
//	x3_param1, x3_param2, y3_param1, y3_param2,
//	x4_param1, x4_param2, y4_param1, y4_param2.
	if (p->d == 2 || p->d == 3 || p->d == 4)
	{
		printf("freeing vector params\n");
		free(vector_params);
	}
}

///////////////////////#endif

////////////////////////////////////////////////////


/*

template <typename T> double sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

*/
/*function to minimize*/
//	v = the argument of the fucntion
//	params contain allocated memory for writing the resulting matrices

double	my_f (const gsl_vector *v, void *params)
{
	t_param *p = (t_param *)params;

	double t = p->t;
	double d = p->d;
	int d2 = p->d2;

	int degrees_in_vector = p->degrees_in_vector;
	int number_of_terms = p->number_of_terms;
	int number_of_weights = p->number_of_weights;
	
	int degrees_of_freedom = p->degrees_of_freedom;
	int line_length = p->line_length; // length of parameters pertaining to the first term (first weight)
	int lines = p->lines;


	double arg[number_of_terms][line_length];
	for (int i = 0; i < number_of_terms; i++) {
		for (int j = 0; j < line_length; j++) {
			arg[i][j] = gsl_vector_get (v, i * line_length + j);
//			printf("arg[%d,%d] is %f\n", i, j, arg[i][j]);
		}
	}


//	the last weight is rewritten to be this,
//	can be problems with this as the absolute value of last weight may be bigger than 1
//	std::cout << "calculating my_f(" << v << ")" << std::endl;
//	std::cout << "t =  " << t << ", d = " << d << std::endl;

//	arg[i][0] i = 0,1,2... number_of_weights - 1 are all number_of_weights weights
//	arg[number_of_weights - 1][0] = 1;

/*	for (int i = 0; i < number_of_weights - 1; i++)
		arg[number_of_weights - 1][0] -= abs(arg[i][0]);
*/

/*	if (arg[number_of_weights - 1][0] < 0)
		printf("POTENTIAL PROBLEM: weight < 0\n");
*/

/*	std::complex<double> Mx[4][2][2];
	std::complex<double> My[4][2][2];
*/
//	std::complex<double> fullx[4];
//	std::complex<double> fully[4];


//	MyMatrixXd Mx;
//	MyMatrixXd My;


//	for (int i = 0; i < number_of_terms; i++)
//	{
//		Mx = MyMatrixXd(d,d);
//		My = MyMatrixXd(d,d);
//	}



	double B, A, Y, X, w;

	//MyMatrixXd MxMy;
//	Eigen::MatrixXcd M(d2,d2);


//	printf("here\n");
	p->M = p->Zero;

//	std::cout << p->M << std::endl;
//	std::cout << p->Identity << std::endl;

	for (int i = 0; i < number_of_terms; i++)
	{
//		|x><x|_i
//		w = abs(arg[i][0]);
		w = 1 / (double)p->number_of_terms;

//		printf("init xy number %d\n", i);

//		here only 1 pair of vectors is calculated from d^2 pairs.
		init_xy(p, v, i);

//		printf("here3\n");
#if 0

/*
		MyMatrixXd x(d,1);
		MyMatrixXd y(d,1);

//		for future transposition
		MyMatrixXd xt(d,1);
		MyMatrixXd yt(d,1);
*/

		int DELETE_ONE_TEST;



		double xal, xbe, xga, xA, xB, xC;

		double yal, ybe, yga, yA, yB, yC;

//		поміняти починаючи з 1;

		xal = arg[i][0];
		xbe = arg[i][1];
		xga = arg[i][2];
		xA = arg[i][3];
		xB = arg[i][4];
		xC = arg[i][5];

		yal = arg[i][6];
		ybe = arg[i][7];
		yga = arg[i][8];
		yA = arg[i][9];
		yB = arg[i][10];
		yC = arg[i][11];




		p->x(0, 0) = cos(xal) * cos(xbe) * cos(xga);
		p->x(1, 0) = polar(1.0, xA) * cos(xal) * sin(xbe) * cos(xga);
		p->x(2, 0) = polar(1.0, xB) * sin(xal) * cos(xga);
		p->x(3, 0) = polar(1.0, xC) * sin(xga);


		p->y(0, 0) = cos(yal) * cos(ybe) * cos(yga);
		p->y(1, 0) = polar(1.0, yA) * cos(yal) * sin(ybe) * cos(yga);
		p->y(2, 0) = polar(1.0, yB) * sin(yal) * cos(yga);
		p->y(3, 0) = polar(1.0, yC) * sin(yga);

# endif
/*
		p->xt = p->x;
		p->yt = p->y;

//		xt.transposeInPlace();
//		xt.conjugate();

		p->xt.adjointInPlace() ;

//		std::cout << "x is: \n" << x << std::endl;
//		std::cout << "transpose conjugate is \n" << xt << std::endl;

//		yt.transposeInPlace();
//		yt.conjugate();
		p->yt.adjointInPlace();
*/

//		Mx = |x><x|
//		My = |y><y|

		p->x.Kron(p->xt, p->Mx);	
		p->y.Kron(p->yt, p->My);

//		std::cout << "Matrix Mx was" << p->Mx << std::endl;
		p->Mx *= w;

//		std::cout << "Matrix Mx is" << p->Mx << std::endl;
/*
		Mx[i](0,0) = w * Mx[i](0,0);
		Mx[i](1,0) = w * Mx[i](1,0);
		Mx[i](0,1) = w * Mx[i](0,1);
		Mx[i](1,1) = w * Mx[i](1,1);
*/

//		std::cout << "Mx[" << i << "] is " << std::endl << Mx[i] << std::endl;

//		fullx[i] = {x[0],x[1]};
//		fully[i] = {y[0], y[1]};



//#if 0
/*
		Mx[i](0, 0) = w * x(0, 0) * x(0, 0);
		Mx[i](0, 1) = w * x(0, 0) * conj(x(1, 0));
		Mx[i](1, 0) = w * x(1, 0) * x(0, 0);
		Mx[i](1, 1) = w * x(1, 0) * conj(x(1, 0));
*/

/*		std::cout << " now Mx is " << std::endl;

		std::cout << Mx[i](0, 0) << Mx[i](0, 1) << endl;
		std::cout << Mx[i](1, 0) << Mx[i](1, 1) << endl;
*/


//		|y><y|_i
/*
		My[i](0, 0) = y(0, 0) * y(0, 0);
		My[i](0, 1) = y(0, 0) * conj(y(1, 0));
		My[i](1, 0) = y(1, 0) * y(0, 0);
		My[i](1, 1) = y(1, 0) * conj(y(1, 0));
*/

//#endif

#if 0
		Mx[i](0, 0) = w * x[0] * x[0];
		Mx[i](1, 0) = w * x[0] * conj(x[1]);
//		Mx[i][1][0] = w * sgn(cos(Y) * sin(Y)) * polar(abs(cos(Y) * sin(Y)), -X);
		Mx[i](0, 1) = w * x[1] * x[0];
		Mx[i](1, 1) = w * x[1] * conj(x[1]);


//		|y><y|_i
		My[i](0, 0) = y[0] * y[0];
		My[i](1, 0) = y[0] * conj(y[1]);
		My[i](0, 1) = y[1] * y[0];
		My[i](1, 1) = y[1] * conj(y[1]);

#endif

//		MyMatrixXd ytn = dynamic_cast<MyMatrixXd *>(&yt);


//		std::cout << "x * x is: " << Mx[i] << std::endl;


/*		printf("<x%d|y%d> = %f\n", i, i, sqrt(norm(x[0] * y[0] + conj(x[1]) * y[1])));
		printf("sqrt((t * (d^2-1) + 1)/d) = %f\n",sqrt((t * (d*d-1) + 1)/d));
*/
/*		std::cout << "x" << i << " is [" << x[0] << ", " << x[1] << "]" << endl;
		std::cout << "y" << i << " is [" << y[0] << ", " << y[1] << "]" << endl;
*/

		p->Mx.Kron(p->My, p->MxMy);
		p->M += p->MxMy;

	}
/*
	printf("|<x0| y1>|^2 is %f\n", norm(conj(fullx[0]) * fully[1]));
	printf("(1 - t) / d = %f\n", (1 - t)/ d);
*/
//	std::cout << endl;

//	std::complex<double> MxMy[4][4][4];
/*
	for (int i = 0; i < number_of_terms; i++)
		MxMy[i] = MyMatrixXd(d2,d2);
*/


/*
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{*/
			/*for (int k = 0; k < d2; k++)
			{
				M += MxMy[k];
			}*/
		/*}
	}*/
//	std::cout << "M now is " << M << endl;

/*	for (int j = 0; j < 4; j++)
	{
		for (int k = 0; k < 4; k++)
		{
			M[j][k] += MxMy[i][j][k];
		}
	}
*/
//	std::cout << MyMatrixXd::Zero(d2,d2);


/*
	MyMatrixXd Usw;

	Usw = MyMatrixXd(d2, d2);
*/

//	printf("______Usw = \n");
//	Usw = Eigen::MatrixXcd::Zero(d2, d2);
//	std::cout << Usw << std::endl;


//	first coordinate is y
	

/*
	Usw(0, 0) = 1;
	Usw(0, 1) = 0;
	Usw(0, 2) = 0;
	Usw(0, 3) = 0;

	Usw(1, 0) = 0;
	Usw(1, 1) = 0;
	Usw(1, 2) = 1;
	Usw(1, 3) = 0;

	Usw(2, 0) = 0;
	Usw(2, 1) = 1;
	Usw(2, 2) = 0;
	Usw(2, 3) = 0;


	Usw(3, 0) = 0;
	Usw(3, 1) = 0;
	Usw(3, 2) = 0;
	Usw(3, 3) = 1;

*/

//	std::cout << "Usw:" << std::endl << Usw << std::endl;

//	for (int j = 0; j < 4; j++)
///	{
//		for (int k = 0; k < 4; k++)
//		{
//			std::cout << p->Usw << std::endl << std::endl;
			p->M -= (t / (double)d) * p->Usw + (1-t) / (double)(d*d) * p->Identity;

//		}
//	}

//	printf(" t / d is %f\n", (t / (double)d));
//	std::cout << "M[0][0] - t / d is " << M[0][0] << endl;


//	for (int k = 0; k < 4; k++)
//	{
//		std::cout << "M was" << std::endl << p->M << std::endl;

//			std::cout << (1 - t) / (double)(d * d) * p->Identity << std::endl << std::endl;
//		p->M-= (1-t) / (double)(d*d) * p->Identity;

//	}
//	std::cout << "Matrix is \n" << M <<std::endl;	
	double n = 0;

		for (int j = 0; j < d2; j++)
		{
			for (int k = 0; k < d2; k++)
			{
	//			printf("M[%d][%d] = ", j, k);
	//			std::cout << M[j][k] << endl;
				n += norm(p->M(j,k));
			}
		}
//	printf("my_f: %f\n", n);
	return (n);
}



void	init_Usw(t_param *param)
{


for (int j = 0; j < param->d2; j++)
{
	for (int i = 0; i < param->d2; i++)
	{
		param->Usw(j, i) = MyMatrixXd::Zero(param->d2, param->d2)(j, i);
	}
}

MyMatrixXd Basis(param->d, param->d);

/*	std::cout << "Identity:" << std::endl;
std::cout << Eigen::MatrixXd::Identity(d, d);
std::cout << std::endl;
*/

//	Basis = MyMatrixXd::Identity(d, d);

for (int j = 0; j < param->d; j++)
{
	for (int i = 0; i < param->d; i++)
	{
		Basis(j, i) = MyMatrixXd::Identity(param->d, param->d)(j, i);
	}
}

/*	std::cout << std::endl;
std::cout << "Basis:" << std::endl;
std::cout << Basis << std::endl;

*/
//#if 0

for (int j = 0; j < param->d; j++)
{
	for (int i = 0; i < param->d; i++)
	{
//			std::cout << "i j =" << std::endl << i << " " << j << std::endl << std::endl;

		MyMatrixXd vi(param->d,1);
		MyMatrixXd vj(param->d,1);
		MyMatrixXd vit(param->d,1);
		MyMatrixXd vjt(param->d,1);

		for (int k = 0; k < param->d; k++)
			vi(k) = Basis.col(i)(k);

		for (int k = 0; k < param->d; k++)
			vj(k) = Basis.col(j)(k);


//			std::cout << "vj" << std::endl << vj << std::endl;

		vit = vi;
		vjt = vj;
		vit.adjointInPlace();
		vjt.adjointInPlace();

		MyMatrixXd vjvi(param->d2, 1);
		vj.Kron(vi, vjvi);

//			std::cout << "vjvi" << std::endl << vjvi << std::endl;
		
		MyMatrixXd vivjt(param->d2, 1);

		vi.Kron(vj, vivjt);
		
		vivjt.adjointInPlace();
		

//			std::cout << "vivjt" << std::endl << vivjt << std::endl;
		

		MyMatrixXd term(param->d2, param->d2);
		vjvi.Kron(vivjt, term);

		param->Usw = param->Usw + term;
	}
}
}

Eigen::MatrixXcd        pow(Eigen::MatrixXcd &M, int power)
{

        Eigen::MatrixXcd ret = Eigen::MatrixXcd(M.cols(), M.rows());
        ret = Eigen::MatrixXcd::Identity(M.cols(), M.rows());
        for (int i = 0; i < power; i++)
        {
                ret = ret * M;
        }
        return (ret);
}


void	init_weyl(t_param *p)
{
//	works for 4x4

//#if 0
        Eigen::MatrixXcd S;
        Eigen::MatrixXcd C;

        S = Eigen::MatrixXcd(p->d,p->d);
        C = Eigen::MatrixXcd(p->d,p->d);

	C = Eigen::MatrixXcd::Zero(p->d, p->d);
	S = Eigen::MatrixXcd::Zero(p->d, p->d);

 
/*
        C(0,0) = 0;
        C(1,0) = 0;
        C(2,0) = 0;
        C(3,0) = 0;


	C(0,1) = 0;
        C(1,1) = 0;
        C(2,1) = 0;
        C(3,1) = 0;

	C(0,2) = 0;
        C(1,2) = 0;
        C(2,2) = 0;
        C(3,2) = 0;

	C(0,3) = 0;
        C(1,3) = 0;
        C(2,3) = 0;
        C(3,3) = 0;

*/
	complex<double>T = polar(1.0,  2 * M_PI * (p->d + 1)/(2 * p->d));
	for (int i = 0; i < p->d; i++)
	{
		C(i ,i) = polar(1.0, 2 * M_PI / p->d * i);
	} 

	for(int i = 0; i < p->d; i++)
	{
			S((i+1) % p->d, i) = 1;
	}
/*

        S(0,0) = 0;
        S(1,0) = 1;
        S(2,0) = 0;
        S(3,0) = 0;


	S(0,1) = 0;
        S(1,1) = 0;
        S(2,1) = 1;
        S(3,1) = 0;

	S(0,2) = 0;
        S(1,2) = 0;
        S(2,2) = 0;
        S(3,2) = 1;

	S(0,3) = 1;
        S(1,3) = 0;
        S(2,3) = 0;
        S(3,3) = 0;
*/
	std::cout << "S is " << endl << S << std::endl;
	std::cout << "C is " << endl << C << std::endl;



        for (int i = 0; i < p->d; i++)
        {
                for (int j = 0; j < p->d; j++)
                {

/*
                        cout << "x is " << endl << x << endl;
                        cout << "y is " << endl << y << endl;
*/

                        p->W[i][j] = Eigen::MatrixXcd(p->d,p->d);
  /*
                      cout << C.pow(i) << endl << endl;
                        cout << C.pow(i) << endl << endl;

                        cout << pow(S, j) << endl << endl;
                        cout << pow(S, j) << endl << endl;

                        cout << pow(C,i) * pow(S,j) << endl << endl;
*/
                        p->W[i][j] = pow(pow(T, i), j) * pow(C, i) * pow(S, j);

			cout << "W" << i << j << endl;
                        cout << p->W[i][j] << endl << endl;

		}
	}
//#endif

}








int	main(void)
{
	srand(time(0));
	int d = 3;
	double eps = /*0.00001 */0.00001;
	double t =  0.5 * 1/( (double) (d + 1))/*0.2*/ /*- 0.0001*/;
//                       t   d
/*
	double par[2] = {t, (double)d};
*/
	printf("d = %d, t = %.20f\n", d, t);

	const gsl_multimin_fminimizer_type *T =
		gsl_multimin_fminimizer_nmsimplex2rand;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;

	int status;
	double size;

	t_param param;


printf("here1\n");
  MyMatrixXd m1(2,2);
  m1(0,0) = 1;
  m1(1,0) = 3;
  m1(0,1) = 2;
  m1(1,1) = 4;
	std::cout << m1 << std::endl << std::endl;
printf("here2\n");
  MyMatrixXd m2(2,2);
  m2(0,0) = 1;
  m2(1,0) = 1;
  m2(0,1) = 1;
  m2(1,1) = 1;

  MyMatrixXd m3(4,4);

printf("here\n");
	m1.Kron(m2, m3);


  std::cout << "kronecker result:" <<std::endl << m3 << std::endl;


MyMatrixXd v1(2,1);
MyMatrixXd v2(2,1);

v1(0,0) = 1;
v1(1,0) = 2;

	v2 = v1;
	v2.adjointInPlace();

	std::cout << "v1" << v1 << "v2" << v2 << std::endl;
MyMatrixXd result(2,2);

v1.Kron(v2,result);

	result = result * 0.25;

  std::cout << "kronecker result:" <<std::endl << result << std::endl;



	/* Starting point */
//	starting point is the setof all parameters that define vectors that define my_f() output
//	for d = 2:
//	4 weights + 4 * 2 complex vectors that each have 2 degrees of freedom. (4 degrees of freedom
//	minus an arbitrary phase multiplication minus the normal requirement.
//	weights are written to the first cooridnate of each of 4 sections of the
//	20-dimensional vector x


	param.degrees_in_vector = 2 * d - 2;
	param.number_of_terms = d * d;
	param.number_of_weights = d * d;
	
	param.degrees_of_freedom = param.number_of_terms * 2 * param.degrees_in_vector;
	param.d = d;
	param.d2 = d * d;

	param.t = t;
//	for d = 4:
//	16 weights + 16 * 2 complex vectors that each have 6 degrees of freedom

	printf("alloc\n");
	printf("degrees of freedom %d\n", param.degrees_of_freedom);
	x = gsl_vector_alloc(param.degrees_of_freedom);
	param.v = x;
//	NOTE: weights are not included now!
//	structure of x: weight1, x1_param1, x1_param2, ..., y1_param1, y1_param2, ...
//			weight2, x2_param1, x2_param2, ..., y2_param1, y2_param2, ...
//			line length = 2 * d - 2
//	for d = 2 structure of x

//	x1_param1, x1_param2, y1_param1, y1_param2,
//	x2_param1, x2_param2, y2_param1, y2_param2,
//	x3_param1, x3_param2, y3_param1, y3_param2,
//	x4_param1, x4_param2, y4_param1, y4_param2.

	param.line_length = 2 * param.degrees_in_vector; // length of parameters pertaining to the first term (first weight)
	printf("line length: %d\n", param.line_length);
	param.lines = param.number_of_terms;

	param.Mx = MyMatrixXd(d,d);
	param.My = MyMatrixXd(d,d);
	param.MxMy = MyMatrixXd(param.d2,param.d2);
	param.M = Eigen::MatrixXcd(param.d2,param.d2);
	param.Usw = MyMatrixXd(param.d2, param.d2);

	param.x = MyMatrixXd(d, 1);
	param.y = MyMatrixXd(d, 1);
	param.xt = MyMatrixXd(d, 1);
	param.yt = MyMatrixXd(d, 1);
	param.Zero = Eigen::MatrixXcd(param.d2, param.d2);
	std::cout << param.Zero << std::endl << std::endl;
	param.Identity = MyMatrixXd(param.d2, param.d2);
	param.Zero = Eigen::MatrixXcd::Zero(param.d2, param.d2);; 


	for (int i = 0; i < param.d2; i++)
		for (int j = 0; j < param.d2; j++)
			param.Identity(i, j) = MyMatrixXd::Identity(param.d2, param.d2)(i, j);

	std::cout << param.Identity << std::endl;

	init_Usw(&param);
	init_weyl(&param);
	std::cout << "Usw := " << std::endl << param.Usw << std::endl << std::endl;
//#endif



//	gsl_vector_set_all (x, 0);


/*
	gsl_vector_set (x, 0, 0.25);
	gsl_vector_set (x, 1, 0.17278);
	gsl_vector_set (x, 2, M_PI / 4.0);
	gsl_vector_set (x, 3, 0.6936);
	gsl_vector_set (x, 4, M_PI / 4.0);

	gsl_vector_set (x, 5, 0.25);
	gsl_vector_set (x, 6, 0.17278);
	gsl_vector_set (x, 7, 5 * M_PI / 4.0);
	gsl_vector_set (x, 8, 0.6936);
	gsl_vector_set (x, 9, 5 * M_PI / 4.0);

	gsl_vector_set (x, 10, 0.25);
	gsl_vector_set (x, 11, 1.398);
	gsl_vector_set (x, 12, -M_PI / 4.0);
	gsl_vector_set (x, 13, 0.8745);
	gsl_vector_set (x, 14, -M_PI / 4.0);


	gsl_vector_set (x, 15, 0.25);
	gsl_vector_set (x, 16, 1.398);
	gsl_vector_set (x, 17, -5 * M_PI / 4.0);
	gsl_vector_set (x, 18, 0.8745);
	gsl_vector_set (x, 19, -5 * M_PI / 4.0);
*/



//	for d = 2
//	for t = 0.5 * 1 / (d + 1) !!!

/*
	gsl_vector_set (x, 0, 0.17278);		// param for x0
	gsl_vector_set (x, 1, M_PI / 4.0);	// param for x0
	gsl_vector_set (x, 2, 0.6936);		// param for y0
	gsl_vector_set (x, 3, M_PI / 4.0);	// param for y0

	gsl_vector_set (x, 4, 0.17278);
	gsl_vector_set (x, 5, 5 * M_PI / 4.0);
	gsl_vector_set (x, 6, 0.6936);
	gsl_vector_set (x, 7, 5 * M_PI / 4.0);

	gsl_vector_set (x, 8, 1.398);
	gsl_vector_set (x, 9, -M_PI / 4.0);
	gsl_vector_set (x, 10, 0.8745);
	gsl_vector_set (x, 11, -M_PI / 4.0);


	gsl_vector_set (x, 12, 1.398);
	gsl_vector_set (x, 13, -5 * M_PI / 4.0);
	gsl_vector_set (x, 14, 0.8745);
	gsl_vector_set (x, 15, -5 * M_PI / 4.0);

*/

//		this function takes all pairs of vectors and extracts angles.
//		i ranges from 0 to 2 * number_of_terms
//		write ith pair here 	//guess in vector form	//number of the pair of vectors


	Eigen::MatrixXcd guess[param.number_of_terms * 2];
	for (int i = 0; i < param.number_of_terms * 2; i++)
		guess[i] = Eigen::MatrixXcd(param.d, 1);


//	for d = 2, t = 0.5 / (d+1)

/*
guess[0] << 1.00000000000000000000 + 0.00000000000000000000 * im, 0.00000000000000000000 + 0.00000000000000000000 * im;
guess[1] << 0.86602540378416570377 + 0.00000000000000000000 * im, 0.50000000000047273296 + 0.00000000000000000000 * im;

guess[2] << 0.22831118948877446306 + 0.00000000000000000000 * im, 0.71111332912290403652 + -0.66497506261360017632 * im;
guess[3] << 0.64549722436762535516 + 0.00000000000000000000 * im, 0.13338316976929021074 + -0.75202544063082621406 * im;

guess[4] << 0.94027515832131147722 + 0.00000000000000000000 * im, -0.33819389032882735124 + -0.03882678441599057145 * im;
guess[5] << 0.64549722436808332215 + 0.00000000000000000000 * im, -0.76368819500914941756 + -0.01066180738747262748 * im;

guess[6] << 0.25250074732206195804 + 0.00000000000000000000 * im, 0.61639494294341312663 + 0.74585564750532440392 * im;
guess[7] << 0.64549722436801748593 + 0.00000000000000000000 * im, -0.04051536800939155303 + 0.76268724801733034369 * im;
*/


//d = 2	// t = 1 / (d + 1) //|x0> = |0>

/*
guess[0] << 1.00000000000000000000 + 0.00000000000000000000 * im, 0.00000000000000000000 + 0.00000000000000000000 * im;
guess[1] << 0.99999999367380554283 + 0.00000000000000000000 * im, 0.00011248283829210840 + 0.00000000000000000000 * im;

guess[2] << 0.57725842734339416484 + 0.00000000000000000000 * im, 0.81656151517263630968 + 0.00000000000000000000 * im;
guess[3] << 0.57735027283339368509 + 0.00000000000000000000 * im, 0.81649657835119304750 + 0.00000000000338500860 * im;

guess[4] << 0.57739618463103348223 + 0.00000000000000000000 * im, -0.40818334848342413546 + -0.70710678118257719316 * im;
guess[5] << 0.57735027284097051314 + 0.00000000000000000000 * im, -0.40834570215653265279 + -0.70705052859095252060 * im;

guess[6] << 0.57739618463403541426 + 0.00000000000000000000 * im, -0.40818334848031445627 + 0.70710678118192105135 * im;
guess[7] << 0.57735027284868301045 + 0.00000000000000000000 * im, -0.40834570215230908685 + 0.70705052858709394048 * im;
*/



//PPPR solution

/*
guess[0] << 0.88807383625838032248 + 0.00000000000000000000 * im, 0.32505758055560213249 + 0.32505758055560207698 * im;
guess[1] << 0.88807383169585019100 + 0.00000000000000000000 * im, 0.32505758678813423401 + 0.32505758678813417850 * im;

guess[2] << 0.45970083897391733618 + -0.00000000000000002776 * im, 0.62796303181265245019 + -0.62796303181265233917 * im;
guess[3] << 0.45970084778804881642 + 0.00000000000000002776 * im, 0.62796302858645647316 + -0.62796302858645636213 * im;

guess[4] << 0.88807383625838032248 + 0.00000000000000000000 * im, -0.32505758055560218800 + -0.32505758055560202147 * im;
guess[5] << 0.88807383169585019100 + 0.00000000000000000000 * im, -0.32505758678813428952 + -0.32505758678813412299 * im;

guess[6] << 0.45970083897391733618 + -0.00000000000000002776 * im, -0.62796303181265233917 + 0.62796303181265245019 * im;
guess[7] << 0.45970084778804881642 + 0.00000000000000002776 * im, -0.62796302858645636213 + 0.62796302858645647316 * im;
*/

//just random solution for t = 1/(d + 1), d = 2

/*
guess[0] << 0.686330 + 0.000000 * im, -0.183671 + 0.703716 * im;
guess[1] << 0.686327 + 0.000000 * im, -0.183678 + 0.703717 * im;

guess[2] << 0.987190 + 0.000000 * im, -0.036161 + -0.155398 * im;
guess[3] << 0.987190 + 0.000000 * im, -0.036155 + -0.155400 * im;

guess[4] << 0.583574 + 0.000000 * im, 0.805660 + -0.101754 * im;
guess[5] << 0.583579 + 0.000000 * im, 0.805656 + -0.101757 * im;

guess[6] << 0.462438 + 0.000000 * im, -0.666911 + -0.584278 * im;
guess[7] << 0.462437 + 0.000000 * im, -0.666920 + -0.584269 * im;
*/


//all of this doesn't work why?
//	d = 2, t = 0.5 / (d+1)


/*
guess[0] << 1.00000000000000000000 + 0.00000000000000000000 * im, 0.00000000000000000000 + 0.00000000000000000000 * im;
guess[1] << 0.86602540378691417189 + 0.00000000000000000000 * im, -0.35258688501728469022 + 0.35451726122929511087 * im;

guess[2] << 0.19060670241694005478 + 0.00000000000000000000 * im, -0.35583460033892377883 + 0.91490481592096706276 * im;
guess[3] << 0.64549722436931178393 + 0.00000000000000000000 * im, -0.00897429900840247124 + 0.76370988947952100911 * im;

guess[4] << 0.90885382340961506920 + 0.00000000000000000000 * im, 0.41620926780469302830 + -0.02746949339186330330 * im;
guess[5] << 0.64549722436621725929 + 0.00000000000000000000 * im, 0.59128703051856745798 + -0.48343870436286351389 * im;

guess[6] << 0.37101726734477979974 + 0.00000000000000000000 * im, -0.83675060956037661253 + -0.40273391306463196537 * im;
guess[7] << 0.64549722436648815371 + 0.00000000000000000000 * im, -0.10926778578411425191 + -0.75590600230782434288 * im;
*/

//	d = 2, t = 0.5 / (d+1)
/*
guess[0] << 1.00000000000000000000 + 0.00000000000000000000 * im, 0.00000000000000000000 + 0.00000000000000000000 * im;
guess[1] << 0.86602540378526204901 + 0.00000000000000000000 * im, -0.27611587396076742174 + 0.41684532400574614286 * im;

guess[2] << 0.30561844789246245258 + 0.00000000000000000000 * im, -0.94839410170145033163 + -0.08453397048347624509 * im;
guess[3] << 0.64549722436874978904 + 0.00000000000000000000 * im, -0.43079357592479100569 + -0.63067442335500745187 * im;

guess[4] << 0.20354056624078911697 + 0.00000000000000000000 * im, -0.02165595409345681330 + 0.97882698039371640597 * im;
guess[5] << 0.64549722436665424308 + 0.00000000000000000000 * im, 0.34297692325718270867 + 0.68242227648867215262 * im;

guess[6] << 0.93014439857530939459 + 0.00000000000000000000 * im, 0.31635367472322617477 + -0.18641821340225075976 * im;
guess[7] << 0.64549722436702650086 + 0.00000000000000000000 * im, 0.45826497097502705280 + -0.61100454148207639093 * im;
*/
//}


//#if 0
//d = 4, t = 1/(d + 1)

/*

guess[0] << 0.48571220566371603455 + 0.00000000000000000000 * im, 0.60043410989507439712 + -0.44989580508021820293 * im, -0.00000022134485245055 + -0.20118866436074883675 * im, -0.39924505078166983019 + -0.03581631296003146697 * im;
guess[1] << 0.48571220121010350024 + 0.00000000000000000000 * im, 0.60043410581348644062 + -0.44989581300127767793 * im, -0.00000022237934508396 + -0.20118866594721568353 * im, -0.39924505331953308307 + -0.03581630508201213980 * im;

guess[2] << 0.40084839004280192754 + 0.00000000000000000000 * im, -0.48376950530221196622 + 0.04339815948262136514 * im, -0.55783371497324985011 + 0.50174587981883667087 * im, 0.01797588775844632253 + 0.20038393068439400158 * im;
guess[3] << 0.40084839076185951878 + 0.00000000000000000000 * im, -0.48376950778052724145 + 0.04339816116647154942 * im, -0.55783370775348195547 + 0.50174588283297227864 * im, 0.01797589266778668385 + 0.20038393500912993206 * im;

guess[4] << 0.20118874233297490139 + 0.00000000000000000000 * im, 0.03581622211567773995 + -0.39924509546129000048 * im, 0.00000002041508061884 + 0.48571219018709288484 * im, 0.44989592627981206396 + 0.60043398118563251487 * im;
guess[5] << 0.20118874107775422400 + 0.00000000000000000000 * im, 0.03581624656015507807 + -0.39924509218375048292 * im, -0.00000001141238429109 + 0.48571218941213500919 * im, 0.44989589033705501553 + 0.60043400988566419940 * im;

guess[6] << 0.75028487689048017906 + 0.00000000000000000000 * im, 0.12063960664885407803 + -0.16100589317735533590 * im, -0.29802940657205773123 + -0.26806351147867168994 * im, 0.38870266112947454706 + 0.29124982238982549676 * im;
guess[7] << 0.75028487820865574776 + 0.00000000000000000000 * im, 0.12063960939808400508 + -0.16100589136249779387 * im, -0.29802940809578004622 + -0.26806350792765271107 * im, 0.38870266238507222845 + 0.29124981889200313168 * im;

guess[8] << 0.48571219380388774844 + 0.00000000000000000000 * im, 0.44989590371547083514 + 0.60043412756094016736 * im, 0.00000004305203688834 + 0.20118835371355753283 * im, -0.03581610949384984866 + 0.39924510228822884805 * im;
guess[9] << 0.48571219400042237035 + 0.00000000000000000000 * im, 0.44989590336042289964 + 0.60043412575145038712 * im, 0.00000004278645691197 + 0.20118835338460328011 * im, -0.03581610840153907904 + 0.39924510543431307452 * im;

guess[10] << 0.40084848092066294178 + 0.00000000000000000000 * im, -0.04339784949562915106 + -0.48376969158366989188 * im, 0.55783375004053603607 + -0.50174560122379563420 * im, 0.20038394399166689630 + -0.01797613602590326237 * im;
guess[11] << 0.40084848492557317856 + 0.00000000000000000000 * im, -0.04339784756032747520 + -0.48376968884190801878 * im, 0.55783375282198799372 + -0.50174559875668145459 * im, 0.20038394181292937257 + -0.01797613201319174167 * im;

guess[12] << 0.20118853364528579974 + 0.00000000000000000000 * im, 0.39924487530314201056 + 0.03581647033650977446 * im, 0.00000008571153665469 + -0.48571222335154873306 * im, 0.60043442196400753641 + -0.44989557114302441976 * im;
guess[13] << 0.20118853178127371617 + 0.00000000000000000000 * im, 0.39924487794770807669 + 0.03581645634914758236 * im, 0.00000007107952503248 + -0.48571222481774622981 * im, 0.60043441163302080366 + -0.44989558294820253259 * im;

guess[14] << 0.75028485626045227086 + 0.00000000000000000000 * im, 0.16100598148416428446 + 0.12063979314791220543 * im, 0.29802912377355567086 + 0.26806342513827419172 * im, 0.29124998594534068364 + -0.38870276031245887260 * im;
guess[15] << 0.75028485920200982395 + 0.00000000000000000000 * im, 0.16100597785790843330 + 0.12063979507336797248 * im, 0.29802912023228111948 + 0.26806343003264371250 * im, 0.29124998573825428005 + -0.38870275503405943285 * im;

guess[16] << 0.48571220269765696953 + 0.00000000000000000000 * im, -0.60043436697082652742 + 0.44989556077439513126 * im, 0.00000009117750315309 + -0.20118855605868582459 * im, 0.39924499496206200533 + 0.03581634285045333377 * im;
guess[17] << 0.48571220505380680033 + 0.00000000000000000000 * im, -0.60043436085072376951 + 0.44989556934372321173 * im, 0.00000008441651774008 + -0.20118855306505578628 * im, 0.39924499397957158298 + 0.03581633362409541999 * im;

guess[18] << 0.40084832354166383128 + 0.00000000000000000000 * im, 0.48376955546817029807 + -0.04339784639564987778 * im, -0.55783373699014804359 + 0.50174587841619999740 * im, -0.01797614518618950719 + -0.20038392953626152360 * im;
guess[19] << 0.40084832273510950795 + 0.00000000000000000000 * im, 0.48376955301968593837 + -0.04339784754631269281 * im, -0.55783373848238715365 + 0.50174587947962356882 * im, -0.01797614128373924494 + -0.20038393034487450572 * im;

guess[20] << 0.20118848068884656599 + 0.00000000000000000000 * im, -0.03581604923688611053 + 0.39924525986914244369 * im, 0.00000031345099763014 + 0.48571222007887515648 * im, -0.44989597021201127580 + -0.60043391274988844319 * im;
guess[21] << 0.20118848197181199433 + 0.00000000000000000000 * im, -0.03581603455022940724 + 0.39924526139794436208 * im, 0.00000031960346425786 + 0.48571221684646787020 * im, -0.44989598155194249207 + -0.60043390629749049392 * im;

guess[22] << 0.75028474675292156082 + 0.00000000000000000000 * im, -0.12063972844241556415 + 0.16100607483953649490 * im, -0.29802921571501156395 + -0.26806344955075805947 * im, -0.38870281269331680152 + -0.29125005678206761228 * im;
guess[23] << 0.75028474365869823881 + 0.00000000000000000000 * im, -0.12063973233339553293 + 0.16100607753154846935 * im, -0.29802921026336487431 + -0.26806345565894890148 * im, -0.38870281080943097640 + -0.29125006412404236267 * im;

guess[24] << 0.48571225269390266854 + 0.00000000000000000000 * im, -0.44989568139662694524 + -0.60043406932338205806 * im, -0.00000007588617055007 + 0.20118877391789016795 * im, 0.03581612501103648843 + -0.39924515560999279673 * im;
guess[25] << 0.48571225353689073678 + 0.00000000000000000000 * im, -0.44989567471439811941 + -0.60043407527298398030 * im, -0.00000008083198722739 + 0.20118877036222962418 * im, 0.03581612188749126041 + -0.39924515523865272693 * im;

guess[26] << 0.40084836728471889833 + 0.00000000000000000000 * im, 0.04339787420051807731 + 0.48376943806260669367 * im, 0.55783392561723277314 + -0.50174579505389627077 * im, -0.20038381895738816008 + 0.01797596820465803893 * im;
guess[27] << 0.40084836924778310951 + 0.00000000000000000000 * im, 0.04339787375129934566 + 0.48376943542497313766 * im, 0.55783392423412647698 + -0.50174579944222597039 * im, -0.20038381444482139537 + 0.01797596723521816334 * im;

guess[28] << 0.20118859529337815295 + 0.00000000000000000000 * im, -0.39924505716996455673 + -0.03581633111751213877 * im, 0.00000002314895719161 + -0.48571222295692917381 * im, -0.60043415991665505249 + 0.44989574342277172114 * im;
guess[29] << 0.20118859448170073323 + 0.00000000000000000000 * im, -0.39924505621703987313 + -0.03581633295774085818 * im, 0.00000002263994770294 + -0.48571222407816527777 * im, -0.60043416240736358080 + 0.44989573995026904685 * im;

guess[30] << 0.75028494390333944075 + 0.00000000000000000000 * im, -0.16100577850080066700 + -0.12063966592206552819 * im, 0.29802919808529487744 + 0.26806336677886088982 * im, -0.29124996985642465086 + 0.38870271003168499480 * im;
guess[31] << 0.75028494303019888090 + 0.00000000000000000000 * im, -0.16100577359455892079 + -0.12063966657294182550 * im, 0.29802920358248946586 + 0.26806337041367173102 * im, -0.29124996718138052021 + 0.38870270883009400142 * im;
*/

//d = 3 t = 1/(d + 1)??


/*
guess[0] << 0.81649658092772603446 + 0.00000000000000000000 * im, -0.40824829046386307274 + 0.00000000000000000000 * im, -0.40824829046386307274 + 0.00000000000000000000 * im;
guess[1] << 0.81649658092772603446 + 0.00000000000000000000 * im, -0.40824829046386307274 + 0.00000000000000000000 * im, -0.40824829046386307274 + 0.00000000000000000000 * im;

guess[2] << 0.40824829046386307274 + 0.00000000000000005000 * im, -0.81649658092772603446 + -0.00000000000000009999 * im, 0.40824829046386307274 + 0.00000000000000005000 * im;
guess[3] << 0.40824829046386307274 + 0.00000000000000005000 * im, -0.81649658092772603446 + -0.00000000000000009999 * im, 0.40824829046386307274 + 0.00000000000000005000 * im;

guess[4] << 0.40824829046386307274 + 0.00000000000000005000 * im, 0.40824829046386307274 + 0.00000000000000005000 * im, -0.81649658092772603446 + -0.00000000000000009999 * im;
guess[5] << 0.40824829046386307274 + 0.00000000000000005000 * im, 0.40824829046386307274 + 0.00000000000000005000 * im, -0.81649658092772603446 + -0.00000000000000009999 * im;

guess[6] << 0.81649658092772603446 + 0.00000000000000000000 * im, 0.20412414523193145310 + -0.35355339059327389739 * im, 0.20412414523193173066 + 0.35355339059327367535 * im;
guess[7] << 0.81649658092772603446 + 0.00000000000000000000 * im, 0.20412414523193145310 + -0.35355339059327389739 * im, 0.20412414523193173066 + 0.35355339059327367535 * im;

guess[8] << 0.40824829046386307274 + 0.00000000000000005000 * im, 0.40824829046386296172 + -0.70710678118654768376 * im, -0.20412414523193167515 + -0.35355339059327367535 * im;
guess[9] << 0.40824829046386307274 + 0.00000000000000005000 * im, 0.40824829046386296172 + -0.70710678118654768376 * im, -0.20412414523193167515 + -0.35355339059327367535 * im;

guess[10] << 0.40824829046386307274 + 0.00000000000000005000 * im, -0.20412414523193150862 + 0.35355339059327389739 * im, 0.40824829046386329479 + 0.70710678118654735069 * im;
guess[11] << 0.40824829046386307274 + 0.00000000000000005000 * im, -0.20412414523193150862 + 0.35355339059327389739 * im, 0.40824829046386329479 + 0.70710678118654735069 * im;

guess[12] << 0.81649658092772603446 + 0.00000000000000000000 * im, 0.20412414523193175842 + 0.35355339059327373086 * im, 0.20412414523193117555 + -0.35355339059327400841 * im;
guess[13] << 0.81649658092772603446 + 0.00000000000000000000 * im, 0.20412414523193175842 + 0.35355339059327373086 * im, 0.20412414523193117555 + -0.35355339059327400841 * im;

guess[14] << 0.40824829046386307274 + 0.00000000000000005000 * im, 0.40824829046386335030 + 0.70710678118654735069 * im, -0.20412414523193123106 + 0.35355339059327400841 * im;
guess[15] << 0.40824829046386307274 + 0.00000000000000005000 * im, 0.40824829046386335030 + 0.70710678118654735069 * im, -0.20412414523193123106 + 0.35355339059327400841 * im;

guess[16] << 0.40824829046386307274 + 0.00000000000000005000 * im, -0.20412414523193170290 + -0.35355339059327373086 * im, 0.40824829046386240661 + -0.70710678118654790580 * im;
guess[17] << 0.40824829046386307274 + 0.00000000000000005000 * im, -0.20412414523193170290 + -0.35355339059327373086 * im, 0.40824829046386240661 + -0.70710678118654790580 * im;
*/


// d = 3, t = 0.5 * 1/(d + 1)




guess[0] << 0.97003018089340831143 + 0.00000000000000000000 * im, -0.17181595990463416346 + 0.00000000000000000000 * im, -0.17181595990463416346 + 0.00000000000000000000 * im;
guess[1] << 0.65173918228522265128 + 0.00000000000000000000 * im, -0.53630030685903462562 + 0.00000000000000000000 * im, -0.53630030685903462562 + 0.00000000000000000000 * im;

guess[2] << 0.17181595990463416346 + 0.00000000000000002104 * im, -0.97003018089340831143 + -0.00000000000000011879 * im, 0.17181595990463416346 + 0.00000000000000002104 * im;
guess[3] << 0.53630030685903462562 + 0.00000000000000006568 * im, -0.65173918228522265128 + -0.00000000000000007982 * im, 0.53630030685903462562 + 0.00000000000000006568 * im;

guess[4] << 0.17181595990463416346 + 0.00000000000000002104 * im, 0.17181595990463416346 + 0.00000000000000002104 * im, -0.97003018089340831143 + -0.00000000000000011879 * im;
guess[5] << 0.53630030685903462562 + 0.00000000000000006568 * im, 0.53630030685903462562 + 0.00000000000000006568 * im, -0.65173918228522265128 + -0.00000000000000007982 * im;

guess[6] << 0.97003018089340831143 + 0.00000000000000000000 * im, 0.08590797995231704010 + -0.14879698605302174585 * im, 0.08590797995231715112 + 0.14879698605302166259 * im;
guess[7] << 0.65173918228522265128 + 0.00000000000000000000 * im, 0.26815015342951720179 + -0.46444968979731388048 * im, 0.26815015342951753485 + 0.46444968979731365843 * im;

guess[8] << 0.17181595990463416346 + 0.00000000000000002104 * im, 0.48501509044670404469 + -0.84007077909130600801 * im, -0.08590797995231713724 + -0.14879698605302166259 * im;
guess[9] << 0.53630030685903462562 + 0.00000000000000006568 * im, 0.32586959114261121462 + -0.56442268850069987618 * im, -0.26815015342951747934 + -0.46444968979731371395 * im;

guess[10] << 0.17181595990463416346 + 0.00000000000000002104 * im, -0.08590797995231705397 + 0.14879698605302174585 * im, 0.48501509044670448878 + 0.84007077909130589699 * im;
guess[11] << 0.53630030685903462562 + 0.00000000000000006568 * im, -0.26815015342951725730 + 0.46444968979731382497 * im, 0.32586959114261154768 + 0.56442268850069965413 * im;

guess[12] << 0.97003018089340831143 + 0.00000000000000000000 * im, 0.08590797995231717887 + 0.14879698605302169034 * im, 0.08590797995231692907 + -0.14879698605302180137 * im;
guess[13] << 0.65173918228522265128 + 0.00000000000000000000 * im, 0.26815015342951759036 + 0.46444968979731371395 * im, 0.26815015342951681321 + -0.46444968979731410252 * im;

guess[14] << 0.17181595990463416346 + 0.00000000000000002104 * im, 0.48501509044670459980 + 0.84007077909130589699 * im, -0.08590797995231694295 + 0.14879698605302180137 * im;
guess[15] << 0.53630030685903462562 + 0.00000000000000006568 * im, 0.32586959114261165871 + 0.56442268850069965413 * im, -0.26815015342951686872 + 0.46444968979731404701 * im;

guess[16] << 0.17181595990463416346 + 0.00000000000000002104 * im, -0.08590797995231716500 + -0.14879698605302169034 * im, 0.48501509044670337856 + -0.84007077909130634108 * im;
guess[17] << 0.53630030685903462562 + 0.00000000000000006568 * im, -0.26815015342951753485 + -0.46444968979731376946 * im, 0.32586959114261082604 + -0.56442268850070009822 * im;


	printf("guess inited\n");



	init_minimizing_guess(&param, guess);

//#endif


/*
	for (int i = 0; i < param.degrees_of_freedom; i++)
		gsl_vector_set (x, i,rand() / (float)RAND_MAX * M_PI * 2);

*/


	
	printf("vector set\n");
//	gsl_vector_set_all (x, 0);
	



	/* Set initial step sizes to 1 */
	ss = gsl_vector_alloc(param.degrees_of_freedom);
	float prec = 1e-8; /*1e-9; *//*1e-11;*/


//	gsl_vector_set_all(ss, prec);
	float step_size = 0.000001;
	gsl_vector_set_all(ss, step_size);

	/* Initialize method and iterate */
	minex_func.n = param.degrees_of_freedom;
	minex_func.f = my_f;
	minex_func.params = &param;


	s = gsl_multimin_fminimizer_alloc (T, param.degrees_of_freedom);
	printf("gsl multimin fminimizer set\n");
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	double arg[param.number_of_terms][param.line_length];

   	FILE *fptr;

   // use appropriate location if you are using MacOS or Linux
	fptr = fopen("./program.txt","w");

	printf("start minimizing\n");
	printf("CHECK t = %.20f, f = %.20f\n", t, my_f(x, &param));
//	int  = lround(1/(eps * (d + 1)));

	for (int i = 0; param.t > 0; i++)
	{
	size_t iter = 0;
/*
	for (int i = 0; i < param.degrees_of_freedom; i++)
		gsl_vector_set (x, i,rand() / (float)RAND_MAX * M_PI * 2);
*/


	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size(s);

		status = gsl_multimin_test_size (size, prec);


		if (status == GSL_SUCCESS)
		{
#if 0
				printf ("converged to minimum at %zu:\n", iter);

//#if 0

//				printf("print params\n");
				for (int i = 0; i < param.degrees_of_freedom; i++)
				{
/*					if (i == param.degrees_of_freedom - param.line_length)
					{
						float ret = 1;
						for (int i = 0; i < param.lines - 1; i++)
							ret -= abs(gsl_vector_get (s->x, i * param.line_length)); 
						printf("\n%f, ", ret);
					}
					else if (i % param.line_length == 0)
						printf("\n%f, ", abs(gsl_vector_get(s->x, i)));
*/
					if (i % param.line_length == 0)
						printf("\n");
//					else
						printf("%f, ", gsl_vector_get(s->x, i));
				}


				printf("\n f() = %f size = %f\n", s->fval, size);
#endif
/*				for (int i = 0; i < param.number_of_terms; i++) {
					for (int j = 0; j < param.line_length; j++) {
						arg[i][j] = gsl_vector_get (s->x, i * param.line_length + j);
					}
				}
*/

//						printf("arg[%d,%d] is %f\n", i, j, arg[i][j]);
//#endif



		int condition;
		condition = 1;

#if 0
		for (int i = 0; i < param.number_of_terms; i++)
		{
			// this checks all pairwise products <xi|yi> if the code found a local minimum
			// to check of it may be global
			init_xy(&param, s->x, i);

/*
			std::complex<double> xy2 = (t * (param.d2 - 1) + 1) / (double)param.d;
			double value = norm(xy2 - (param.xt * param.y)(0,0).real() * (param.xt * param.y)(0,0) );
			std::cout << std::endl << "(...) - |<xi|yi>|^2 = " << value;
*/

///////////////////////////////////////////////////////////////////
			double xy2 = (t * (param.d2 - 1) + 1) / (double)param.d;
			std::complex<double> xty2 = std::norm((param.xt * param.y)(0,0)); 
			double value = xy2 - xty2.real();
			std::cout << "(...) - |<xi|yi>|^2 = " << value<< std::endl ;
			std::cout  << "|<xi|yi|^2 = " << xty2.real()<< std::endl ;
////////////////////////////////
			if (value > 0.1)
				condition = 0;



		/*		
			std::cout << "<xi|yi> = " << param.xt * param.y << std::endl;
			std::cout << std::endl << "|sqrt(...) - <xi|yi>|" << std::endl << value << std::endl;
		*/
			
		}
			

		if (condition == 0)
		{
			printf("local minimum. restart\n");
			status = GSL_CONTINUE; 		
			for (int i = 0; i < param.degrees_of_freedom; i++)
			{
				gsl_vector_set (x, i,rand() / (float)RAND_MAX * M_PI * 2);
			}
//			gsl_vector_set_all (ss, 0.0001);
			gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
		}
#endif
	


		}

//#if 0

		if (iter % 4000 == 0)
			printf("iteration: %zu\n f() = %.30f\n size = %.30f\n", iter, s->fval, size);	

//#endif
/*
		if ((iter > 100000 && s->fval > 0.1) || (iter > 200000 && s->fval > 0.01))
			break;
*/

	}
	while (status == GSL_CONTINUE && iter < 200000);

//#if 0

//		print paramterisation params
#if 0

				for (int i = 0; i < param.degrees_of_freedom; i++)
				{
/*					if (i == param.degrees_of_freedom - param.line_length)
					{
						float ret = 1;
						for (int i = 0; i < param.lines - 1; i++)
							ret -= abs(gsl_vector_get (s->x, i * param.line_length)); 
						printf("\n%f, ", ret);
					}
					else if (i % param.line_length == 0)
						printf("\n%f, ", abs(gsl_vector_get(s->x, i)));
*/
					if (i % param.line_length == 0)
						printf("\n");
//					else
						printf("%f, ", gsl_vector_get(s->x, i));
				}
#endif

/*				printf("\n f() = %f size = %f\n", s->fval, size);
				for (int i = 0; i < param.number_of_terms; i++) {
					for (int j = 0; j < param.line_length; j++) {
						arg[i][j] = gsl_vector_get (s->x, i * param.line_length + j);
					}
				}
*/
//						printf("arg[%d,%d] is %f\n", i, j, arg[i][j]);

//		print condition
/*
		std::cout << endl << "t = " << param.t << endl << "d = " << d << endl;
		std::cout << std::endl << "M is" << std::endl << param.M << std::endl;
*/
		MyMatrixXd x1[param.number_of_terms];
		MyMatrixXd y1[param.number_of_terms];
		MyMatrixXd xt1[param.number_of_terms];
		MyMatrixXd yt1[param.number_of_terms];

//		this function gets the vectors solutions from parameters
		for (int i = 0; i < param.number_of_terms; i++)
		{
			init_xy(&param, s->x, i);
			
			x1[i] =  MyMatrixXd(d, 1);
			y1[i] =  MyMatrixXd(d, 1);
			xt1[i] =  MyMatrixXd(d, 1);
			yt1[i] =  MyMatrixXd(d, 1);

			x1[i] = param.x;
			y1[i] = param.y;
			xt1[i] = param.xt;
			yt1[i] = param.yt;

//			std::cout << "x[" << i << "] = " << endl;
//			for (int j = 0; j < param.d; j++)			
			{
		//		fprintf(fptr, "%.21f\t %.30f\t %.10f\n", t, s->fval * 1e24, size);
			}


#if 0

//			print vectors my format // guess forma lower
			std::cout << "x[" << i << "] = " << endl;
			for (int j = 0; j < param.d; j++)			
			{
				std::cout << x1[i](j, 0).real() << " + " << x1[i](j, 0).imag() << " * im" << std::endl;
			}

			std::cout << "y[" << i << "] = " << endl;
			for (int j = 0; j < param.d; j++)			
			{
				std::cout << y1[i](j, 0).real() << " + " << y1[i](j, 0).imag() << " * im" << std::endl;
			}
#endif


			// Yakymenko format
#if 0
			std::cout << "x[" << i << "] = [";
			for (int j = 0; j < param.d; j++)			
			{
				std::cout << x1[i](j, 0).real() << " + " << x1[i](j, 0).imag() << " * im";
				if (j != param.d - 1)
					std::cout << "; ";
			}
			std::cout << "]" << endl;

			std::cout << "y[" << i << "] = [";
			for (int j = 0; j < param.d; j++)			
			{
				std::cout << y1[i](j, 0).real() << " + " << y1[i](j, 0).imag() << " * im";
				if (j != param.d - 1)
					std::cout << "; ";
			}
			std::cout << "]" << endl;


//			std::cout << "norm(x" << i << ") = " << std::endl << sqrt((x1[i] * xt1[i])(0,0)) << std::endl;


		//	std::cout << "y[" << i << "] = " << std::endl << y1[i].real() << " + " << y[i].imag() << " * im" << std::endl << std::endl;

//			std::cout << "norm(y" << i << ") = " << std::endl << sqrt((yt1[i] * y1[i])(0,0)) << std::endl;

//			std::cout << "norm(y" << i << ") = " << std::endl << norm(y1[i]) << std::endl;



			std::cout << std::endl;
			
#endif
		}

#if 0
		printf("bloch sphere parameters\n");
		if (d == 2)
		{
			Eigen::MatrixXd x[param.d2];

		for (int i = 0; i < param.number_of_terms; i++)
		{
			x[i] = Eigen::MatrixXd(3, 1);
			double Tx = 2 * acos(x1[i](0,0).real());
			double Px = std::arg(x1[i](1, 0));

			double Ty = 2 * acos(y1[i](0,0).real());
			double Py = std::arg(y1[i](1, 0));

			printf("Tx%d = %f, Px%d = %f\tTy%d = %f, Py%d = %f\n", i, Tx, i, Px, i, Ty, i, Py);
			printf("check if tetra\n");
		
/*			Eigen::MatrixXd x01;
			x01 = Eigen::matrixxd(3,1);

			Eigen::MatrixXd x01(3,1);
*/
			x[i] = Eigen::MatrixXd(3, 1);
			x[i](0, 0) = sin(Tx) * cos(Px);
			x[i](1, 0) = sin(Tx) * sin(Px);
			x[i](2, 0) = cos(Tx);
		}
		for (int i = 0; i < param.d2; i++)
			for (int j = 0; j < param.d2; j++)
				printf("3D vector product x[%d] * x[%d] = %f\n", i, j, (x[i].adjoint() * x[j])(0, 0));



		Eigen::MatrixXd x01(3,1);
		x01 = x[1] + 0.5 * (x[0] - x[1]);
		Eigen::MatrixXd x23(3, 1);
		x23 = x[3] + 0.5 * (x[2] - x[3]);

		Eigen::MatrixXd x01t(3, 1);
		Eigen::MatrixXd x23t(3, 1);
		
		x01t = x01;
		x01t.adjointInPlace();

		x23t = x23;
		x23t.adjointInPlace();

		Eigen::MatrixXd x02(3,1);
		x02 = x[2] + 0.5 * (x[0] - x[2]);
		Eigen::MatrixXd x13(3, 1);
		x13 = x[3] + 0.5 * (x[1] - x[3]);


/*
		Eigen::MatrixXd S0123 = (x01 - x23);
		Eigen::MatrixXd S0213 = (x02 - x13);
*/
		Eigen::MatrixXd S0123 = 0.5 * (x[1] + x[0] - (x[2] + x[3]));
		Eigen::MatrixXd S0213 = 0.5 * (x[2] + x[0] - (x[1] + x[3]));

		double a = (S0213.adjoint() * (x[0] - x[2]))(0, 0);
		printf("a = %f\n", a);

		double b = (S0213.adjoint() * (x[0] - x[3]))(0, 0);
		printf("b = %f\n", b);


		Eigen::MatrixXd S0213t = S0213;
		S0213t.adjointInPlace();

		double value2 = 0.25 * (2 * x[1].adjoint() * x[2] - 2 * x[0].adjoint() * x[3])(0,0);

		double value3 = 0.25 * (2 * x[0].adjoint() * x[2] - 2 * x[1].adjoint() * x[3])(0,0);

		double value4 = 0.25 * (2 * x[0].adjoint() * x[3] - 2 * x[1].adjoint() * x[2])(0,0);

//		double value = (S0213t * S0123)(0, 0);
//		double value = ((0.5 * (x[1] + x[0] - (x[2] + x[3]))).adjoint() *  (0.5 * (x[2] + x[0] - (x[1] + x[3]))))(0, 0);

		std::cout << "x[1] * x[1]" << (x[1].adjoint() * x[1])(0, 0) << std::endl; 
/*
		double value = 0.25 * ((x[1].adjoint() * x[2])(0, 0) + (x[1].adjoint() *x[0])(0, 0) - (x[1].adjoint() *x[1])(0, 0) - (x[1].adjoint() *x[3])(0, 0))
					+( (x[0].adjoint() *  x[2])(0, 0) +  (x[0].adjoint() *  x[0])(0, 0) -  (x[0].adjoint() *  x[1])(0, 0) -  (x[0].adjoint() *  x[3])(0, 0))
					-  ((x[2].adjoint() * x[2])(0, 0) + (x[2].adjoint() * x[0])(0, 0) - (x[2].adjoint() * x[1])(0, 0) - (x[2].adjoint() * x[3])(0, 0))
					- ( (x[3].adjoint()  * x[2])(0, 0) + ( x[3].adjoint()  *x[0])(0, 0) -  (x[3].adjoint()  *x[1])(0, 0) -  (x[3].adjoint() * x[3])(0, 0));
*/
		double value = 0.25 * (x[1].adjoint() * x[2] + x[1].adjoint() *x[0] - x[1].adjoint() *x[3]
					+ x[0].adjoint() *  x[2] -  x[0].adjoint() *  x[1] -  x[0].adjoint() *  x[3]
					-  ( x[2].adjoint() * x[0] - x[2].adjoint() * x[1] - x[2].adjoint() * x[3])
					- ( x[3].adjoint()  * x[2] +  x[3].adjoint()  * x[0] -  x[3].adjoint()  * x[1]))(0, 0);



//		double value = 0.25 * (2 * (x[1].adjoint() * x[2])(0, 0) - 2 * (x[0].adjoint() * x[3])(0, 0));


		printf("tetraeder axis product: %.10f, %.10f, %.10f, %.10f\n", value, value2, value3, value4);
		}


#endif

#if 0
		printf("next guess format \n");
		for (int i = 0; i < param.number_of_terms; i++)
		{

			std::cout << "guess[" << 2 * i << "] << ";
			for (int j = 0; j < param.d; j++)			
			{
				printf("%.20f + %.20f * im", x1[i](j, 0).real(), x1[i](j, 0).imag());
				if (j != param.d - 1)
					cout << ", ";
			}
			cout << ";" << endl;
			std::cout << "guess[" << 2 * i + 1 << "] << ";
			for (int j = 0; j < param.d; j++)			
			{
//				std::cout << y1[i](j, 0).real() << " + " << y1[i](j, 0).imag() << " * im";

				printf("%.20f + %.20f * im", y1[i](j, 0).real(), y1[i](j, 0).imag());	
				if (j != param.d - 1)
					cout << ", ";
			}
			cout << ";" << endl;
	

//			std::cout << "norm(x" << i << ") = " << std::endl << sqrt((x1[i] * xt1[i])(0,0)) << std::endl;


		//	std::cout << "y[" << i << "] = " << std::endl << y1[i].real() << " + " << y[i].imag() << " * im" << std::endl << std::endl;

//			std::cout << "norm(y" << i << ") = " << std::endl << sqrt((yt1[i] * y1[i])(0,0)) << std::endl;

//			std::cout << "norm(y" << i << ") = " << std::endl << norm(y1[i]) << std::endl;


			std::cout << std::endl;
			

		}
#endif

#if 0
		printf("weyl generated vectors from x0, y0\n");
		for (int i = 0; i < param.number_of_terms; i++)
		{
//			Weyl generated vector
			MyMatrixXd xW = MyMatrixXd(param.d, 1);
			MyMatrixXd yW = MyMatrixXd(param.d, 1);


			int j = i / param.d;
			int k = i % param.d;

	complex<double>X[param.d];
	complex<double>Y[param.d];

	for (int l = 0; l < param.d; l++)
		X[l] = (param.W[j][k] * x1[0])(l, 0);
	for (int l = 0; l < param.d; l++)
		xW(l, 0) = X[l];
	for (int l = 0; l < param.d; l++)
		Y[l] = (param.W[j][k] * y1[0])(l, 0);
	for (int l = 0; l < param.d; l++)
		yW(l, 0) = Y[l];




/*
			complex<double>x_0 = (param.W[j][k] * x1[0])(0,0);
			complex<double>x_1 = (param.W[j][k] * x1[0])(1,0);
			xW(0,0) = x_0;
			xW(1,0) = x_1;

			complex<double>y_0 = (param.W[j][k] * y1[0])(0,0);
			complex<double>y_1 = (param.W[j][k] * y1[0])(1,0);
			yW(1,0) = y_1;
			yW(0,0) = y_0;
*/
			std::cout << "xW[" << i << "] = " << std::endl << xW << std::endl;
//			std::cout << "norm(x" << i << ") = " << std::endl << sqrt((x1[i] * xt1[i])(0,0)) << std::endl;


			std::cout << "yW[" << i << "] = " << std::endl << yW << std::endl << std::endl;


			
		}

#endif
//		print conditions
//#if 0
	int check = 1;
//	std::cout << std::endl << "Condition 3.2 (1):" << std::endl;
	for (int i = 0; i < param.number_of_terms; i++)
		{
		//	init_xy(&param, s->x, i);

/*
			std::complex<double> xy = sqrt((t * (param.d2 - 1) + 1) / (double)param.d);
			double value = norm(xy - (param.xt * param.y)(0,0));
			std::cout << std::endl << "|sqrt(...) - <xi|yi>| = " << value;
*/
			double xy2 = (t * (param.d2 - 1) + 1) / (double)param.d;
			std::complex<double> xty2 = std::norm((xt1[i] * y1[i])(0,0)); 
			double value = xy2 - xty2.real();
/*
			std::complex<double> xty = (xt1[i] * y1[i])(0, 0);
			std::cout << std::endl << "<xi|yi> = r * e^i" << std::arg(xty);
*/
//			get vectors for future check of conditions

//#if 0
			if (i == 0)
				std::cout << std::endl << "Condition 3.2 (1):" << std::endl;

			std::cout << std::endl << "(...) - |<xi|yi>|^2 = " << value;
			std::cout << std::endl << "<xi|yi> = " << xt1[i] * y1[i];
//#endif
			if (value > prec * 10)
				check = 0;
		}

//		std::cout << std::endl;
//		std::cout << std::endl << "Condition 3.2 (2) and 3.3:" << std::endl;

		for (int i = 0; i < param.number_of_terms; i++)
		{
			for (int j = 0; j < param.number_of_terms; j++)
			{
				if (i != j)
				{	
					double xy = (1 - t) / (double)param.d;
					std::complex<double> xty = (xt1[i] * y1[j])(0,0);
//						squared <x|y>
					xty = conj(xty) * xty;
					double sxty = xty.real();
					double value1 = sxty - xy;
					if (value1 > prec * 10)
						check = 0;
					complex<double> xtxyty = ((xt1[j] * x1[i])(0,0)) * ((yt1[i] * y1[j])(0,0));
					xtxyty = xtxyty * conj(xtxyty);
					double sxtxyty = xtxyty.real();
					double value2 = sxtxyty - t * t; 
					if (value2 > prec * 10)
						check = 0;

/*					if (i == 0 && j == 0)
						std::cout << endl << std::endl << "Condition 3.2 (2) and 3.3:" << std::endl;
					std::cout << "|x" << i << " * y" << j << "|^2 - (1-t)/d = " << value1 << std::endl;
					std::cout << "|<x" << j << "|x" << i << "> * <y" << i << "|y" << j << ">| - t =  " << value2 << std::endl;
*/
					if (i < j)
					{
					printf("abs product of x': |<x%d|x%d>| = ", i, j);

					std::cout <<  std::abs((x1[i].adjoint() * x1[j])(0, 0)) << std::endl;
					}


/*
					printf("product of x': <x%d|x%d> = ", i, j);
					std::cout <<  (x1[i].adjoint() * x1[j])(0, 0) << std::endl;
*/

/*
					printf("arg og product of x': arg(<x%d|x%d>) = ", i, j);
					std::cout <<  std::arg((x1[i].adjoint() * x1[j])(0, 0)) << std::endl;
*/
/*					printf("<x2|x0><x0|x1> = ");
					std::cout << (x1[2].adjoint() * x1[0]) * (x1[0].adjoint() * x1[1]) << endl;

					printf("<x2|(x0><x0)|x1> = ");
					std::cout << x1[2].adjoint() * (x1[0] * x1[0].adjoint()) * x1[1] << endl;

					printf("x0><x0|> = ");
					std::cout << x1[0] * x1[0].adjoint() << endl;


					printf("<x2|x1> = ");
					std::cout << x1[2].adjoint() * x1[1] << endl;

					printf("x20|x10 = ");
					std::cout << (x1[2](0, 0) * x1[1](0, 0)) << endl;
*/


				}
			}
		}
//#endif
//		printf("f() = %.30f\n size = %.30f\n", (double)s->fval, (double)size);
//		fprintf(fptr, "%.21f\t %.30f\t %.10f\n", t, s->fval * 1e24, size);
//		printf("%.10f\t %.30f\t %.20f\n", t, s->fval, size);

//		if (s->fval < prec * 10 && check == 1)
		{
			fprintf(fptr, "%f\t %f\t %.30f\t %.20f\n", param.t, x1[0](0, 0).real(), s->fval, size);
			float x, y, z;
			x = (x1[0].adjoint() * (param.W[0][1] * x1[0]))(0, 0).real();
			y = (x1[0].adjoint() * (param.W[1][0] * x1[0]))(0, 0).real();
			z = (x1[0].adjoint() * (param.W[1][1] * x1[0]))(0, 0).real();
//			std::cout << "z = " << (x1[0].adjoint() * (param.W[1][1] * x1[0])) << endl;
			printf("%f\t  %f\t  %f\t %f\t %.30f\t %.20f\n", param.t, x1[0](0, 0).real(),  x1[0](1, 0).real(), x1[0](2, 0).real(), s->fval, size);

//			printf("%f\t %f\t %f\t %f\t %.30f\t %.20f\n", param.t, x, y, z, s->fval, size);


		}
		t = t - eps;
		param.t = t;
		gsl_vector_set_all (ss, step_size);
		gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
	}

///////////////////////////////////////////////////

//	printf("CHECK t = %.30f, f = %.30f\n", t, my_f(param.v, &param));


/*
	for (int i = 0; i < 1000; i++)
	{	
		t = t - 0.001;
		param.t = t;
//                       t   d
		printf("%.30f\t %.30f\n", t, my_f(param.v, &param));

	}
*/
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return status;
}

