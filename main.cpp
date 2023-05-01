
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

	complex<double>x0 = (param->W[j][k] * param->x)(0,0);
	complex<double>x1 = (param->W[j][k] * param->x)(1,0);
	param->x(0,0) = x0;
	param->x(1,0) = x1;

	complex<double>y0 = (param->W[j][k] * param->y)(0,0);
	complex<double>y1 = (param->W[j][k] * param->y)(1,0);
	param->y(1,0) = y1;
	param->y(0,0) = y0;




//	std::cout << "x is " << std::endl << param->x << std::endl;


	}

	else if (param->d == 3)
	{
		
		param->x(0,0) = sqrt(2/3.0);
		param->x(1, 0) = -1 / sqrt(6);
		param->x(2, 0) = -1/sqrt(6);

		param->y(0,0) = sqrt(2/3.0);
		param->y(1, 0) = -1 / sqrt(6);
		param->y(2, 0) = -1/sqrt(6);


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

	}

	param->xt = param->x;
	param->yt = param->y;
	param->xt.adjointInPlace() ;
	param->yt.adjointInPlace();



}

 

//	inits i'th vector
void		init_xy(t_param *param, const gsl_vector *x, int i)
{

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

		param->x(0, 0) = cos(Y); //					cos(Y)
//					magnitude (>0 always)
//		x(1, 0) = sgn(sin(Y)) * polar(abs(sin(Y)), X);//	sin(Y) * e^iX
		param->x(1, 0) = sin(Y) * polar(1.0, X);//	sin(Y) * e^iX
		param->y(0, 0) = cos(B); //					cos(B)
//		y(1, 0) = sgn(sin(B)) * polar(abs(sin(B)), A);//	sin(B) * e^iA
		param->y(1, 0) = sin(B) * polar(1.0, A);//	sin(B) * e^iA

//		std::cout << "vector is " << param->x << std::endl;
/*
////		EXPERIMENT

		double B, A, Y, X, w;

		Y = gsl_vector_get(x, 0 + param->line_length * 0);
		X = gsl_vector_get(x, 1 + param->line_length * 0);
		B = gsl_vector_get(x, 2 + param->line_length * 0);
		A = gsl_vector_get(x, 3 + param->line_length * 0);

		param->x(0, 0) = cos(Y);
		param->x(1, 0) = sin(Y) * polar(1.0, X);
		param->y(0, 0) = cos(B);
		param->y(1, 0) = sin(B) * polar(1.0, A);


	complex<double>x0 = (param->W[j][k] * param->x)(0,0);
	complex<double>x1 = (param->W[j][k] * param->x)(1,0);
	param->x(0,0) = x0;
	param->x(1,0) = x1;

	complex<double>y0 = (param->W[j][k] * param->y)(0,0);
	complex<double>y1 = (param->W[j][k] * param->y)(1,0);
	param->y(1,0) = y1;
	param->y(0,0) = y0;
*/

		
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
	if (p->d == 2)
	{

		double *ret = (double *)malloc(sizeof(double) * 2 * 2);
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

		return (ret);
	}
	if (p->d == 4)
	{
		double *ret = (double *)malloc(sizeof(double) * 6 * 2);
	double sinxga = std::abs(x(3, 0));
	double sinxga2 = std::norm(x(3, 0));
	double xga = asin(sinxga);
	double Xc = std::arg(x(3, 0));

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
	ret[5] = Xc;

//	now y
	ret[6] = yal;
	ret[7] = ybe;
	ret[8] = yga;
	ret[9] = yA;
	ret[10] = yB;
	ret[11] = Xc;

	return (ret);
	}
}


//		this function takes all pairs of vectors written to v and extracts angles.
//		write ith pair here 	//guess in vector form	//number of the pair of vectors
void	init_minimizing_guess(t_param *p, Eigen::MatrixXcd *v)
{
	int ll = p->line_length;
	double *vector_params;
//	for (int j = 0; j < ll; j++)
//		vector_params[j] = get_angles_form_

//	line_length  = 2 * degrees
//	ith x vector
	int number_of_lines = p->number_of_terms;
	for (int i = 0; i <  number_of_lines; i++)
	{
	vector_params = get_angles_from_vector(p, v[2 * i], v[2 * i + 1]);
	for (int j = 0; j < ll; j++)
	{
//								ith line
		printf("setting param %f\n", vector_params[j]);
		gsl_vector_set(p->v, j + ll * i, vector_params[j]);
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
	if (p->d == 2 || p->d == 4)
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



/*	for (int i = 0; i < d2; i++)
	{
		for (int j = 0; j < d2; j++)
		{
			M(i, j) = 0;
		}
	}
*/

	p->M = p->Zero;

//	std::cout << p->M << std::endl;
//	std::cout << p->Identity << std::endl;

	for (int i = 0; i < number_of_terms; i++)
	{
//		|x><x|_i
//		w = abs(arg[i][0]);
		w = 1 / (double)p->number_of_terms;

//		printf("init xy\n");

//		here only 1 pair of vectors is calculated from d^2 pairs.
		init_xy(p, v, i);

//		init_with_known_solution(p, v, i);
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
                        p->W[i][j] = pow(C, i) * pow(S, j);

			cout << "W" << i << j << endl;
                        cout << p->W[i][j] << endl << endl;

		}
	}
//#endif

}








int	main(void)
{
	srand(time(0));
	int d = 2;
	double t = 0.5 * 1/( (double) (d + 1));
//                       t   d
	double par[2] = {t, (double)d};
	printf("d = %d, t = %f\n", d, t);

	const gsl_multimin_fminimizer_type *T =
		gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;

	size_t iter = 0;
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



//		this function takes all pairs of vectors and extracts angles.
//		i ranges from 0 to 2 * number_of_terms
//		write ith pair here 	//guess in vector form	//number of the pair of vectors


	Eigen::MatrixXcd guess[param.number_of_terms * 2];
	for (int i = 0; i < param.number_of_terms * 2; i++)
		guess[i] = Eigen::MatrixXcd(param.d, 1);


/*
	guess[0] << 0.9851137757355578, 0.12155420366246719 + 0.12155420366246718 * im; //x0
	guess[1] << 0.7671817537135385, 0.4535593438399321 + 0.45355934383993207 * im; //y0

	guess[2] << 0.9851137757355578, -0.12155420366246719 - 0.12155420366246718 * im; //x1
	guess[3] << 0.7671817537135385 + 0.0* im, -0.4535593438399321 - 0.45355934383993207* im; //y1

	guess[4] << 0.171938 + 0.0 * im, 0.696576 + -0.696576 * im; //x2
	guess[5] << 0.641381 + 0.0 * im, 0.542509 + -0.542509 * im; //y2

	guess[6] << 0.171938 + 0.0 * im, -0.696576 + 0.696576 * im; //x3
	guess[7] << 0.641381 + 0.0 * im, -0.542509 + 0.542509 * im; //y3



	printf("guess inited\n");


	init_minimizing_guess(&param, guess);
*/


/*
	for (int i = 0; i < param.degrees_of_freedom; i++)
		gsl_vector_set (x, i,rand() / (float)RAND_MAX * M_PI * 2);
*/


	
	printf("vector set\n");
//	gsl_vector_set_all (x, 0);
	



	/* Set initial step sizes to 1 */
	ss = gsl_vector_alloc (param.degrees_of_freedom);
	gsl_vector_set_all (ss, 0.001);

	/* Initialize method and iterate */
	minex_func.n = param.degrees_of_freedom;
	minex_func.f = my_f;
	minex_func.params = &param;


	s = gsl_multimin_fminimizer_alloc (T, param.degrees_of_freedom);
	printf("gsl multimin fminimizer set\n");
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	double arg[param.number_of_terms][param.line_length];

	printf("start minimizing\n");

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-7);

		if (status == GSL_SUCCESS)
		{
				printf ("converged to minimum at %zu:\n", iter);

//#if 0
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
			std::complex<double> xy = sqrt((t * (param.d2 - 1) + 1) / (double)param.d);
			double value = norm(xy - (param.xt * param.y)(0,0));
	*/
			std::complex<double> xy2 = (t * (param.d2 - 1) + 1) / (double)param.d;
			double value = norm(xy2 - (param.xt * param.y)(0,0).real() * (param.xt * param.y)(0,0) );
			std::cout << std::endl << "(...) - |<xi|yi>|^2 = " << value;


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
		if (iter % 2000 == 0)
			printf("iteration: %zu\n f() = %f\n size = %f\n", iter, s->fval, size);	

/*		if ((iter > 100000 && s->fval > 0.1) || (iter > 200000 && s->fval > 0.01))
			break;
*/
	}
	while (status == GSL_CONTINUE && iter < 20000);

//#if 0

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

/*				printf("\n f() = %f size = %f\n", s->fval, size);
				for (int i = 0; i < param.number_of_terms; i++) {
					for (int j = 0; j < param.line_length; j++) {
						arg[i][j] = gsl_vector_get (s->x, i * param.line_length + j);
					}
				}
*/
//						printf("arg[%d,%d] is %f\n", i, j, arg[i][j]);

//		print condition
		std::cout << endl << "t = " << param.t << endl << "d = " << d << endl;
		std::cout << std::endl << "M is" << std::endl << param.M << std::endl;

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

			std::cout << "x[" << i << "] = " << std::endl;
			for (int j = 0; j < param.d; j++)			
			{
				std::cout << x1[i](j, 0).real() << " + " << x1[i](j, 0).imag() << " * im" << std::endl;
			}

			std::cout << "x[" << i << "] = " << std::endl;
			for (int j = 0; j < param.d; j++)			
			{
				std::cout << y1[i](j, 0).real() << " + " << y1[i](j, 0).imag() << " * im" << std::endl;
			}


//			std::cout << "norm(x" << i << ") = " << std::endl << sqrt((x1[i] * xt1[i])(0,0)) << std::endl;


		//	std::cout << "y[" << i << "] = " << std::endl << y1[i].real() << " + " << y[i].imag() << " * im" << std::endl << std::endl;

//			std::cout << "norm(y" << i << ") = " << std::endl << sqrt((yt1[i] * y1[i])(0,0)) << std::endl;

//			std::cout << "norm(y" << i << ") = " << std::endl << norm(y1[i]) << std::endl;


			std::cout << std::endl;
			

		}
		printf("next guess form\n");
		for (int i = 0; i < param.number_of_terms; i++)
		{

			std::cout << "guess[" << 2 * i << "] = " << std::endl;
			for (int j = 0; j < param.d; j++)			
			{
				std::cout << x1[i](j, 0).real() << " + " << x1[i](j, 0).imag() << " * im";
				if (j != param.d - 1)
					cout << ", ";
			}
			cout << ";" << endl;
			std::cout << "guess[" << 2 * i + 1 << "] = " << std::endl;
			for (int j = 0; j < param.d; j++)			
			{
				std::cout << y1[i](j, 0).real() << " + " << y1[i](j, 0).imag() << " * im" << ", ";
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
/*
		printf("weyl generated vectors from x0, y0\n");
		for (int i = 0; i < param.number_of_terms; i++)
		{
//			Weyl generated vector
			MyMatrixXd xW = MyMatrixXd(param.d, 1);
			MyMatrixXd yW = MyMatrixXd(param.d, 1);


			int j = i / param.d;
			int k = i % param.d;

			complex<double>x_0 = (param.W[j][k] * x1[0])(0,0);
			complex<double>x_1 = (param.W[j][k] * x1[0])(1,0);
			xW(0,0) = x_0;
			xW(1,0) = x_1;

			complex<double>y_0 = (param.W[j][k] * y1[0])(0,0);
			complex<double>y_1 = (param.W[j][k] * y1[0])(1,0);
			yW(1,0) = y_1;
			yW(0,0) = y_0;

			std::cout << "xW[" << i << "] = " << std::endl << xW << std::endl;
//			std::cout << "norm(x" << i << ") = " << std::endl << sqrt((x1[i] * xt1[i])(0,0)) << std::endl;


			std::cout << "yW[" << i << " ]= " << std::endl << yW << std::endl << std::endl;


			
		}
*/
	std::cout << std::endl << "Condition 3.2 (1):" << std::endl;
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
			std::cout << std::endl << "(...) - |<xi|yi>|^2 = " << value;
//			std::cout << std::endl << "|<xi|yi|^2 = " << xty2.real();
//			get vectors for future check of conditions
		}

		std::cout << std::endl;
		std::cout << std::endl << "Condition 3.2 (2) and 3.3:" << std::endl;

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
//					std::cout << "xti is " << xty << std::endl;
					double sxty = xty.real();
//					std::cout << "(1 - t) / d = " << xy << std::endl;
					std::cout << "|x" << i << " * y" << j << "|^2 - (1-t)/d = " << sxty - xy << std::endl;
					complex<double> xtxyty = ((xt1[j] * x1[i])(0,0)) * ((yt1[i] * y1[j])(0,0));
					xtxyty = xtxyty * conj(xtxyty);
					double sxtxyty = xtxyty.real();
					std::cout << "|<x" << j << "|x" << i << "> * <y" << i << "|y" << j << ">| - t =  " << sxtxyty - t * t << std::endl;
				}
			}
		}
		printf("f() = %f\n size = %f\n", s->fval, size);


///////////////////////////////////////////////////




	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return status;
}



