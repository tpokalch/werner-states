
#include </usr/local/Cellar/gsl/2.7.1/include/gsl/gsl_multimin.h>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <Eigen/Dense>

using namespace std;



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
			std::cout << "assigned matrix is" << std::endl << M2;
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
}	t_param;


void		init_xy(t_param *param, double *arg, const gsl_vector *x, int i)
{


//#if 0

		int DELETE_ONE_TEST;

		double B, A, Y, X, w;
/*
		Y = arg[0 + param->line_length * i];
		X = arg[1 + param->line_length * i];
		B = arg[2 + param->line_length * i];
		A = arg[3 + param->line_length * i];
*/
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

//#endif

#if 0
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

# endif

		param->xt = param->x;
		param->yt = param->y;
		param->xt.adjointInPlace() ;
		param->yt.adjointInPlace();

//		std::cout << sqrt((t * (param->d2 - 1) + 1) / param->d) << std::endl;
	}
//#endif

////////////////////////////////////////////////////




template <typename T> double sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
/*function to minimize*/
	//v = the argument of the fucntion
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

	p->M = /*Eigen::MatrixXcd::*/p->Zero;

//	std::cout << p->M << std::endl;
//	std::cout << p->Identity << std::endl;

	for (int i = 0; i < number_of_terms; i++)
	{
//		|x><x|_i
//		w = abs(arg[i][0]);
		w = 1 / (double)p->number_of_terms;

//		printf("init xy\n");
		init_xy(p, arg[0], v, i);
//		printf("done\n");
#if 0

		int DELETE_ONE_TEST;


		Y = arg[i][0];
		X = arg[i][1];
		B = arg[i][2];
		A = arg[i][3];

//		std::complex<double> x[2];
//		std::complex<double> y[2];

//		MyMatrixXd x(d,1);
//		MyMatrixXd y(d,1);

//		for future transposition
//		MyMatrixXd xt(d,1);
//		MyMatrixXd yt(d,1);


//		init |x> and |y> from parameters A, B, X, Y, alway explicitly
//		passed to my_f() so that <x|x> = <y|y> = 1

		p->x(0, 0) = cos(Y); //					cos(Y)
//					magnitude (>0 always)
//		x(1, 0) = sgn(sin(Y)) * polar(abs(sin(Y)), X);//	sin(Y) * e^iX
		p->x(1, 0) = sin(Y) * polar(1.0, X);//	sin(Y) * e^iX
		p->y(0, 0) = cos(B); //					cos(B)
//		y(1, 0) = sgn(sin(B)) * polar(abs(sin(B)), A);//	sin(B) * e^iA
		p->y(1, 0) = sin(B) * polar(1.0, A);//	sin(B) * e^iA

#endif

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

		p->Mx *= w;

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
			p->M -= (t / (double)d) * p->Usw;
//		}
//	}

//	printf(" t / d is %f\n", (t / (double)d));
//	std::cout << "M[0][0] - t / d is " << M[0][0] << endl;


//	for (int k = 0; k < 4; k++)
//	{
//		std::cout << "M was" << std::endl << p->M << std::endl;

//			std::cout << (1 - t) / (double)(d * d) * p->Identity << std::endl << std::endl;
		p->M-= (1-t) / (double)(d*d) * /* MyMatrixXd::*/p->Identity;

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











int	main(void)
{
	srand(time(0));
	int d = 2;
	double t = 1/( (double)(d + 1));
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
//	structure of x: weight1, x1_param1, x1_param2, ..., y1_param1, y1_param2, ...
//			weight2, x2_param1, x2_param2, ..., y2_param1, y2_param2, ...
//			line length = 2 * d - 2
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

	for (int j = 0; j < param.d2; j++)
	{
		for (int i = 0; i < param.d2; i++)
		{
			param.Usw(j, i) = MyMatrixXd::Zero(param.d2, param.d2)(j, i);
		}
	}

	MyMatrixXd Basis(param.d, param.d);

/*	std::cout << "Identity:" << std::endl;
	std::cout << Eigen::MatrixXd::Identity(d, d);
	std::cout << std::endl;
*/

//	Basis = MyMatrixXd::Identity(d, d);

	for (int j = 0; j < d; j++)
	{
		for (int i = 0; i < d; i++)
		{
			Basis(j, i) = MyMatrixXd::Identity(d, d)(j, i);
		}
	}

/*	std::cout << std::endl;
	std::cout << "Basis:" << std::endl;
	std::cout << Basis << std::endl;

*/
//#if 0

	for (int j = 0; j < d; j++)
	{
		for (int i = 0; i < d; i++)
		{
//			std::cout << "i j =" << std::endl << i << " " << j << std::endl << std::endl;

			MyMatrixXd vi(d,1);
			MyMatrixXd vj(d,1);
			MyMatrixXd vit(d,1);
			MyMatrixXd vjt(d,1);

			for (int k = 0; k < d; k++)
				vi(k) = Basis.col(i)(k);

			for (int k = 0; k < d; k++)
				vj(k) = Basis.col(j)(k);


//			std::cout << "vj" << std::endl << vj << std::endl;

			vit = vi;
			vjt = vj;
			vit.adjointInPlace();
			vjt.adjointInPlace();

			MyMatrixXd vjvi(param.d2, 1);
			vj.Kron(vi, vjvi);

//			std::cout << "vjvi" << std::endl << vjvi << std::endl;
			
			MyMatrixXd vivjt(param.d2, 1);

			vi.Kron(vj, vivjt);
			
			vivjt.adjointInPlace();
			

//			std::cout << "vivjt" << std::endl << vivjt << std::endl;
			

			MyMatrixXd term(param.d2, param.d2);
			vjvi.Kron(vivjt, term);

			param.Usw = param.Usw + term;
		}
	}
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

	gsl_vector_set (x, 0, 0.17278);
	gsl_vector_set (x, 1, M_PI / 4.0);
	gsl_vector_set (x, 2, 0.6936);
	gsl_vector_set (x, 3, M_PI / 4.0);

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





	for (int i = 0; i < param.degrees_of_freedom; i++)
		gsl_vector_set (x, i,rand() / (float)RAND_MAX * M_PI * 2);

	printf("vector set\n");
//	gsl_vector_set_all (x, 0);
	



	/* Set initial step sizes to 1 */
	ss = gsl_vector_alloc (param.degrees_of_freedom);
	gsl_vector_set_all (ss, 0.01);

	/* Initialize method and iterate */
	minex_func.n = param.degrees_of_freedom;
	minex_func.f = my_f;
	minex_func.params = &param;

	s = gsl_multimin_fminimizer_alloc (T, param.degrees_of_freedom);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	double arg[param.number_of_terms][param.line_length];

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-8);

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
		for (int i = 0; i < param.number_of_terms; i++)
		{
			init_xy(&param, arg[0], s->x, i);

			std::complex<double> xy = sqrt((t * (param.d2 - 1) + 1) / (double)param.d);
			double value = norm(xy - (param.xt * param.y)(0,0));
			if (value > 0.1)
				condition = 0;
			std::cout << std::endl << "|sqrt(...) - <xi|yi>|" << std::endl << value << std::endl;
		
			
		}
			

		if (condition == 0)
		{
			printf("local minimum. restart\n");
			status = GSL_CONTINUE; 		
			for (int i = 0; i < param.degrees_of_freedom; i++)
			{
				gsl_vector_set (x, i,rand() / (float)RAND_MAX * M_PI * 2);
			}
			gsl_vector_set_all (ss, 0.01);
			gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
		}

	


		}
		if (iter % 1000 == 0)
			printf("iteration: %zu\n f() = %f\n size = %f\n", iter, s->fval, size);	

/*		if ((iter > 100000 && s->fval > 0.1) || (iter > 200000 && s->fval > 0.01))
			break;
*/
	}
	while (status == GSL_CONTINUE && iter < 100000);

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
		
		std::cout << std::endl << "M is" << std::endl << param.M << std::endl;

		for (int i = 0; i < param.number_of_terms; i++)
		{
			init_xy(&param, arg[0], s->x, i);

			std::complex<double> xy = sqrt((t * (param.d2 - 1) + 1) / (double)param.d);
			double value = norm(xy - (param.xt * param.y)(0,0));
			std::cout << std::endl << "|sqrt(...) - <xi|yi>| = " << value;
		}
		std::cout << std::endl;


///////////////////////////////////////////////////




	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return status;
}



