
#include </usr/local/Cellar/gsl/2.7.1/include/gsl/gsl_multimin.h>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <Eigen/Dense>


//typedef Eigen::Matrix<std::complex<double>, 2, 2> ComplexMatrix2d;
//typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> Eigen::MatrixXcd;



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
			MyMatrixXd result(colb,rowb);

			for (int i = 0; i < colb; i++)
			{
				for (int j = 0; j < rowb; j++)
				{
					result(i,j) = M2(i,j);
				}
			}			
			return (result);	
		}
		virtual ~MyMatrixXd()
		{
//			std::cout<<"\n Destructor executed";
		}
};

template <typename T> double sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

using namespace std;

/*function to minimize*/
	//v = the argument of the fucntion
double	my_f (const gsl_vector *v, void *params)
{
	double *p = (double *)params;
	double t = p[0];
	double d = p[1];

	double arg[4][5];
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 5; j++) {
			arg[i][j] = gsl_vector_get (v, i * 5 + j);
	//		printf("arg[%d,%d] is %f\n", i, j, arg[i][j]);
		}
	}
	arg[3][0] = 1 - abs(arg[0][0]) -  abs(arg[1][0]) -  abs(arg[2][0]);

/*	std::complex<double> Mx[4][2][2];
	std::complex<double> My[4][2][2];
*/
	std::complex<double> fullx[4];
	std::complex<double> fully[4];

	MyMatrixXd Mx[] = {MyMatrixXd(2,2), MyMatrixXd(2,2), MyMatrixXd(2,2), MyMatrixXd(2,2)};
	MyMatrixXd My[] = {MyMatrixXd(2,2), MyMatrixXd(2,2), MyMatrixXd(2,2), MyMatrixXd(2,2)}; 


	double B, A, Y, X, w;
	for (int i = 0; i < 4; i++)
	{
//		w_i * |x><x|_i
		w = abs(arg[i][0]);
		Y = arg[i][1];
		X = arg[i][2];
		B = arg[i][3];
		A = arg[i][4];

//		std::complex<double> x[2];
//		std::complex<double> y[2];

		MyMatrixXd x(2,1);
		MyMatrixXd y(2,1);

		MyMatrixXd xt(2,1);
		MyMatrixXd yt(2,1);

		x(0, 0) = cos(Y);
		x(1, 0) = sgn(sin(Y)) * polar(abs(sin(Y)), X);
		y(0, 0) = cos(B); 
		y(1, 0) = sgn(sin(B)) * polar(abs(sin(B)), A);

		xt = x;
		yt = y;

//		xt.transposeInPlace();
//		xt.conjugate();

		xt.adjointInPlace() ;

//		std::cout << "x is: \n" << x << std::endl;
//		std::cout << "transpose conjugate is \n" << xt << std::endl;

//		yt.transposeInPlace();
//		yt.conjugate();
		yt.adjointInPlace() ;

		x.Kron(xt, Mx[i]);	
		y.Kron(yt, My[i]);

		Mx[i] *= w;


/*
		Mx[i](0,0) = w * Mx[i](0,0);
		Mx[i](1,0) = w * Mx[i](1,0);
		Mx[i](0,1) = w * Mx[i](0,1);
		Mx[i](1,1) = w * Mx[i](1,1);
*/

		std::cout << "Mx is " << std::endl << Mx[i] << std::endl;

//		fullx[i] = {x[0],x[1]};
//		fully[i] = {y[0], y[1]};



//#if 0

/*		Mx[i](0, 0) = w * x(0, 0) * x(0, 0);
		Mx[i](0, 1) = w * x(0, 0) * conj(x(1, 0));
		Mx[i](1, 0) = w * x(1, 0) * x(0, 0);
		Mx[i](1, 1) = w * x(1, 0) * conj(x(1, 0));
*/

		std::cout << " now Mx is " << std::endl;

		std::cout << Mx[i](0, 0) << Mx[i](0, 1) << endl;
		std::cout << Mx[i](1, 0) << Mx[i](1, 1) << endl;



//		|y><y|_i
		My[i](0, 0) = y(0, 0) * y(0, 0);
		My[i](0, 1) = y(0, 0) * conj(y(1, 0));
		My[i](1, 0) = y(1, 0) * y(0, 0);
		My[i](1, 1) = y(1, 0) * conj(y(1, 0));

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
	}
/*
	printf("|<x0| y1>|^2 is %f\n", norm(conj(fullx[0]) * fully[1]));
	printf("(1 - t) / d = %f\n", (1 - t)/ d);
*/
//	std::cout << endl;

//	std::complex<double> MxMy[4][4][4];
	MyMatrixXd MxMy[] = {MyMatrixXd(4, 4), MyMatrixXd(4, 4), MyMatrixXd(4, 4), MyMatrixXd(4, 4)};
	for (int i = 0; i < 4; i++)
	{
	/*
		xx00
		xx00
		0000
		0000
	*/

		Mx[i].Kron(My[i], MxMy[i]);
#if 0

		MxMy[i][0][0] = Mx[i][0][0] * My[i][0][0];
		MxMy[i][0][1] = Mx[i][0][0] * My[i][0][1];
		MxMy[i][1][0] = Mx[i][0][0] * My[i][1][0];
		MxMy[i][1][1] = Mx[i][0][0] * My[i][1][1];

	/*
		00xx
		00xx
		0000
		0000
	*/

		MxMy[i][2][0] = Mx[i][1][0] * My[i][0][0];
		MxMy[i][2][1] = Mx[i][1][0] * My[i][0][1];
		MxMy[i][3][0] = Mx[i][1][0] * My[i][1][0];
		MxMy[i][3][1] = Mx[i][1][0] * My[i][1][1];


	/*
		0000
		0000
		xx00
		xx00
	*/

		MxMy[i][0][2] = Mx[i][0][1] * My[i][0][0];
		MxMy[i][0][3] = Mx[i][0][1] * My[i][0][1];
		MxMy[i][1][2] = Mx[i][0][1] * My[i][1][0];
		MxMy[i][1][3] = Mx[i][0][1] * My[i][1][1];

	/*
		0000
		0000
		00xx
		00xx
	*/

		MxMy[i][2][2] = Mx[i][1][1] * My[i][0][0];
		MxMy[i][2][3] = Mx[i][1][1] * My[i][0][1];
		MxMy[i][3][2] = Mx[i][1][1] * My[i][1][0];
		MxMy[i][3][3] = Mx[i][1][1] * My[i][1][1];

#endif
//	printf("MxMy[0] is");
//	std::cout <<  MxMy[i][0][0] << endl << endl;
	}


	Eigen::MatrixXcd M(4,4);
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			M(i, j) = 0;
		}
	}

/*
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{*/
			for (int k = 0; k < 4; k++)
			{
				M += MxMy[k];
			}
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
	Eigen::Matrix<double,4,4> Usw;

//	first coordinate is y
	
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


//	for (int j = 0; j < 4; j++)
///	{
//		for (int k = 0; k < 4; k++)
//		{
			M -= (t / (double)d) * Usw;
//		}
//	}

//	printf(" t / d is %f\n", (t / (double)d));
//	std::cout << "M[0][0] - t / d is " << M[0][0] << endl;


//	for (int k = 0; k < 4; k++)
//	{
		M-= (1-t) / (double)(d*d) *  Eigen::Matrix4d :: Identity ();
//	}
//	std::cout << "Matrix is \n" << M <<std::endl;	
	double n = 0;

		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
	//			printf("M[%d][%d] = ", j, k);
	//			std::cout << M[j][k] << endl;
				n += norm(M(j,k));
			}
		}
	return (n);
}


int	main(void)
{
	srand(time(0));
	double par[2] = {1/9.0, 2};

	const gsl_multimin_fminimizer_type *T =
		gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;

	size_t iter = 0;
	int status;
	double size;
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
	x = gsl_vector_alloc (20);
//	gsl_vector_set_all (x, 0);




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



	for (int i = 0; i < 20; i++)
		gsl_vector_set (x, i,rand() / (float)RAND_MAX * M_PI * 2);






	printf("vector set\n");
//	gsl_vector_set_all (x, 0);
	



	/* Set initial step sizes to 1 */
	ss = gsl_vector_alloc (20);
	gsl_vector_set_all (ss, 0.1);

	/* Initialize method and iterate */
	minex_func.n = 20;
	minex_func.f = my_f;
	minex_func.params = par;

	s = gsl_multimin_fminimizer_alloc (T, 20);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-10);

		if (status == GSL_SUCCESS)
		{
			printf ("converged to minimum at\n");

		printf ("%d:\n%f,%f,%f,%f,%f\n%f,%f,%f,%f,%f\n%f,%f,%f,%f,%f\n%f,%f,%f,%f,%f\n f() = %f size = %f\n",
						iter,
						abs(gsl_vector_get (s->x, 0)),
						gsl_vector_get (s->x, 1),
						gsl_vector_get (s->x, 2),
						gsl_vector_get (s->x, 3),
						gsl_vector_get (s->x, 4),

						abs(gsl_vector_get (s->x, 5)),
						gsl_vector_get (s->x, 6),
						gsl_vector_get (s->x, 7),
						gsl_vector_get (s->x, 8),
						gsl_vector_get (s->x, 9),

						abs(gsl_vector_get (s->x, 10)),
						gsl_vector_get (s->x, 11),
						gsl_vector_get (s->x, 12),
						gsl_vector_get (s->x, 13),
						gsl_vector_get (s->x, 14),

						1 - abs(gsl_vector_get (s->x, 0)) - abs(gsl_vector_get (s->x, 5))  - abs(gsl_vector_get (s->x, 10)),
						gsl_vector_get (s->x, 16),
						gsl_vector_get (s->x, 17),
						gsl_vector_get (s->x, 18),
						gsl_vector_get (s->x, 19),

						s->fval, size);

		}
		}
	while (status == GSL_CONTINUE && iter < 100000);

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return status;
}


