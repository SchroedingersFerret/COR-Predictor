/*
 * COR-cor.hpp
 *
 *  Copyright 2019
 *      J. Ball (SchroedingersFerret)
 */
 
//This file is part of COR-Predictor.
//
//   COR-Predictor is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   COR-Predictor is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with COR-Predictor.  If not, see <https://www.gnu.org/licenses/>.
    
#ifndef COR_COR_HPP_
#define COR_COR_HPP_

//generates a random double between -1.0 and 1.0
double cor::RandInit()
{
	double r = rand() % RAND_MAX + (RAND_MAX/2 - 1);
	r /= RAND_MAX;
	return r;
};

//generates a random double between 0 and 1.0
double cor::rand_double()
{
	return ((double) rand() / (RAND_MAX));
};

//generates a random bool
bool cor::rand_bool()
{
	int rnd = rand() % (RAND_MAX-1);
	rnd /= (RAND_MAX-1)/2;
	if (rnd==1)
		return true;
	if (rnd==0)
		return false;
	else
		return rand_bool();
}

//Evaluates Chebyshev approximation at x with coefficients from param[]
double cor::Chebyshev(double x, std::vector<double> param)
{
	int ni = param.size();
	double b1 = 0.f, b2 = 0.f;
	for (int i=ni-1; i>0; --i)
	{
		double temp = b1;
		b1 = 2.f*x*b1-b2+param[i];
		b2 = temp;
	}
	return x*b1-b2+param[0];
}

//combines material properties
double cor::combine(double x, double y)
{
	return sqrt(0.5*(x*x + y*y));
}

//returns the approximate COR with independent variables x[] and coefficients parameters[][]
double cor::f(std::vector<double> x, parameters param)
{
	double y1 = Chebyshev(x[0],param.c[0]);
	y1 /= Chebyshev(x[2],param.c[1]);
	y1 /= Chebyshev(x[4],param.c[2]);
	double y2 = Chebyshev(x[1],param.c[0]);
	y2 /= Chebyshev(x[3],param.c[1]);
	y2 /= Chebyshev(x[5],param.c[2]);
	return 3.1*combine(y1,y2)/Chebyshev(x[6],param.c[3]);
}

//gets the residual of each datapoint
std::vector<double> cor::GetResiduals(std::vector<double> y, std::vector<std::vector<double> > x, parameters param)
{
	std::vector<double> residuals(n_data);
	for (int i=0; i<n_data; ++i)
	{
		double yi = f(x[i],param);
		residuals[i] = y[i]-yi;
	}
	return residuals;
}

//returns the sum of the square of each residual
double cor::SumOfSquares(std::vector<double> residuals)
{
	double sum = 0;
	int ni = residuals.size();
	for (int i=0; i<ni; ++i)
	{
		sum += residuals[i]*residuals[i];
	}
	return sum;
}

#endif /* COR_COR_HPP_ */
