/*
 * COR-optimization.hpp
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
    
#ifndef COR_OPTIMIZATION_HPP_
#define COR_OPTIMIZATION_HPP_

class optimization : public COR_predictor
{
	public:	
		static bool rand_bool();
		static double RandInit();
		static double combine(double x, double y);
		static double Chebyshev(double x, std::vector<double> &param, double b, double a);
		static double f(std::vector<double> &x, std::vector<std::vector<double> > &param);	
		static std::vector<double> GetResiduals(std::vector<double> &y, std::vector<std::vector<double> > &x, std::vector<std::vector<double> > &param);
		static double Mean_square_error(std::vector<double> residuals);
};

class genetic : public optimization
{
	private:
		static std::vector<std::bitset<384> > encode(std::vector<std::vector<double> > &param);
		static std::vector<std::vector<double> > decode(std::vector<std::bitset<384> > &w);
		static int partition(std::vector<double> &cost, std::vector<int> &index, int low, int high);
		static void quicksort_index(std::vector<double> &cost, std::vector<int> &index, int low, int high);
		static std::vector<std::vector<double> > Get_random_parameters();
		static void Initiate(std::vector<std::vector<std::bitset<384> > > &population,std::vector<double> &mean_squared);
		static void shuffle(std::vector<int> &index);
		static void tournament(std::vector<std::vector<std::bitset<384> > > &population,std::vector<double> &mean_squared);
		static void reproduction(std::vector<std::vector<std::bitset<384> > > &population,std::vector<double> &mean_squared);
		static void mutate(std::vector<std::vector<std::bitset<384> > > &population,std::vector<double> &mean_squared);
		static void rankChromosomes(std::vector<std::vector<std::bitset<384> > > &population,std::vector<double> &mean_squared);
		static double percentDifference(std::vector<std::bitset<384> > &individual1, std::vector<std::bitset<384> > &individual2);
		static double getDiversity(std::vector<std::vector<std::bitset<384> > > &population);
		static void DivergenceError();
		static void BottleneckError();
		static void CheckDiversity(std::vector<std::vector<std::bitset<384> > > &population);
		static void show_mean_squared(double mean_squared);
	public:
		static void run();
};

class anneal : private optimization
{
	private:
		static int partition(std::vector<double> &value, int low, int high);
		static void quicksort_x(std::vector<double> &value, int low, int high);
		static double Gaussian_move(double mean, double std_dev,int accepted);
		static std::vector<std::vector<double> > neighbor(std::vector<std::vector<double> > &state0,double error,int accepted);
		static double Temperature(double new_energy, int accepted);
		static double rand_double();
	public:
		static void run(std::vector<std::vector<double> > &old_state);
};

//generates a random double between -1.0 and 1.0
double optimization::RandInit()
{
	double r = rand() % RAND_MAX + (RAND_MAX/2 - 1);
	r /= RAND_MAX;
	return r;
};

//generates a random bool
bool optimization::rand_bool()
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
double optimization::Chebyshev(double x, std::vector<double> &param, double a, double b)
{
	int ni = param.size();
	double b1 = 0.f, b2 = 0.f;
	double y = (2.f*x-a-b)/(b-a);
	double y2 = 2.f*y;
	for (int i=ni-1; i>0; --i)
	{
		double temp = b1;
		b1 = y2*b1-b2+param[i];
		b2 = temp;
	}
	return y*b1-b2+0.5*param[0];
}

//combines material properties
double optimization::combine(double x, double y)
{
	return sqrt(0.5*(x*x + y*y));
}

//returns the approximate COR with independent variables x[] and coefficients parameters[][]
double optimization::f(std::vector<double> &x, std::vector<std::vector<double> > &param)
{
	double y1 = pow(Chebyshev(x[0],param[0],0.00018,1.03),0.5);
	y1 /= pow(Chebyshev(x[2],param[1],0.002701,370.f),0.5);
	y1 /= pow(Chebyshev(x[4],param[2],1.9,8.553),0.5);
	double y2 = pow(Chebyshev(x[1],param[0],0.00018,1.03),0.5);
	y2 /= pow(Chebyshev(x[3],param[1],0.002701,370.f),0.5);
	y2 /= pow(Chebyshev(x[5],param[2],1.9,8.553),0.5);
	double E = pow(1.f/(0.5/(y1*y1)+0.5/(y2*y2)),0.5);
	double v = pow(Chebyshev(x[6],param[3],0.f,6.f),0.5);
	double e = E*v*v;
	if (e > 1.f)
		return 1.f;
	return e;
}

//gets the residual of each datapoint
std::vector<double> optimization::GetResiduals(std::vector<double> &y, std::vector<std::vector<double> > &x, std::vector<std::vector<double> > &param)
{
	std::vector<double> residuals(n_data);
	for (int i=0; i<n_data; ++i)
	{
		double yi = f(x[i],param);
		residuals[i] = y[i]-yi;
		if (isnan(residuals[i]))
			residuals[i] = 1.f;
	}
	return residuals;
}

//returns the mean of the square of each residual
double optimization::Mean_square_error(std::vector<double> residuals)
{
	double sum = 0;
	int ni = residuals.size();
	for (int i=0; i<ni; ++i)
	{
		sum += residuals[i]*residuals[i];
	}
	return sum/n_data;
}

#endif /* COR_OPTIMIZATION_HPP_ */
