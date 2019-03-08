/*
 * COR-anneal.hpp
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
    
#ifndef COR_ANNEAL_HPP_
#define COR_ANNEAL_HPP_

//partition function for quicksort_x
int anneal::partition(std::vector<double> &value, int low, int high)
{
	double pivot = value[low];
	int i = low;
	for (int j=low+1; j<high; ++j)
	{
		if (value[j] <= pivot)
		{
			i++;
			std::swap(value[i],value[j]);
		}
	}
	std::swap(value[i],value[low]);
	return i;
}

//quicksorts values by cost
void anneal::quicksort_x(std::vector<double> &value, int low, int high)
{
	if (low < high)
	{
		int pi = partition(value, low, high);
		quicksort_x(value, low, pi);
		quicksort_x(value, pi+1, high);
	}
}

//returns the range of the indendent variables
std::vector<std::vector<double> > anneal::range(std::vector<std::vector<double> > x)
{
	std::vector<std::vector<double> > range(nx, std::vector<double> (2));
	std::vector<int> indices(n_data);
	for (int i=0; i<nx; ++i)
	{
		std::vector<double> value(n_data);
		for (int j=0; j<n_data; ++j)
			value[j] = x[j][i];
		quicksort_x(value,0,n_data);
		range[i][0] = value[0];
		range[i][1] = value[n_data-1];
	}
	return range;
}

//returns a set of random training points
std::vector<std::vector<double> > anneal::random_points(std::vector<std::vector<double> > range)
{
	int ni = 10*n_data;
	std::vector<std::vector<double> > points(ni, std::vector<double> (nx));
	for (int i=0; i<ni; ++i)
	{
		for (int j=0; j<nx; ++j)
		{
			double r = rand_double();
			points[i][j] = range[j][0]+(range[j][1]-range[j][0])*r;
		}
	}
	return points;
}

//returns a random number from an gaussian distribution
double anneal::Gaussian_move(double mean, double T)
{
	double u,v,x,xx;
	do
	{
		u = RandInit();
		v = RandInit();
		x = 1.71552776992141*(v-0.5)/u;
		xx = x*x;
	}while(xx >= 5.f-5.13610166675097*u && 
		(xx <= 1.036961042583566/u+1.4 || xx <= -4*log(u)));
	double g = 1.f/(2.f*(v/u+T)*(1.f+1.f/(T)));
	return mean + g;
}
	
//returns neighboring state
parameters anneal::neighbor(parameters state0,double temperature)
{
	parameters state1;
	int ni = state0.c.size();
	int nj = state0.c[0].size();
	for (int i=0; i<ni; ++i)
	{
		for (int j=0; j<nj; ++j)
		{
			state1.c[i][j] = Gaussian_move(state0.c[i][j],temperature);
		}
	}
	return state1;
}

//gets the residual of each datapoint
std::vector<double> anneal::GetResiduals(std::vector<std::vector<double> > x_rand, parameters state)
{
	int ni = x_rand.size();
	std::vector<double> residuals(ni);
	for (int i=0; i<ni; ++i)
	{
		double yi = f(x_rand[i],state);
		residuals[i] = 0.f;
		if ((int) yi >= 1)
			residuals[i] = yi - 1.f;
		if (yi < 0 )
			residuals[i] = yi;
	}
	
	return residuals;
}

//returns the sum of the square of each residual
double anneal::SumOfSquares(std::vector<double> residuals)
{
	double sum;
	int ni = residuals.size();
	for (int i=0; i<ni; ++i)
		sum += residuals[i]*residuals[i];
	return sum;
}

//returns temperature given a change in energy and entropy
double anneal::Temperature(double new_energy,int accepted)
{
	return FLT_MAX*new_energy*exp(-accepted);
}

//runs simulated annealing to make sure predictions are in the accepted range

//This one is iffy. It takes the parameters returned by the genetic algorithm
//and uses the objective function to test whether a set of random points 
//give an answer outside [0,1]. So far, this has never happened beyond sets of
//random unoptimized parameters. Still, I feel there needs to be some safeguard
//against unphysical results, so I want the annealing to tweak the results into 
//an acceptable range without destroying what the genetic algorithm does. 
//We'll see if it actually does that if it ever needs to be called.
void anneal::run(parameters old_state)
{
	std::vector<std::vector<double> > x_rand = random_points(range(x));
	double old_energy = SumOfSquares(GetResiduals(x_rand,old_state));
	double old_temperature = FLT_MAX*old_energy;
	int accepted = 0;
	int iterations = 0;
	while (old_temperature > 0 && iterations < 1000)
	{
		parameters new_state = neighbor(old_state,old_temperature);
		double new_energy = SumOfSquares(GetResiduals(x_rand,new_state));
		double delta_energy = new_energy-old_energy;
		double new_temperature = Temperature(new_energy,accepted);
		double P = rand_double();
		double probability;
		if (delta_energy < 0)
			probability = 1.f;
		else
			probability = exp(-delta_energy/new_temperature);
		if (P < probability)
		{
			old_state = new_state;
			old_energy = new_energy;
			old_temperature = new_temperature;
			accepted++;
		}
		iterations++;
	}
	if (old_energy > 0)
		run(old_state);
	else
		param = old_state;
}

#endif /* COR_ANNEAL_HPP_ */
