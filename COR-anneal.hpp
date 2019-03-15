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

//returns a random number from an gaussian distribution
double anneal::Gaussian_move(double mean, double error)
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
	return mean + mean*error*v/u;
}
	
//returns neighboring state
parameters anneal::neighbor(parameters state0,double error)
{
	parameters state1;
	int ni = state0.c.size();
	int nj = state0.c[0].size();
	for (int i=0; i<ni; ++i)
	{
		for (int j=0; j<nj; ++j)
		{
			state1.c[i][j] = Gaussian_move(state0.c[i][j],error);
		}
	}
	return state1;
}

//returns temperature given a change in energy and entropy
double anneal::Temperature(double initial_temperature, int accepted)
{
	return initial_temperature*exp(-accepted/4);
}

//runs simulated annealing to make sure predictions are in the accepted range
parameters anneal::run(parameters old_state)
{
	std::vector<std::vector<double> > x_rand = random_points(range(independent));
	double old_energy = SumOfSquares(GetResiduals(dependent,independent,old_state));
	double initial_temperature = FLT_MAX*old_energy;
	double old_temperature = initial_temperature;
	int accepted = 0;
	int iterations = 0;
	while (old_temperature > 0 && iterations < 1000)
	{
		parameters new_state = neighbor(old_state,error);
		double new_energy = SumOfSquares(GetResiduals(dependent,independent,new_state));
		
		double delta_energy = new_energy-old_energy;
		double new_temperature = Temperature(initial_temperature,accepted);
		
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
	return old_state;
}

#endif /* COR_ANNEAL_HPP_ */
