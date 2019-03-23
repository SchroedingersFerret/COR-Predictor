/*
 * COR-anneal.cpp
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
double anneal::Gaussian_move(double mean, double error, int accepted)
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
	return mean + error*v/u*1.f/(n_data*(1+accepted));
}
	
//returns neighboring state
std::vector<std::vector<double> > anneal::neighbor(std::vector<std::vector<double> > state0,double error,int accepted)
{
	int ni = state0.size();
	int nj = state0[0].size();
	std::vector<std::vector<double> > state1 (ni,std::vector<double> (nj));
	for (int i=0; i<ni; ++i)
	{
		for (int j=0; j<nj; ++j)
		{
			state1[i][j] = Gaussian_move(state0[i][j],error,accepted);
		}
	}
	return state1;
}

//returns temperature given a change in energy and entropy
double anneal::Temperature(double initial_temperature, int accepted)
{
	return initial_temperature*exp(-sqrt(accepted));
}

//generates a random double between 0 and 1.0
double anneal::rand_double()
{
	return ((double) rand() / (RAND_MAX));
}

//runs simulated annealing to make sure predictions are in the accepted range
std::vector<std::vector<double> > anneal::run(std::vector<std::vector<double> > old_state)
{
	double old_energy = Mean_square_error(GetResiduals(dependent,independent,old_state));
	double initial_temperature = FLT_MAX*old_energy;
	double old_temperature = initial_temperature;
	int accepted = 0;
	int iterations = 0;
	
	while (old_temperature > 0 && iterations < 10000)
	{
		
		std::vector<std::vector<double> > new_state = neighbor(old_state,error,accepted);
		double new_energy = Mean_square_error(GetResiduals(dependent,independent,new_state));
		
		double delta_energy = new_energy-old_energy;
		double new_temperature = Temperature(old_temperature,accepted);
		
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
