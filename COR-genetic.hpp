/*
 * COR-genetic.hpp
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
    
#ifndef COR_GENETIC_HPP_
#define COR_GENETIC_HPP_

//encodes the parameters into an offset binary array
genome genetic::encode(parameters param)
{
	genome w;
	int ni = param.c.size();
	int nj = param.c[0].size();
	for (int i=0; i<ni; ++i)
	{
		for (int j=0; j<nj; ++j)
		{
			double sum = param.c[i][j];
			w.chromosome[i][j*64] = 1;
			if (param.c[i][j] < 0)
			{
				w.chromosome[i][j*64] = 0;
				sum *= -1;
			}
			w.chromosome[i][j*64+1] = 0;
			if ((int)(0.5+sum)==1)
				w.chromosome[i][j*64+1] = 1;
			double d = 2.f;
			for (int k=2; k<64; ++k)
			{
				if (w.chromosome[i][j*64+k-1])
					sum -= 1.f/d;
				w.chromosome[i][j*64+k] = 0;
				if ((int)(0.5+d*sum) == 1)
					w.chromosome[i][j*64+k] = 1;
				d *= 2;
			}
		}
	}
	return w;
}

//recovers the parameters from a binary chromosome
parameters genetic::decode(genome w)
{
	parameters param;
	int ni = param.c.size();
	int nj = param.c[0].size();
	for (int i=0; i<ni; ++i)
	{
		for (int j=0; j<nj; ++j)
		{
			double d = 2.f;
			double sum = 0;
			for (int k=1; k<64; ++k)
			{
				if (w.chromosome[i][j*64+k])
					sum += 1.f/d;
				d *= 2;
			}
			param.c[i][j] = sum + 1.f/d;
			if (!w.chromosome[i][j*64])
				param.c[i][j] *= -1;
		}
	}

	return param;
}

//fills a parameter array with random doubles
parameters genetic::Get_random_parameters()
{
	parameters param;
	int ni = param.c.size();
	int nj = param.c[0].size();
	for (int i=0; i<ni; ++i)
		for (int j=0; j<nj; ++j)
			param.c[i][j] = RandInit();
		
	return param;
}

//partition function for quicksort_index
int genetic::partition(std::vector<double> &cost, std::vector<int> &index, int low, int high)
{
	double pivot = cost[low];
	int i = low;
	for (int j=low+1; j<high; ++j)
	{
		if (cost[j] <= pivot)
		{
			i++;
			std::swap(cost[i],cost[j]);
			std::swap(index[i],index[j]);
		}
	}
	std::swap(cost[i],cost[low]);
	std::swap(index[i],index[low]);
	return i;
}

//quicksorts indices by cost
void genetic::quicksort_index(std::vector<double> &cost, std::vector<int> &index, int low, int high)
{
	if (low < high)
	{
		int pi = partition(cost, index, low, high);
		quicksort_index(cost, index, low, pi);
		quicksort_index(cost, index, pi+1, high);
	}
}

void genetic::Initiate(std::vector<genome> &population,std::vector<double> &squareSums)
{
	std::vector<genome> bin (n_initial);
	parameters param;
	std::vector<double> cost (n_initial);
	std::vector<int> index (n_initial);
	//fills the elite population with the parameters read from file unless user specifies an entirely random population
	if (!random_parameters)
		param = parameters_global;
	
	for (int i=0; i<n_elite; ++i)
	{
		if (random_parameters)
			param = Get_random_parameters();
		
		std::vector<double> residuals = GetResiduals(dependent,independent,param);
		
		cost[i] = SumOfSquares(residuals);
		index[i] = i;
		bin[i] = encode(param);
	}
	random_parameters = false;
	//The remaining population is initialized randomly
	for (int i=n_elite; i<n_initial; ++i)
	{
		param = Get_random_parameters();
		
		std::vector<double> residuals = GetResiduals(dependent,independent,param);
		cost[i] = SumOfSquares(residuals);
		index[i] = i;
		bin[i] = encode(param);
	}
	//sorts population by cost
	quicksort_index(cost,index,0,cost.size());
	for (int i=0; i<n_gpool; ++i)
	{
		squareSums[i] = cost[i];
		population[i] = bin[index[i]];
	}
	
}

//shuffles the indices into a random configuration
void genetic::shuffle(std::vector<int> &index)
{
	int k = index.size();
	for (int i = 0; i<k; ++i)
	{
		int j = rand() % k;
		if (j!= i)
			std::swap(index[i],index[j]);
	}
}

//performs tournament selection on the chromosome population
void genetic::tournament(std::vector<genome> &population,std::vector<double> &squareSums)
{
	std::vector<int> index(n_gpool);
	std::vector<genome> bin = population;
	std::vector<double> cost = squareSums;
	for (int i=0; i<n_gpool; ++i)
		index[i] = i;
	shuffle(index);
	int k = 0;
	for (int i=0; i<n_repro; ++i)
	{
		if (cost[index[k]] < cost[index[k+1]])
		{
			population[i] = bin[index[k]];
			squareSums[i] = cost[index[k]];
		}
		else
		{
			population[i] = bin[index[k+1]];
			squareSums[i] = cost[index[k+1]];
		}
		k += 2;
	}
}



//performs uniform crossover reproduction on the chromosomes
void genetic::reproduction(std::vector<genome> &population)
{
	int k = 0;
	for (int i=n_repro; i<n_repro+n_repro/2; ++i)
	{
		//perform reproduction ( ͡° ͜ʖ ͡°)
		int nl = population[i].chromosome.size();
		int nm = population[i].chromosome[0].size();
		for (int l=0; l<nl; ++l)
		{
			for (int m=0; m<nm; ++m)
			{
				bool parent = rand_bool();
				if (parent)
				{
					population[i].chromosome[l][m] = population[k].chromosome[l][m];
					population[i+n_repro/2].chromosome[l][m] = population[k+1].chromosome[l][m];
				}
				if (!parent)
				{
					population[i].chromosome[l][m] = population[k+1].chromosome[l][m];
					population[i+n_repro/2].chromosome[l][m] = population[k].chromosome[l][m];
				}
			}
		}
	}
}

//ranks the chromosomes by cost
void genetic::rankChromosomes(std::vector<genome> &population, std::vector<double> &squareSums)
{
	std::vector<genome> bin = population;
	parameters param;
	std::vector<double> cost(n_gpool);
	std::vector<int> index(n_gpool);
	
	for (int i=0; i<n_gpool; ++i)
	{
		param = decode(bin[i]);
		std::vector<double> residuals = GetResiduals(dependent,independent,param);
		cost[i] = SumOfSquares(residuals);
		index[i] = i;
	}
	
	quicksort_index(cost,index,0,cost.size());
	for (int i=0; i<n_gpool; ++i)
	{
		squareSums[i] = cost[i];
		population[i] = bin[index[i]];
	}
}

//performs mutations on the chromosome population
void genetic::mutate(std::vector<genome> &population,std::vector<double> &squareSums)
{
	//elite population remains unchanged if mutation increases the cost
	int ni = population[0].chromosome.size();
	int nj = population[0].chromosome[0].size();
	int m_elite = (int) (n_elite*ni*nj*pm+1);
	int m_normal = (int) ((n_gpool-n_elite)*ni*nj*pm+1);
	for (int i=0; i<m_elite; ++i)
	{
		int i_elite = rand() % n_elite;
		int i_mutate = rand() % ni;
		int j_mutate = rand() % nj;
		genome bin = population[i_elite];
		bin.chromosome[i_mutate].flip(j_mutate);
		parameters param = decode(bin);
		std::vector<double> residuals = GetResiduals(dependent,independent,param);
		double cost = SumOfSquares(residuals);
		if (cost < squareSums[i_elite])
		{
			population[i_elite] = bin;
			squareSums[i_elite] = cost;
		}
	}
	//cost increases are accepted in the remaining population
	for (int i=0; i<m_normal; ++i)
	{
		int i_gpool = rand() % (n_gpool-n_elite) + n_elite;
		int i_mutate = rand() % ni;
		int j_mutate = rand() % nj;
		population[i_gpool].chromosome[i_mutate].flip(j_mutate);
	}
	rankChromosomes(population,squareSums);
}

//returns the percentage of differing bits between two chromosomes
double genetic::percentDifference(genome individual1, genome individual2)
{
	int ni = individual1.chromosome.size();
	int nj = individual1.chromosome[0].size();
	double pd = 0.f;
	for (int i=0; i<ni; ++i)
	{
		for (int j=0; j<nj; ++j)
		{
			if (individual1.chromosome[i][j]!=individual2.chromosome[i][j])
				pd++;
		}
	}
	pd /= (ni*nj);
	return pd;
}

//returns the diversity of the population
double genetic::getDiversity(std::vector<genome> &population)
{
	double diversity = 0.f;
	for (int i=1; i<n_gpool; ++i)
	{
		diversity += percentDifference(population[0],population[i]);
	}
	diversity /= n_gpool-1;
	return diversity;
}

//aborts the program if the population diverges
void genetic::DivergenceError()
{
	std::cout << "\nError: Population divergence.\n";
	abort();
}

//aborts the program if the population bottlenecks
void genetic::BottleneckError()
{
	std::cout << "\nError: Population bottleneck.\n";
	abort();
}

//checks the diversity of the population and aborts if it is too large or small
void genetic::CheckDiversity(std::vector<genome> &population)
{
	double diversity = getDiversity(population);
	if (diversity > 0.75)
		DivergenceError();
	if (diversity < 0.001)
		BottleneckError();
}

//displays least squares regression
void genetic::show_least_squares(double S)
{
	std::cout << std::scientific << "\rS = " << S << std::flush;
}

//execute genetic algorithm
void genetic::run()
{
	//population stores each binary genome
	std::vector<genome> population(n_gpool);
	//sum of the square of each residual
	std::vector<double> squareSums(n_gpool);
	anneal ANNEAL;
	Initiate(population,squareSums);
	std::cout << "Running genetic algorithm...\n\n";
	int iterations = 0;
	double old_S = 0;
	double new_S = 0;
	while(squareSums[0] > error)
	{
		tournament(population,squareSums);
		
		reproduction(population);
		
		rankChromosomes(population,squareSums);
		
		mutate(population,squareSums);
		
		new_S = squareSums[0];
		if (new_S != old_S)
		{
			show_least_squares(squareSums[0]);
			old_S = new_S;
		}
		
		iterations++;
		if (iterations >= 200)
		{
			CheckDiversity(population);
			iterations = 0;
			parameters param = ANNEAL.run(decode(population[0]));
			double cost = SumOfSquares(GetResiduals(dependent,independent,param));
			population[0] = encode(param);
			squareSums[0] = cost;
		}
	};
	show_least_squares(squareSums[0]);
	parameters_global = decode(population[0]);
}

#endif /* COR_GENETIC_HPP_ */
