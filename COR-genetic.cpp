/*
 * COR-genetic.cpp
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
    
//encodes the parameters into an offset binary array
std::vector<std::bitset<64*5> > genetic::encode(std::vector<std::vector<double> > param)
{
	std::vector<std::bitset<64*5> > w(4);
	int ni = param.size();
	int nj = param[0].size();
	for (int i=0; i<ni; ++i)
	{
		for (int j=0; j<nj; ++j)
		{
			double sum = param[i][j];
			w[i][j*64] = 1;
			if (param[i][j] < 0)
			{
				w[i][j*64] = 0;
				sum *= -1;
			}
			w[i][j*64+1] = 0;
			if ((int)(0.5+sum)==1)
				w[i][j*64+1] = 1;
			double d = 2.f;
			for (int k=2; k<64; ++k)
			{
				if (w[i][j*64+k-1])
					sum -= 1.f/d;
				w[i][j*64+k] = 0;
				if ((int)(0.5+d*sum) == 1)
					w[i][j*64+k] = 1;
				d *= 2;
			}
		}
	}
	return w;
}

//recovers the parameters from a binary chromosome
std::vector<std::vector<double> > genetic::decode(std::vector<std::bitset<64*5> > w)
{
	std::vector<std::vector<double> > param(4,std::vector<double> (5));
	int ni = param.size();
	int nj = param[0].size();
	for (int i=0; i<ni; ++i)
	{
		for (int j=0; j<nj; ++j)
		{
			double d = 2.f;
			double sum = 0;
			for (int k=1; k<64; ++k)
			{
				if (w[i][j*64+k])
					sum += 1.f/d;
				d *= 2;
			}
			param[i][j] = sum + 1.f/d;
			if (!w[i][j*64])
				param[i][j] *= -1;
		}
	}

	return param;
}

//fills a parameter array with random doubles
std::vector<std::vector<double> > genetic::Get_random_parameters()
{
	std::vector<std::vector<double> > param(4,std::vector<double> (5));
	int ni = param.size();
	int nj = param[0].size();
	for (int i=0; i<ni; ++i)
		for (int j=0; j<nj; ++j)
			param[i][j] = RandInit();
		
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

//initiates genetic algorithm
void genetic::Initiate(std::vector<std::vector<std::bitset<64*5> > > &population,std::vector<double> &mean_squared)
{
	std::vector<std::vector<std::bitset<64*5> > > bin (n_initial, std::vector<std::bitset<64*5> > (4));
	std::vector<std::vector<double> > param(4,std::vector<double> (5));
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
		
		cost[i] = Mean_square_error(residuals);
		index[i] = i;
		bin[i] = encode(param);
	}
	random_parameters = false;
	//The remaining population is initialized randomly
	for (int i=n_elite; i<n_initial; ++i)
	{
		param = Get_random_parameters();
		
		std::vector<double> residuals = GetResiduals(dependent,independent,param);
		cost[i] = Mean_square_error(residuals);
		index[i] = i;
		bin[i] = encode(param);
	}
	//sorts population by cost
	quicksort_index(cost,index,0,cost.size());
	for (int i=0; i<n_gpool; ++i)
	{
		mean_squared[i] = cost[i];
		population[i] = bin[index[i]];
	}
	
}

//shuffles the indices into a random configuration
void genetic::shuffle(std::vector<int> &index)
{
	for (int i = 0; i<n_gpool; ++i)
	{
		int j = rand() % n_gpool;
		if (j!= i)
			std::swap(index[i],index[j]);
	}
}

//performs tournament selection on the chromosome population
void genetic::tournament(std::vector<std::vector<std::bitset<64*5> > > &population,std::vector<double> &mean_squared)
{
	std::vector<int> index(n_gpool);
	std::vector<std::vector<std::bitset<64*5> > > bin = population;
	std::vector<double> cost = mean_squared;
	for (int i=0; i<n_gpool; ++i)
		index[i] = i;
	shuffle(index);
	int k = 0;
	for (int i=0; i<n_repro; ++i)
	{
		if (cost[index[k]] < cost[index[k+1]])
		{
			population[i] = bin[index[k]];
			mean_squared[i] = cost[index[k]];
		}
		else
		{
			population[i] = bin[index[k+1]];
			mean_squared[i] = cost[index[k+1]];
		}
		k += 2;
	}
}

//performs uniform crossover reproduction on the chromosomes
void genetic::reproduction(std::vector<std::vector<std::bitset<64*5> > > &population, std::vector<double> &mean_squared)
{
	int k = 0;
	int n_repro2 = n_repro/2;
	int nl = population[0].size();
	int nm = population[0][0].size();
	for (int i=n_repro; i<n_repro+n_repro2; ++i)
	{
		//perform reproduction ( ͡° ͜ʖ ͡°)
		for (int l=0; l<nl; ++l)
		{
			for (int m=0; m<nm; ++m)
			{
				bool parent = rand_bool();
				if (parent)
				{
					population[i][l][m] = population[k][l][m];
					population[i+n_repro2][l][m] = population[k+1][l][m];
				}
				if (!parent)
				{
					population[i][l][m] = population[k+1][l][m];
					population[i+n_repro2][l][m] = population[k][l][m];
				}
			}
		}
	}
	rankChromosomes(population,mean_squared);
}

//ranks the chromosomes by cost
void genetic::rankChromosomes(std::vector<std::vector<std::bitset<64*5> > > &population, std::vector<double> &mean_squared)
{
	std::vector<std::vector<std::bitset<64*5> > > bin = population;
	std::vector<std::vector<double> > param(4,std::vector<double> (5));
	std::vector<double> cost(n_gpool);
	std::vector<int> index(n_gpool);
	
	for (int i=0; i<n_gpool; ++i)
	{
		param = decode(bin[i]);
		std::vector<double> residuals = GetResiduals(dependent,independent,param);
		cost[i] = Mean_square_error(residuals);
		index[i] = i;
	}
	
	quicksort_index(cost,index,0,cost.size());
	for (int i=0; i<n_gpool; ++i)
	{
		mean_squared[i] = cost[i];
		population[i] = bin[index[i]];
	}
}

//mutates bits in the elite population
void genetic::mutate_elite(std::vector<std::vector<std::bitset<64*5> > > &population,std::vector<double> &mean_squared)
{
	//elite population remains unchanged if mutation increases the cost
	int ni = population[0].size();
	int nj = population[0][0].size();
	int m_elite = (int) (n_elite*ni*nj*pm+1);
	for (int i=0; i<m_elite; ++i)
	{
		int i_elite = rand() % n_elite;
		int i_mutate = rand() % ni;
		int j_mutate = rand() % nj;
		std::vector<std::bitset<64*5> > bin = population[i_elite];
		bin[i_mutate].flip(j_mutate);
		std::vector<std::vector<double> > param = decode(bin);
		std::vector<double> residuals = GetResiduals(dependent,independent,param);
		double cost = Mean_square_error(residuals);
		if (cost < mean_squared[i_elite])
		{
			population[i_elite] = bin;
			mean_squared[i_elite] = cost;
		}
	}
}

//mutates bits in the remainder of the population
void genetic::mutate_normal(std::vector<std::vector<std::bitset<64*5> > > &population,std::vector<double> &mean_squared)
{
	int ni = population[0].size();
	int nj = population[0][0].size();
	int m_normal = (int) (n_normal*ni*nj*pm+1);
	//cost increases are accepted in the remaining population
	for (int i=0; i<m_normal; ++i)
	{
		int i_gpool = rand() % n_normal + n_elite;
		int i_mutate = rand() % ni;
		int j_mutate = rand() % nj;
		population[i_gpool][i_mutate].flip(j_mutate);
	}
}

//returns the percentage of differing bits between two chromosomes
double genetic::percentDifference(std::vector<std::bitset<64*5> > individual1, std::vector<std::bitset<64*5> > individual2)
{
	int ni = individual1.size();
	int nj = individual1[0].size();
	double pd = 0.f;
	for (int i=0; i<ni; ++i)
	{
		for (int j=0; j<nj; ++j)
		{
			if (individual1[i][j]!=individual2[i][j])
				pd++;
		}
	}
	pd /= (ni*nj);
	return pd;
}

//returns the diversity of the population
double genetic::getDiversity(std::vector<std::vector<std::bitset<64*5> > > &population)
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
void genetic::CheckDiversity(std::vector<std::vector<std::bitset<64*5> > > &population)
{
	double diversity = getDiversity(population);
	if (diversity > 0.60)
		DivergenceError();
	if (diversity < 0.0001)
		BottleneckError();
}

//displays least squares regression
void genetic::show_mean_squared(double S)
{
	std::cout << std::scientific << "\rmean square error = " << S << std::flush;
}

//execute genetic algorithm
void genetic::run()
{
	//population stores each binary genome
	std::vector<std::vector<std::bitset<64*5> > > population(n_gpool,std::vector<std::bitset<64*5> > (4));
	//sum of the square of each residual
	std::vector<double> mean_squared(n_gpool);
	std::thread init(&Initiate,std::ref(population),std::ref(mean_squared));
	init.join();
	std::cout << "Running genetic algorithm...\n\n";
	int iterations = 0;
	double old_S = 0;
	double new_S = 0;
	while(mean_squared[0] > error)
	{
		std::thread tour(&tournament,std::ref(population),std::ref(mean_squared));
		tour.join();
		
		std::thread repro(&reproduction,std::ref(population),std::ref(mean_squared));
		repro.join();
		
		new_S = mean_squared[0];
		if (new_S != old_S)
		{
			show_mean_squared(mean_squared[0]);
			old_S = new_S;
		}
		
		iterations++;
		if (iterations >= 50)
		{
			std::thread check([](std::vector<std::vector<std::bitset<64*5> > > &population,std::vector<double> &mean_squared)
			{
				CheckDiversity(population);
				std::vector<std::vector<double> > param = anneal::run(decode(population[0]));
				double cost = Mean_square_error(GetResiduals(dependent,independent,param));
				population[n_elite-1] = encode(param);
				mean_squared[n_elite-1] = cost;
				rankChromosomes(population,mean_squared);
			},std::ref(population),std::ref(mean_squared));
			check.join();
			iterations = 0;
		}
		
		std::thread mu1(&mutate_elite,std::ref(population),std::ref(mean_squared));
		std::thread mu2(&mutate_normal,std::ref(population),std::ref(mean_squared));
		mu1.join();
		mu2.join();
	};
	show_mean_squared(mean_squared[0]);
	parameters_global = decode(population[0]);
}
