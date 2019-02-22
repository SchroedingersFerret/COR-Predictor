/*
 * COR-lib.hpp
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
    
#ifndef COR_LIB_HPP_
#define COR_LIB_HPP_

//generates a random double between 0 and 1
double cor::NextRand()
{
	double r = rand() % RAND_MAX + (RAND_MAX/2 - 1);
	r /= RAND_MAX;
	return r;
}

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

//Returns a nearby point
std::vector<double> cor::genetic::nearbyPoint(std::vector<double> x)
{
	std::vector<double> newx(nx);
	for (int i=0; i<nx; ++i)
	{
		newx[i] = x[i]+x[i]*COR.NextRand();
	}
	return newx;
}

//Returns a modifier if nearby points are not in the acceptible range for the COR
double cor::genetic::nearbyPointTest(std::vector<double> x, parameters param)
{
	double y = COR.f(nearbyPoint(x),param);
	if (y > 1.f)
		return y-1.f;
	if (y < 0)
		return y;
	return 0;
}

//gets the residual of each datapoint
std::vector<double> cor::genetic::GetResiduals(std::vector<double> y, std::vector<std::vector<double> > x, parameters param)
{
	std::vector<double> residuals(2*n_data);
	for (int i=0; i<n_data; ++i)
	{
		residuals[i] = y[i]-COR.f(x[i],param);
		//residuals[n_data+i] = nearbyPointTest(x[i],param);
	}
	return residuals;
}

//returns the sum of the square of each residual
double cor::genetic::SumOfSquares(std::vector<double> residuals)
{
	double sum;
	int ni = residuals.size();
	for (int i=0; i<ni; ++i)
		sum += residuals[i]*residuals[i];
	return sum;
}

//encodes the parameters into an offset binary array
genome cor::genetic::encode(parameters param)
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
parameters cor::genetic::decode(genome w)
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

//partition function for quicksort_index
int cor::genetic::partition(std::vector<double> &cost, std::vector<int> &index, int low, int high)
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
void cor::genetic::quicksort_index(std::vector<double> &cost, std::vector<int> &index, int low, int high)
{
	if (low < high)
	{
		int pi = partition(cost, index, low, high);
		quicksort_index(cost, index, low, pi);
		quicksort_index(cost, index, pi+1, high);
	}
}

//reads and applies settings
void cor::Get_settings()
{
	std::ifstream fin;
	fin.open("settings.txt");
	if (fin.fail())
	{
		std::cout << "Error: File 'settings.txt' not found.\n";
		abort();
	}
	char ch;
	do
	{
		fin.get(ch);
	}while(ch!='=');
	fin >> n_initial;

	do
	{
		fin >> ch;
	}while(ch!='=');
	fin >> n_gpool;
	population.resize(n_gpool);
	squareSums.resize(n_gpool);
	n_repro = n_gpool/2;

	do
	{
		fin >> ch;
	}while(ch!='=');
	fin >> n_elite;

	do
	{
		fin >> ch;
	}while(ch!='=');
	fin >> pm;

	do
	{
		fin >> ch;
	}while(ch!='=');
	fin >> error;

	fin.close();
	
}

//reads 1d csv files
std::vector<double> cor::read_csv1d(const char * filename)
{
	std::ifstream fin;
	double input;
	std::vector<double> output;
	fin.open(filename);
	while(!fin.eof())
	{
		if (fin.peek()==',')
		{
			char ch;
			fin.get(ch);
		}
		if (fin >> input)
			output.push_back(input);
	}
	fin.close();
	return output;
}

//reads 2d csv files
std::vector<std::vector<double> > cor::read_csv2d(const char * filename)
{
	std::ifstream fin;
	double input;
	std::vector<double> datapoint;
	std::vector<std::vector<double> > output;
	fin.open(filename);
	while(!fin.eof())
	{
		
		if (fin.peek() == ',' || fin.peek() == '\n')
		{
			char ch;
			fin.get(ch);
			if (ch == '\n')
			{
				output.push_back(datapoint);
				datapoint.clear();
			}
		}
		if (fin >> input)
		{
			datapoint.push_back(input);
		}
	}
	fin.close();
	return output;
}

//reads the independent variables of the training datapoints
void cor::Get_x()
{
	std::ifstream fin;
	fin.open("cor_x.csv");
	fin.close();
	if (fin.fail())
	{
		std::cout << "Error: File 'cor_x.csv' not found.\n";
		abort();
	}
	
	x = read_csv2d("cor_x.csv");
	
	if (x[0].size() != 7)
	{
		std::cout << "Error: File 'cor_x.csv' must be of dimension n*7.\n";
		abort();
	}
	n_data = x.size();
	fin.close();
}
	
//reads the dependent variables of the training datapoints
void cor::Get_y()
{
	std::ifstream fin;
	fin.open("cor_y.csv");
	fin.close();
	if (fin.fail())
	{
		std::cout << "Error: File 'cor_y.csv' not found.\n";
		abort();
	}
	
	y = read_csv1d("cor_y.csv");

	if (y.size() != x.size())
	{
		std::cout << "Error: Files 'cor_x.csv' and 'cor_y.csv' must have the same number of entries.\n";
		abort();
	}
}

//asks user whether to begin optimization with completely random population
bool cor::genetic::Query_random()
{
	std::cout << "Initiate with random values? (Convergence will take longer)\nEnter y/n: ";
	char input;
	std::cin >> input;
	switch(input)
	{
		case 'y':
		case 'Y':	return true;
					break;

		case 'n':
		case 'N':	abort();
					break;

		default	:	return Query_random();
					break;
	}
}

//asks user whether to initialize elite population with parameters read from file
bool cor::genetic::Query_initiate()
{
	std::cout << "Initiate with these values?\nEnter y/n: ";
	char input;
	std::cin >> input;
	switch(input)
	{
		case 'y':
		case 'Y':	return false;
					break;

		case 'n':
		case 'N':	return Query_random();
					break;

		default	:	return Query_initiate();
					break;
	}
}

//returns a boolean operator to instruct program whether to randomize elite population
bool cor::genetic::Use_random()
{
	std::ifstream fin;
	fin.open("cor_parameters.csv");
	fin.close();
	if (fin.fail())
	{
		std::cout << "File: 'cor_parameters.csv' not found.\n";
		return Query_random();
	}
	
	std::cout << "File: 'cor_parameters.csv' found.\n";
	return Query_initiate();
}

//reads parameter array from file
parameters cor::Get_parameters()
{
	parameters param;

	
		std::vector<std::vector<double> > temp = read_csv2d("cor_parameters.csv");
		if (temp.size() != param.c.size() || temp[0].size() != param.c[0].size())
		{

			std::cout << "Error: File 'cor_parameters.csv' must be of dimensions " << param.c.size() << "*" << param.c[0].size() << ".\n";
			abort();
		}
		param.c = temp;
	
	return param;
}

//fills a parameter array with random doubles
parameters cor::Get_random_parameters()
{
	parameters param;
	int ni = param.c.size();
	int nj = param.c[0].size();
	for (int i=0; i<ni; ++i)
		for (int j=0; j<nj; ++j)
			param.c[i][j] = NextRand();
		
	return param;
}

void cor::genetic::Initiate()
{
	std::vector<genome> bin (n_initial);
	parameters param;
	std::vector<double> cost (n_initial);
	std::vector<int> index (n_initial);
	bool random = Use_random();
	//fills the elite population with the parameters read from file unless user specifies an entirely random population
	if (!random)
		param = COR.Get_parameters();
	
	for (int i=0; i<n_elite; ++i)
	{
		if (random)
			param = COR.Get_random_parameters();
		
		std::vector<double> residuals = GetResiduals(y,x,param);
		
		cost[i] = SumOfSquares(residuals);
		index[i] = i;
		bin[i] = encode(param);
	}

	//The remaining population is initialized randomly
	for (int i=n_elite; i<n_initial; ++i)
	{
		param = COR.Get_random_parameters();
		
		std::vector<double> residuals = GetResiduals(y,x,param);
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
void cor::genetic::shuffle(std::vector<int> &index)
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
void cor::genetic::tournament()
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
void cor::genetic::reproduction()
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
				bool parent = COR.rand_bool();
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
void cor::genetic::rankChromosomes()
{
	std::vector<genome> bin = population;
	parameters param;
	std::vector<double> cost(n_gpool);
	std::vector<int> index(n_gpool);
	
	for (int i=0; i<n_gpool; ++i)
	{
		param = decode(bin[i]);
		std::vector<double> residuals = GetResiduals(y,x,param);
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
void cor::genetic::mutate()
{
	
	std::vector<genome> bin = population;
	//elite population remains unchanged if mutation increases the cost
	int ni = bin[0].chromosome.size();
	int nj = bin[0].chromosome[0].size();
	int mmax = (int) (n_gpool*ni*nj*64*pm+1);
	for (int i=0; i<n_elite; ++i)
	{
		parameters param;
		
		int n_mutate = (int) (ni*nj*pm+1);
		for (int l=0; l<n_mutate; ++l)
		{
			
			int i_mutate = rand() % ni;
			int j_mutate = rand() % nj;
			bin[i].chromosome[i_mutate].flip(j_mutate);
		}
		param = decode(bin[i]);
		std::vector<double> residuals = GetResiduals(y,x,param);
		double cost = SumOfSquares(residuals);
		if (cost < squareSums[i])
		{
			population[i] = bin[i];
			squareSums[i] = cost;
		}
	}
	//cost increases are accepted in the remaining population
	for (int i=0; i<mmax; ++i)
	{
		int i_gpool = rand() % (n_gpool-n_elite) + n_elite;
		int i_mutate = rand() % ni;
		int j_mutate = rand() % nj;
		population[i_gpool].chromosome[i_mutate].flip(j_mutate);
	}
	rankChromosomes();
}

//returns the percentage of differing bits between two chromosomes
double cor::genetic::percentDifference(genome individual1, genome individual2)
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
double cor::genetic::getDiversity()
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
void cor::genetic::DivergenceError()
{
	std::cout << "Error: Population divergence.\n";
	abort();
}

//aborts the program if the population bottlenecks
void cor::genetic::BottleneckError()
{
	std::cout << "Error: Population bottleneck.\n";
	abort();
}

//checks the diversity of the population and aborts if it is too large or small
void cor::genetic::CheckDiversity()
{
	double diversity = getDiversity();
	if (diversity > 0.50)
		DivergenceError();
	if (diversity < 0.01)
		BottleneckError();
}

//prints the parameters in the terminal
void cor::Print_parameters(parameters param)
{
	int ni = param.c.size();
	int nj = param.c[0].size();
	for (int i=0; i<ni; ++i)
	{
		for (int j=0; j<nj; ++j)
		{
			std::cout << param.c[i][j] << ",";
		}
		std::cout << "\n";
	}
	std::cout << "\n\n";
}

//asks the user whether to write the new parameters to file
bool cor::Query_write()
{
	std::cout << "Write these parameters to 'cor_parameters.csv'?\nPrevious values will be overwritten.\nEnter y/n: ";
	char input;
	std::cin >> input;
	switch(input)
	{
		case 'y':
		case 'Y':	return true;
					break;

		case 'n':
		case 'N':	return false;
					break;

		default	:	return Query_write();
					break;
	}
}
	
//runs Query_write() and writes parameters to file depending on user input
void cor::Write_parameters(parameters param)
{
	bool write = Query_write();
	if (write)
	{
		std::ofstream fout;
		fout.open("cor_parameters.csv");
		int ni = param.c.size();
		int nj = param.c[0].size();
		for (int i=0; i<ni; ++i)
		{
			for (int j=0; j<nj; ++j)
			{
				fout << param.c[i][j];
				if (j != nj-1)
					fout << ",";
			}
			fout << "\n";
		}
		fout.close();
		
		std::cout << "Parameters written.\n";
	}
	else
		std::cout << "Program terminated without writing new parameters.\n";
}

#endif /* COR_LIB_HPP_ */
