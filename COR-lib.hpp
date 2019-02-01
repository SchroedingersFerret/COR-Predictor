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

//generates a random double between 0 and 1.0
double cor::NextRand()
{
	double r = rand() % RAND_MAX + 1;
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

//returns the taylor approximation of x with coefficients from parameters[]
double cor::taylor(double x, std::vector<double> parameters)
{
	double sum;
	for (int i=0; i<nx; ++i)
	{
		double product = parameters[i];
		for (int j=0; j<i; ++j)
			product *= x;
		sum += product;
	}
	return sum;
}
	
//returns the approximate COR with independent variables x[] and coefficients parameters[][]
double cor::f(std::vector<double> x, std::vector<std::vector<double> > parameters)
{
	double y = 3.1;
	y *= taylor(x[0],parameters[0]);
	y *= taylor(x[1],parameters[0]);
	y /= taylor(x[2],parameters[1]);
	y /= taylor(x[3],parameters[1]);
	y /= taylor(x[4],parameters[2]);
	y /= taylor(x[5],parameters[2]);
	y /= taylor(x[6],parameters[3]);
	return y;
}

//Performs a ratio test 
double cor::genetic::convergence(double x, std::vector<double> parameters)
{
	return (parameters[7]/parameters[6])*x;
}

//Returns a modifier if a the product of the taylor series convergences is divergent
double cor::genetic::convergenceTest(std::vector<double> x, std::vector<std::vector<double> > parameters)
{
	double r = 3.1;
	r *= convergence(x[0],parameters[0]);
	r *= convergence(x[1],parameters[0]);
	r /= convergence(x[2],parameters[1]);
	r /= convergence(x[3],parameters[1]);
	r /= convergence(x[4],parameters[2]);
	r /= convergence(x[5],parameters[2]);
	r /= convergence(x[6],parameters[3]);
	if (r*r > 1.f)
		return r;
	return 1.f;
}

//Returns a nearby point
std::vector<double> cor::genetic::nearbyPoint(std::vector<double> x)
{
	std::vector<double> newx(nx);
	for (int i=0; i<nx; ++i)
	{
		newx[i] = x[i]+x[i]*COR.NextRand()*0.01;
	}
	return newx;
}

//Returns a modifier if nearby points are not in the acceptible range for the COR
double cor::genetic::nearbyPointTest(std::vector<double> x, std::vector<std::vector<double> > parameters)
{
	double y = COR.f(nearbyPoint(x),parameters);
	if (y > 1.f)
		return y;
	if (y < 0)
		return y - 1.f;
	return 1.f;
}

//gets the residual of each datapoint
std::vector<double> cor::genetic::GetResiduals(std::vector<double> y, std::vector<std::vector<double> > x, std::vector<std::vector<double> > parameters)
{
	std::vector<double> residuals(n_data);
	for (int i=0; i<n_data; ++i)
		residuals[i] = y[i]-COR.f(x[i],parameters)*convergenceTest(x[i],parameters)*nearbyPointTest(x[i],parameters);
	return residuals;
}

//returns the sum of the square of each residual
double cor::genetic::SumOfSquares(std::vector<double> residuals)
{
	double sum;
	for (int i=0; i<n_data; ++i)
		sum += residuals[i]*residuals[i];
	return sum;
}

//encodes the parameters into an offset binary array
std::vector<std::bitset<64*nx> > cor::genetic::encode(std::vector<std::vector<double> > parameters)
{
	std::vector<std::bitset<64*nx> > w(n_par);
	for (int i=0; i<n_par; ++i)
	{
		for (int j=0; j<nx; ++j)
		{
			double sum = parameters[i][j];
			w[i][j*64] = 1;
			if (parameters[i][j] < 0)
			{
				w[i][j*64] = 0;
				sum *= -1;
			}
			w[i][j*64+1] = 0;
			if ((int)(0.5+sum)==1)
				w[i][j*64+1] = 1;
			double d = 2.0;
			for (int k=2; k<64; ++k)
			{
				if (w[i][j*64+k-1])
					sum -= 1.0/d;
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
std::vector<std::vector<double> > cor::genetic::decode(std::vector<std::bitset<64*nx> > w)
{
	std::vector<std::vector<double> > parameters(n_par,std::vector<double> (nx));
	for (int i=0; i<n_par; ++i)
	{
		for (int j=0; j<nx; ++j)
		{
			double d = 2.0;
			double sum = 0;
			for (int k=1; k<64; ++k)
			{
				if (w[i][j*64+k])
					sum += 1.0/d;
				d *= 2;
			}
			parameters[i][j] = sum + 1.0/d;
			if (w[i][j*64])
				parameters[i][j] *= -1;
		}
	}
	return parameters;
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
	chromosomes.resize(n_gpool);
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



//reads the independent variables of the training datapoints
void cor::Get_x()
{
	std::ifstream fin;
	fin.open("cor_x.csv");
	if (fin.fail())
	{
		std::cout << "Error: File 'cor_x.csv' not found.\n";
		abort();
	}
	double input;
	int in_i = 0;
	std::vector<double> datapoint(nx);

	while(!fin.eof())
	{
		if(fin.peek() == ',' || fin.peek() == '\n')
		{
			char ch;
			fin.get(ch);
			switch(ch)
			{
				case ','	:
								in_i++;
								if(in_i > nx-1)
								{

									std::cout << "Error: File 'cor_x.csv' must be of dimensions n*8.\n";
									fin.close();
									abort();
								}
								break;

				case '\n'	:	if (in_i == nx-1)
								{
									x.push_back(datapoint);
									in_i = 0;
									n_data++;
								}
								else
								{
									std::cout << "Error: File 'cor_x.csv' must be of dimensions n*8.\n";
									fin.close();
									abort();
								}
								break;

			}
		}
		if(fin >> input)
		{
			datapoint[in_i] = input;
		}

	}
	fin.close();
}

//reads the dependent variables of the training datapoints
void cor::Get_y()
{
	std::ifstream fin;
	fin.open("cor_y.csv");
	if (fin.fail())
	{
		std::cout << "Error: File 'cor_y.csv' not found.\n";
		abort();
	}
	double input;
	int in_i = 0;
	int in_j = 0;
	while(!fin.eof())
	{
		if(fin.peek() == ',')
		{
			char ch;
			fin.get(ch);
			switch(ch)
			{
				case ','	:
								in_i++;
								if(in_i > 0)
								{

									std::cout << "Error: File 'cor_y.csv' must be of dimension n*1'.\n";
									fin.close();
									abort();
								}
								break;
				case '\n'	:	if (in_i == 0)
								{
									in_i = 0;
									in_j++;
								}
								else
								{
									std::cout << "Error: File 'cor_y.csv' must be of dimensions n*1.\n";
									fin.close();
									abort();
								}
								break;
			}
		}
		if(fin >> input)
		{
			y.push_back(input);
		}
	}
	if (y.size() != x.size())
	{
		std::cout << "Error: Files 'cor_x.csv' and 'cor_y.csv' must have the same number of entries.\n";
		fin.close();
		abort();
	}
	fin.close();
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
	bool random;
	std::ifstream fin;
	fin.open("cor_parameters.csv");
	fin.close();
	if (fin.fail())
	{
		std::cout << "File: 'cor_parameters.csv' not found.\n";
		random = Query_random();
	}
	else
	{
		std::cout << "File: 'cor_parameters.csv' found.\n";
		random = Query_initiate();
	}

	return random;
}

//reads parameter array from file
std::vector<std::vector<double> > cor::Get_parameters()
{
	std::vector<std::vector<double> > parameters;
	std::ifstream fin;
	fin.open("cor_parameters.csv");
	double input;
	int in_i;
	std::vector<double> datapoint(n_par);
	while(!fin.eof())
	{
		if(fin.peek() == ',' || fin.peek() == '\n')
		{
			char ch;
			fin.get(ch);
			switch(ch)
			{
				case ','	:
								in_i++;
								if(in_i > nx-1)
								{

									std::cout << "Error: File 'cor_parameters.csv' must be of dimensions 8*8.\n";
									fin.close();
									abort();
								}
								break;

				case '\n'	:	if (in_i == nx-1)
								{
									parameters.push_back(datapoint);
									in_i = 0;

								}
								else
								{
									std::cout << "Error: File 'cor_parameters.csv' must be of dimensions 8*8.\n";
									fin.close();
									abort();
								}
								break;

			}
		}
		if(fin >> input)
		{
			datapoint[in_i] = input;
		}
	}
	fin.close();
	return parameters;
}

//fills a parameter array with random doubles
std::vector<std::vector<double> > cor::Get_random_parameters()
{
	std::vector<std::vector<double> > parameters(n_par,std::vector<double>(nx));
	for (int i=0; i<n_par; ++i)
		for (int j=0; j<nx; ++j)
			parameters[i][j] = NextRand();
	return parameters;
}

void cor::genetic::Initiate()
{
	std::vector<std::vector<std::bitset<64*nx> > > bin (n_initial, std::vector<std::bitset<64*nx> > (n_par));
	std::vector<std::vector<double> > parameters (n_par,std::vector<double> (nx));
	std::vector<double> cost (n_initial);
	std::vector<int> index (n_initial);
	bool random = Use_random();

	//fills the elite population with the parameters read from file unless user specifies an entirely random population
	if (!random)
		parameters = COR.Get_parameters();
	for (int i=0; i<n_elite; ++i)
	{
		if (random)
			parameters = COR.Get_random_parameters();
		std::vector<double> residuals = GetResiduals(y,x,parameters);
		cost[i] = SumOfSquares(residuals);
		index[i] = i;
		bin[i] = encode(parameters);
	}

	//The remaining population is initialized randomly
	for (int i=n_elite; i<n_initial; ++i)
	{
		parameters = COR.Get_random_parameters();
		std::vector<double> residuals = GetResiduals(y,x,parameters);
		cost[i] = SumOfSquares(residuals);
		index[i] = i;
		bin[i] = encode(parameters);
	}
	//sorts population by cost
	quicksort_index(cost,index,0,cost.size());
	for (int i=0; i<n_gpool; ++i)
	{
		squareSums[i] = cost[i];
		chromosomes[i] = bin[index[i]];
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
		{
			int i_tmp = index[i];
			index[i] = index[j];
			index[j] = i_tmp;
		}
	}
}

//performs tournament selection on the chromosome population
void cor::genetic::tournament()
{
	std::vector<int> index(n_gpool);
	std::vector<std::vector<std::bitset<64*nx> > > bin = chromosomes;
	std::vector<double> cost = squareSums;
	for (int i=0; i<n_gpool; ++i)
		index[i] = i;
	shuffle(index);
	int k = 0;
	for (int i=0; i<n_repro; ++i)
	{
		if (cost[index[k]] < cost[index[k+1]])
		{
			chromosomes[i] = bin[index[k]];
			squareSums[i] = cost[index[k]];
		}
		else
		{
			chromosomes[i] = bin[index[k+1]];
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
		for (int l=0; l<4; ++l)
		{
			for (int m=0; m<nx*64; ++m)
			{
				bool parent = COR.rand_bool();
				if (parent)
				{
					chromosomes[i][l][m] = chromosomes[k][l][m];
					chromosomes[i+n_repro/2][l][m] = chromosomes[k+1][l][m];
				}
				if (!parent)
				{
					chromosomes[i][l][m] = chromosomes[k+1][l][m];
					chromosomes[i+n_repro/2][l][m] = chromosomes[k][l][m];
				}
			}
		}
		k += 2;
	}
}

//ranks the chromosomes by cost
void cor::genetic::rankChromosomes()
{
	std::vector<std::vector<std::bitset<64*nx> > > bin = chromosomes;
	std::vector<std::vector<double> > parameters(n_par,std::vector<double> (nx));
	std::vector<double> cost(n_gpool);
	std::vector<int> index(n_gpool);
	for (int i=0; i<n_gpool; ++i)
	{
		parameters = decode(bin[i]);
		std::vector<double> residuals = GetResiduals(y,x,parameters);
		cost[i] = SumOfSquares(residuals);
		index[i] = i;
	}
	quicksort_index(cost,index,0,cost.size());
	for (int i=0; i<n_gpool; ++i)
	{
		squareSums[i] = cost[i];
		chromosomes[i] = bin[index[i]];
	}
}

//performs mutations on the chromosome population
void cor::genetic::mutate()
{
	int mmax = (int) (n_gpool*n_par*nx*64*pm+1);
	std::vector<std::vector<std::bitset<64*nx> > > bin = chromosomes;
	//elite population remains unchanged if mutation increases the cost
	for (int i=0; i<n_elite; ++i)
	{
		std::vector<std::vector<double> > parameters(n_par,std::vector<double> (nx));
		int n_mutate = (int) (n_par*nx*64*pm+1);
		for (int l=0; l<n_mutate; ++l)
		{
			int i_mutate = rand() % n_par;
			int j_mutate = rand() % (nx*64);
			bin[i][i_mutate].flip(j_mutate);
		}
		parameters = decode(bin[i]);
		std::vector<double> residuals = GetResiduals(y,x,parameters);
		double cost = SumOfSquares(residuals);
		if (cost < squareSums[i])
		{
			chromosomes[i] = bin[i];
			squareSums[i] = cost;
		}
	}
	//cost increases are accepted in the remaining population
	for (int i=0; i<mmax; ++i)
	{
		int i_gpool = rand() % (n_gpool-n_elite) + n_elite;
		int i_mutate = rand() % n_par;
		int j_mutate = rand() % (nx*64);
		chromosomes[i_gpool][i_mutate].flip(j_mutate);
	}
	rankChromosomes();
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
void cor::Write_parameters()
{
	bool write = Query_write();
	if (write)
	{
		std::ofstream fout;
		fout.open("cor_parameters.csv");
		std::vector<std::vector<double> > parameters = GENETIC.decode(chromosomes[0]);
		int ni = parameters.size();
		int nj = parameters[0].size();
		for (int i=0; i<ni; ++i)
		{
			for (int j=0; j<nj; ++j)
			{
				fout << parameters[i][j];
				if (j != nj-1)
					fout << ",";
			}
			fout << "\n";
		}
		fout.close();
		std::cout << "Parameters written to 'cor_parameters.csv'.\n";
	}
	else
		std::cout << "Program terminated without writing new parameters.\n";
}

#endif /* COR_LIB_HPP_ */
