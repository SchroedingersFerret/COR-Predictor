/*
 * COR-common.hpp
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
    
#ifndef COR_COMMON_HPP_
#define COR_COMMON_HPP_

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
	
	for (int i=0; i<n_data; ++i)
	{
		if (y[i]<0.f || y[i]>1.f)
		{
			std::cout << "Error: The dependent variables must be between 0.0 and 1.0.\n";
			abort();
		}
	}
}

//asks user whether to begin optimization with completely random population
bool cor::Query_random()
{
	std::cout << "Initiate with random values? (Convergence will take longer)\nEnter y/n: ";
	char input;
	std::cin >> input;
	std::cout << "\n";
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
bool cor::Query_initiate()
{
	std::cout << "Initiate with these values?\nEnter y/n: ";
	char input;
	std::cin >> input;
	std::cout << "\n";
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
bool cor::Use_random()
{
	std::ifstream fin;
	fin.open("cor_parameters.csv");
	fin.close();
	if (fin.fail())
	{
		std::cout << "File: 'cor_parameters.csv' not found.\n\n";
		return Query_random();
	}
	
	std::cout << "File: 'cor_parameters.csv' found.\n\n";
	return Query_initiate();
}

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
	std::cout << "\n";
}

//asks the user whether to write the new parameters to file
bool cor::Query_write()
{
	std::cout << "Write these parameters to 'cor_parameters.csv'?\nPrevious values will be overwritten.\nEnter y/n: ";
	char input;
	std::cin >> input;
	std::cout << "\n";
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

#endif /* COR_COMMON_HPP_ */
