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
void Get_settings()
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
std::vector<double> read_csv1d(const char * filename)
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
std::vector<std::vector<double> > read_csv2d(const char * filename)
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

//writes 1d vector to csv file
void write_csv1d(std::vector<double> a, const char * filename)
{
	std::ofstream fout;
	fout.open(filename);
	int ni = a.size();
	for (int i=0; i<ni; ++i)
		fout << a[i] << "\n";
	fout.close();
}

//writes 2d vector to csv file
void write_csv2d(std::vector<std::vector<double> > a, const char * filename)
{
	std::ofstream fout;
	fout.open(filename);
	int ni = a.size();
	int nj = a[0].size();
	for (int i=0; i<ni; ++i)
	{
		for (int j=0; j<nj; ++j)
		{
			fout << a[i][j];
			if (j != nj-1)
				fout << ",";
		}
		fout << "\n";
	}
	fout.close();
}

//reads the independent variables of the training datapoints
void Get_x()
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
void Get_y()
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

//user enters new training datapoints
void Point_entry()
{
	double input;
	std::vector<double> new_point(nx);
	std::cout << "Enter the yield strength of the first object.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	new_point[0] = input;
	
	std::cout << "Enter the yield strength of the second object.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	new_point[1] = input;
	
	std::cout << "Enter the Young's modulus of the first object.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	new_point[2] = input;
	
	std::cout << "Enter the Young's modulus of the second object.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	new_point[3] = input;
	
	std::cout << "Enter the density of the first object.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	new_point[4] = input;
	
	std::cout << "Enter the density of the second object.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	new_point[5] = input;
	
	std::cout << "Enter the objects' relative velocity.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	new_point[6] = input;
	
	x.push_back(new_point);
	
	std::cout << "Enter the coefficient of restitution.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	
	y.push_back(input);
}

//asks if user wants to enter more points
bool Enter_more()
{
	std::cout << "Enter more datapoints?\nEnter y/n: ";
	char input;
	std::cin.get(input);
	std::cin.get();
	std::cout << "\n";
	switch(input)
	{
		case 'y':
		case 'Y':	return true;
					break;

		case 'n':
		case 'N':	abort();
					break;

		default	:	return Enter_more();
					break;
	}
}

//asks user to return to main menu or quit program
bool Return_quit()
{
	int input;
	std::cout << "Enter '1' to return to the main menu. Enter '2' to quit.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	switch(input)
	{
		case 1:	return false;
					break;

		case 2:	return true;
					break;

		default	:	return Return_quit();
					break;
	}
}
	
//training datapoint entry
void Enter()
{
	bool enter_point = true;
	while (enter_point)
	{
		Point_entry();
		enter_point = Enter_more();
	}
	write_csv2d(x,"cor_x.csv");
	write_csv1d(y,"cor_y.csv");
	quit_cor = Return_quit();
}
	
//asks user whether to initialize elite population with parameters read from file
bool Query_initiate()
{
	std::cout << "Initiate with these values?\nEnter y/n: ";
	char input;
	std::cin.get(input);
	std::cin.get();
	std::cout << "\n";
	switch(input)
	{
		case 'y':
		case 'Y':	return false;
					break;

		case 'n':
		case 'N':	std::cout << "Initiating with random values. (Convergence will take longer)\n\n";
					return true;
					break;

		default	:	return Query_initiate();
					break;
	}
}

//returns a boolean operator to instruct program whether to randomize elite population
bool Use_random()
{
	std::ifstream fin;
	fin.open("cor_parameters.csv");
	fin.close();
	if (fin.fail())
	{
		std::cout << "File: 'cor_parameters.csv' not found.\n";
		std::cout << "Initiating with random values. (Convergence will take longer)\n\n";
		return true;
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
void Print_parameters(parameters param)
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
bool Query_write()
{
	std::cout << "Write these parameters to 'cor_parameters.csv'?\nPrevious values will be overwritten.\nEnter y/n: ";
	char input;
	std::cin.get(input);
	std::cin.get();
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
void Write_parameters(parameters param)
{
	bool write = Query_write();
	if (write)
	{
		write_csv2d(param.c,"cor_parameters.csv");
		std::cout << "Parameters written.\n\n";
	}
	else
		std::cout << "Optimization terminated without writing new parameters.\n\n";
}

void Optimize()
{
	random_parameters = Use_random();
	clock_t tStart = clock();
	std::thread t1(&GENETIC.run);
	t1.join();
	ANNEAL.run(param);
	std::cout << "\n\nParameters found:\n\n" ;
	Print_parameters(param);
	std::cout << "Execution time: " << ( (double) clock()-tStart)/CLOCKS_PER_SEC << " s\n\n";
	Write_parameters(param);
	quit_cor = Return_quit();
}

void Show_main_menu()
{
	std::cout << "Enter '1' to enter a new training datapoint.\n\n";
	std::cout << "Enter '2' to optimize the parameters.\n\n";
	std::cout << "Enter '3' to predict a coefficient of restitution.\n\n";
	std::cout << "Enter '4' to quit.\n\n";
}

struct Mode
{
	public:
	enum Enum
	{
		ENTER = 1,
		OPTIMIZE = 2,
		PREDICT = 3,
		QUIT = 4
	};
	Enum e;
	void Set_mode();
};	

void Mode::Set_mode()
{
	char input;
	std::cin.get(input);
	std::cin.get();
	std::cout << "\n";
	switch(input)
	{
		case '1' :  e = ENTER;
					break;
		case '2' : e = OPTIMIZE;
					break;
		case '3' : e = PREDICT;
					break;
		case '4' : e = QUIT;
					break;
		default : std::cout << "Invalid input. Please enter '1','2','3', or '4'.\n\n";
					Set_mode();
					break;
	}
}

void Main_menu()
{
	Show_main_menu();
	Mode mode;
	mode.Set_mode();
	switch(mode.e)
	{
		case mode.Enum::ENTER 	 :  Enter();
									break;
		case mode.Enum::OPTIMIZE :	Optimize();
									break;
		case mode.Enum::PREDICT  : ;
									break;
		case mode.Enum::QUIT   	 : 	quit_cor = true;
									break;	
		default : break;
	}
}

#endif /* COR_COMMON_HPP_ */
