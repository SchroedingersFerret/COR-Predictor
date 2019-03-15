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
			output.push_back((double)input);
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
			datapoint.push_back((double)input);
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
void Get_independent()
{
	std::ifstream fin;
	fin.open("cor_independent.csv");
	fin.close();
	if (fin.fail())
	{
		std::cout << "Error: File 'cor_independent.csv' not found.\n";
		abort();
	}
	
	independent = read_csv2d("cor_independent.csv");
	
	if (independent[0].size() != 7)
	{
		std::cout << "Error: File 'cor_independent.csv' must be of dimension n*7.\n";
		abort();
	}
	n_data = independent.size();
	fin.close();
}
	
//reads the dependent variables of the training datapoints
void Get_dependent()
{
	std::ifstream fin;
	fin.open("cor_dependent.csv");
	fin.close();
	if (fin.fail())
	{
		std::cout << "Error: File 'cor_dependent.csv' not found.\n";
		abort();
	}
	
	dependent = read_csv1d("cor_dependent.csv");

	if (dependent.size() != independent.size())
	{
		std::cout << "Error: Files 'cor_independent.csv' and 'cor_dependent.csv' must have the same number of entries.\n";
		abort();
	}
	
	for (int i=0; i<n_data; ++i)
	{
		if (dependent[i]<0.f || dependent[i]>1.f)
		{
			std::cout << "Error: The dependent variables must be between 0.0 and 1.0.\n";
			abort();
		}
	}
}

//reads parameter array from file
void Get_parameters()
{
	std::ifstream fin;
	fin.open("cor_parameters.csv");
	fin.close();
	if (!fin.fail())
	{
		std::vector<std::vector<double> > temp = read_csv2d("cor_parameters.csv");
		if (temp.size() != parameters_global.c.size() || temp[0].size() != parameters_global.c[0].size())
		{
			std::cout << "Error: File 'cor_parameters.csv' must be of dimensions " << parameters_global.c.size() << "*" << parameters_global.c[0].size() << ".\n";
			abort();
		}
	parameters_global.c = temp;	
	}
}

//user enters new training datapoints
void Point_entry()
{
	double input;
	std::vector<double> new_point(nx);
	std::cout << "Enter the yield strength of the first object in GPa.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	new_point[0] = input;
	
	std::cout << "Enter the yield strength of the second object in GPa.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	new_point[1] = input;
	
	std::cout << "Enter the Young's modulus of the first object in GPa.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	new_point[2] = input;
	
	std::cout << "Enter the Young's modulus of the second object in GPa.\n";
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
	
	std::cout << "Enter the objects' relative velocity in m/s.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	new_point[6] = input;
	
	independent.push_back(new_point);
	
	std::cout << "Enter the coefficient of restitution.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	
	dependent.push_back(input);
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
		case 'N':	return false;
					break;

		default	:	return Enter_more();
					break;
	}
}

//asks user to return to main menu or quit program
bool Return_quit()
{
	char input;
	std::cout << "Enter '1' to return to the main menu. Enter '2' to quit.\n";
	std::cin.get(input);
	std::cin.get();
	std::cout << "\n";
	switch(input)
	{
		case '1':	return false;
					break;

		case '2':	return true;
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
		write_csv2d(independent,"cor_independent.csv");
		write_csv1d(dependent,"cor_dependent.csv");
		enter_point = Enter_more();
	}
	
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
	genetic GENETIC;
	random_parameters = Use_random();
	clock_t tStart = clock();
	GENETIC.run();
	std::cout << "\n\nParameters found:\n\n" ;
	Print_parameters(parameters_global);
	std::cout << std::fixed << "Execution time: " << ( (double) clock()-tStart)/CLOCKS_PER_SEC << " s\n\n";
	Write_parameters(parameters_global);
	quit_cor = Return_quit();
}

std::vector<double> pGet_independent()
{
	std::vector<double> x(nx);
	double input;
	std::cout << "Enter the yield strength of the first object in GPa.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	x[0] = input;
	
	std::cout << "Enter the yield strength of the second object in GPa.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	x[1] = input;
	
	std::cout << "Enter the Young's modulus of the first object in GPa.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	x[2] = input;
	
	std::cout << "Enter the Young's modulus of the second object in GPa.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	x[3] = input;
	
	std::cout << "Enter the density of the first object.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	x[4] = input;
	
	std::cout << "Enter the density of the second object.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	x[5] = input;
	
	std::cout << "Enter the objects' relative velocity in m/s.\n";
	std::cin >> input;
	std::cin.get();
	std::cout << "\n";
	x[6] = input;
	
	return x;
}	

void Predict()
{
	std::vector<double> pred_x = pGet_independent();
	cor COR;
	double pred_y = COR.f(pred_x,parameters_global);
	std::cout << "e = " << pred_y << "\n\n";
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
		case mode.Enum::PREDICT  :  Predict();
									break;
		case mode.Enum::QUIT   	 : 	quit_cor = true;
									break;	
		default : break;
	}
}

#endif /* COR_COMMON_HPP_ */
