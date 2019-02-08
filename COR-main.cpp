/*
 * COR-main.cpp
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

#include "COR-predictor.hpp"
#include "COR-lib.hpp"

int main()
{
	std::cout << "Welcome to COR Predictor 0.3\n";
	std::cout << "Copyright 2019, J. Ball (SchroedingersFerret)\n\n";
	srand((unsigned int)time(NULL));

	COR.Get_settings();
	COR.Get_x();
	COR.Get_y();
	
	clock_t tStart = clock();
	//execute genetic algorithm
	GENETIC.Initiate();
	while(squareSums[0] > error)
	{
		GENETIC.tournament();
		
		GENETIC.reproduction();
		
		GENETIC.rankChromosomes();
		
		GENETIC.mutate();
		
		std::cout << "S = " << squareSums[0] << "\n";
		iterations++;
		if (iterations >= 100)
		{
			GENETIC.CheckDiversity();
			iterations = 0;
		}
	};
	std::cout << "S = " << squareSums[0] << "\n";
	std::vector<std::vector<double> > parameters = GENETIC.decode(chromosomes[0]);
	std::cout << "\nParameters found:\n" ;

	for (int i=0; i<n_par; ++i)
	{
		for (int j=0; j<nx; ++j)
		{
			std::cout << parameters[i][j] << ",";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
	std::cout << "Execution time: " << ( (double) clock()-tStart)/CLOCKS_PER_SEC << " s\n\n";
	COR.Write_parameters();
	return 0;
}



