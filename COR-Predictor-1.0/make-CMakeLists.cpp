/*
 * make-CMakeLists.cpp
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
 
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include "COR-predictor.hpp"
#include "COR-lib.hpp"

int main()
{
	std::vector<std::vector<double> > parameters = COR.Get_parameters();
	std::ofstream fout;
	fout.open("function/CMakeLists.txt");
	fout << "cmake_minimum_required (VERSION 3.1)\n";
	fout << "project (Restitution)\n\n";
	
	for (int i=0; i<nx; ++i)
	{
		for (int j=0; j<nx; ++j)
			fout << "set (P" << i << j << " " << parameters[i][j] << ")\n";
	}
	fout << "\n";
	fout << "configure_file(\n";
	fout << "\"${PROJECT_SOURCE_DIR}/Restitution.hpp.in\"\n";
	fout << "\"${PROJECT_BINARY_DIR}/Restitution.hpp\"\n)";
	std::cout << "'CMakeLists.txt' written in '~/COR-Predictor-1.0/function'\n";
	return 0;
}

