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
#include "COR-cor.hpp"
#include "COR-anneal.hpp"
#include "COR-genetic.hpp"
#include "COR-common.hpp"

int main()
{
	srand((unsigned int)time(NULL));
	Get_settings();
	Get_independent();
	Get_dependent();
	Get_parameters();
	std::cout << "Welcome to COR Predictor 0.4\n";
	std::cout << "Copyright 2019, J. Ball (SchroedingersFerret)\n\n";
	while (!quit_cor)
		Main_menu();
	return 0;
}
