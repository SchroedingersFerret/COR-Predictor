/*
 * COR-predictor.hpp
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

#ifndef COR_PREDICTOR_HPP_
#define COR_PREDICTOR_HPP_

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <time.h>
#include <bitset>
#include <stdlib.h>

//initial population size
static int n_initial;
//size of gene pool for selection
static int n_gpool;
//number of chromosomes to be selected for reproduction
static int n_repro;
//number of independent variables
static const int nx = 7;
//percentage selected for mutation
static double pm;
//number of elites
static int n_elite;
//number of datapoints
static int n_data;
//least squares error tolerance
static double error;
//independent variables array size n_data*nx
static std::vector<std::vector<double> > x;
//dependent variable array size n_data
static std::vector<double> y;
//chromosome stores parameters as a binary array
static std::vector<std::vector<std::bitset<64*nx> > > chromosomes(0, std::vector<std::bitset<64*nx> > (nx));
//sum of the square of each residual
static std::vector<double> squareSums;

class cor
{
	private:
		double NextRand();
		bool rand_bool();
		double taylor(double x, std::vector<double> parameters);
		double f(std::vector<double> x, std::vector<std::vector<double> > parameters);
		bool Query_write();
	public:
		void Get_settings();
		void Get_x();
		void Get_y();
		std::vector<std::vector<double> > Get_parameters();
		std::vector<std::vector<double> > Get_random_parameters();
		class genetic
		{
			private:
				bool Query_random();
				bool Query_initiate();
				bool Use_random();
				std::vector<std::bitset<64*nx> > encode(std::vector<std::vector<double> > parameters);
				std::vector<double> GetResiduals(std::vector<double> y, std::vector<std::vector<double> > x, std::vector<std::vector<double> > parameters);
				int partition(std::vector<double> &cost, std::vector<int> &index, int low, int high);
				void quicksort_index(std::vector<double> &cost, std::vector<int> &index, int low, int high);
				void shuffle(std::vector<int> &index);
			public:
				std::vector<std::vector<double> > decode(std::vector<std::bitset<64*nx> > w);
				double SumOfSquares(std::vector<double> residuals);
				void Initiate();
				void tournament();
				void reproduction();
				void rankChromosomes();
				void mutate();
		};
		void Write_parameters();
};

static cor COR;
static cor::genetic GENETIC;
#endif /* COR_PREDICTOR_HPP_ */
