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
//sum of the square of each residual
static std::vector<double> squareSums;

class genome
{
	public:
		std::vector<std::bitset<64*10> > chromosome;
		genome()
		{
			chromosome.resize(4);
		};
};

//population stores each binary genome
static std::vector<genome> population;

class parameters
{
	public:
		std::vector<std::vector<double> > c;
		parameters()
		{
			c.resize(4,std::vector<double> (10));
		};
};

class cor
{
	private:
		double NextRand();
		bool rand_bool();
		double combine(double x, double y);
		double Chebyshev(double x, std::vector<double> param);
		std::vector<double> read_csv1d(const char * filename);
		std::vector<std::vector<double> > read_csv2d(const char * filename);
		bool Query_write();
	public:
		double f(std::vector<double> x, parameters param);
		void Get_settings();
		void Get_x();
		void Get_y();
		parameters Get_parameters();
		parameters Get_random_parameters();
		class genetic
		{
			private:
				bool Query_random();
				bool Query_initiate();
				bool Use_random();
				genome encode(parameters param);
				std::vector<double> GetResiduals(std::vector<double> y, std::vector<std::vector<double> > x, parameters param);
				int partition(std::vector<double> &cost, std::vector<int> &index, int low, int high);
				void quicksort_index(std::vector<double> &cost, std::vector<int> &index, int low, int high);
				void shuffle(std::vector<int> &index);
				std::vector<double> nearbyPoint(std::vector<double> x);
				double nearbyPointTest(std::vector<double> x, parameters param);
				double percentDifference(genome individual1, genome individual2);
				double getDiversity();
				void DivergenceError();
				void BottleneckError();
			public:
				parameters decode(genome w);
				double SumOfSquares(std::vector<double> residuals);
				void Initiate();
				void tournament();
				void reproduction();
				void rankChromosomes();
				void mutate();
				void CheckDiversity();
		};
		void Print_parameters(parameters param);
		void Write_parameters(parameters param);
};

static cor COR;
static cor::genetic GENETIC;
#endif /* COR_PREDICTOR_HPP_ */
