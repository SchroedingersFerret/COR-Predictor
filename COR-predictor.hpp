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
#include <float.h>
#include <time.h>
#include <bitset>
#include <stdlib.h>
#include <thread>

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
//determines whether initial population contains entirely random parameters
bool random_parameters = false;
//quits the program
bool quit_cor = false;

class genome
{
	public:
		std::vector<std::bitset<64*10> > chromosome;
		genome()
		{
			chromosome.resize(4);
		};
};

class parameters
{
	public:
		std::vector<std::vector<double> > c;
		parameters()
		{
			c.resize(4,std::vector<double> (10));
		};
};

//parameter array
static parameters param;

class cor
{
	public:	
		double rand_double();
		bool rand_bool();
		double RandInit();
		double combine(double x, double y);
		double Chebyshev(double x, std::vector<double> param);
		double f(std::vector<double> x, parameters param);	
};

class genetic : public cor
{
	
	private:
		std::vector<double> GetResiduals(std::vector<double> y, std::vector<std::vector<double> > x, parameters param);
		double SumOfSquares(std::vector<double> residuals);
		genome encode(parameters param);
		parameters decode(genome w);
		int partition(std::vector<double> &cost, std::vector<int> &index, int low, int high);
		void quicksort_index(std::vector<double> &cost, std::vector<int> &index, int low, int high);
		parameters Get_parameters();
		parameters Get_random_parameters();
		void Initiate(std::vector<genome> &population,std::vector<double> &squareSums);
		void shuffle(std::vector<int> &index);
		void tournament(std::vector<genome> &population,std::vector<double> &squareSums);
		void reproduction(std::vector<genome> &population);
		void rankChromosomes(std::vector<genome> &population,std::vector<double> &squareSums);
		void mutate(std::vector<genome> &population,std::vector<double> &squareSums);
		double percentDifference(genome individual1, genome individual2);
		double getDiversity(std::vector<genome> &population);
		void DivergenceError();
		void BottleneckError();
		void CheckDiversity(std::vector<genome> &population);
		void show_least_squares(double S);
	public:
		static void run();
};

class anneal : public cor
{
	private:
		int partition(std::vector<double> &value, int low, int high);
		void quicksort_x(std::vector<double> &value, int low, int high);
		std::vector<std::vector<double> > range(std::vector<std::vector<double> > x);
		std::vector<std::vector<double> > random_points(std::vector<std::vector<double> > range);
		double Gaussian_move(double mean, double std_dev);
		parameters neighbor(parameters state0,double energy);
		std::vector<double> GetResiduals(std::vector<std::vector<double> > x_rand, parameters state);
		double SumOfSquares(std::vector<double> residuals);
		double Temperature(double new_energy, int accepted);
	public:
		void run(parameters old_state);
};

void Get_settings();
std::vector<double> read_csv1d(const char * filename);
std::vector<std::vector<double> > read_csv2d(const char * filename);
void write_csv1d(std::vector<double> a, const char * filename);
void write_csv2d(std::vector<std::vector<double> > a, const char * filename);
void Get_x();
void Get_y();
void Point_entry();
bool Enter_more();
bool Return_quit();
void Enter();
bool Query_initiate();
bool Use_random();
bool Query_write();
void Print_parameters(parameters param);
void Write_parameters(parameters param);
void Optimize();
void Show_main_menu();
void Main_menu();

static genetic GENETIC;
static anneal ANNEAL;

#endif /* COR_PREDICTOR_HPP_ */

