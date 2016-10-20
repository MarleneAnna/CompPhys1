// CP1_1_1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <complex>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <random>
#include <time.h>
#include <stdio.h>

using namespace std;

#define pi 3.14159265358979323846

const int n = 4;

//global vectors containing data
vector<vector<double>> *ordata = new vector<vector<double>>();

//global vectors containing fourier (back)transformed data
vector<vector<double>> *fourierdata = new vector<vector<double>>();
vector<vector<double>> *backfourier = new vector<vector<double>>();

void DiscreteFourierTrafo(vector<vector<double>> *data, vector<vector<double>> *fourier, int n_input, int n_output, bool backtrafo);
void GenerateData(string file, int num);
void ReadData(string file, vector<vector<double>> *data, int num);
void WriteToFile(string file, vector<vector<double>> *fourier, int num);


int main()
{

	GenerateData("DataFile", 2*n);
	ReadData("DataFile", ordata, 2*n);

	//Calculations - Discrete Fourier Transform
	cout << "DATEN" << endl;

	for (int k = 0; k < n; k++)
	{
		cout << ordata->at(0).at(k) << " + i * " << ordata->at(1).at(k) << endl;
	}

	cout << endl;
	cout << "DFT Forward" << endl;
	DiscreteFourierTrafo(ordata, fourierdata, n, 8, 0);

	cout << endl;
	cout << "DFT Backward" << endl;
	DiscreteFourierTrafo(fourierdata, backfourier, n, 8, 1);

	WriteToFile("DFTForward", fourierdata, 2*n);
	WriteToFile("DFTBackward", backfourier, 2*n);

	delete ordata;
	delete fourierdata;
	delete backfourier;


	system("pause");
    return 0;
}

//Discrete Fourier Transformation (Forward or Backward)
void DiscreteFourierTrafo(vector<vector<double>> *data, vector<vector<double>> *fourier, int n_input, int n_output, bool backtrafo)
{
	fourier->resize(2);


	int t = -2;
	if (backtrafo) { t = 2; }

	for (int i = 0; i < n_output; i++)
	{
		double sumreal = 0.;
		double sumim = 0.;
		sumreal = 0.;
		sumim = 0.;

		for (int j = 0; j < n_input; j++)
		{
			sumreal += data->at(0).at(j)*cos(t * pi*i*j / n_input) - data->at(1).at(j)*sin(t * pi*i*j / n_input);
			sumim += data->at(0).at(j)*sin(t * pi*i*j / n_input) + data->at(1).at(j)*cos(t * pi*i*j / n_input);

		}

		if (backtrafo)
		{
			fourier->at(0).push_back(sumreal / n_input);
			fourier->at(1).push_back(sumim / n_input);
		}
		else
		{
			fourier->at(0).push_back(sumreal);
			fourier->at(1).push_back(sumim);
		}

		cout << fourier->at(0).at(i) << " + i * " << fourier->at(1).at(i) << endl;
	}
}

//Generate Data and write it to File (real \n im \n ... )
void GenerateData(string file, int num)
{
	mt19937 randgen;
	uniform_real_distribution<double> randdist(0., 1.);
	randgen.seed(time(NULL));

	ofstream datafile;
	datafile.open(file + ".txt");

	for (int i = 0; i < num; i++)
	{
		datafile << randdist(randgen) << endl;
	}

	datafile.close();
}

//Read Data from file
void ReadData(string file, vector<vector<double>> *data, int num)
{
	if (num % 2 != 0) throw invalid_argument("Odd number of data points!");
	if (num == 0) throw invalid_argument("Invalid number of data points!");

	data->resize(2);

	ifstream datafile;
	datafile.open(file + ".txt");
	string line;
	double entry;

	for(int i=0; i<num; i++)
	{
		datafile >> line;
		entry = atof(line.c_str());

		if (i % 2 == 0)
		{
			data->at(0).push_back(entry);
		}
		else
		{
			data->at(1).push_back(entry);
		}
		
	}

}

//Write transformed Data to file ( real + i * im )
void WriteToFile(string file, vector<vector<double>> *fourier, int num)
{
	ofstream outputfile;
	outputfile.open(file + ".txt");

	for (int i = 0; i < num; i++)
	{
		outputfile << fourier->at(0).at(i) << " + i * " << fourier->at(1).at(i) << "\n";
	}

	outputfile.close();
}