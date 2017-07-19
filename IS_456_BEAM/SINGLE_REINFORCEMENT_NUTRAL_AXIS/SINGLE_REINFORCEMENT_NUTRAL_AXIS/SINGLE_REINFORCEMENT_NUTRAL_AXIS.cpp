// SINGLE_REINFORCEMENT_NUTRAL_AXIS.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<stdio.h>
#include <iostream>
#include <vector>

using namespace std;

class BEAM
{
public:
double nutral_axis(int fst, int Ast, int fck, int b) 
{
		return (0.87*fst*Ast) / (0.36*fck*b);
	}
};


int main()
{
	vector<vector<int>> v2d;


	vector<double>v1;
	vector<double>::iterator it;
	it = v1.begin();
	cout << &it ;
	int k = 0;

	vector<vector<vector<int>>> matrix3d(3, vector<vector<int>>(6, vector<int>(10, 0)));
	vector<vector<int>> matirx2d(3, vector<int>(6, 0));


	

	BEAM beam1;
	for (int j = 100; j <= 500; j = j + 100)
	{
		cout << endl << "Breadth: " << j;
		for (int i = 100; i <= 1000; i = i + 100)
		{
			v1.push_back(beam1.nutral_axis(415, i, 25, j));
		}
	}
	
	getchar();
    return 0;
}

