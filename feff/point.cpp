/*
 * point.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: zubair
 */
#include "point.h"

point::point()	{
	x = 0.0;
	y = 0.0;
	z = 0.0;
	type = -1;
}

point::point(double X, double Y, double Z, int Type)	{
	x = X;
	y = Y;
	z = Z;
	type = Type;
}

void point::print(ostream &out)	{
	out << x << " " << y << " " << z << endl;
}

point::~point()	{

}



