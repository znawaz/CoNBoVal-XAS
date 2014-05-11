/*
 * point.h
 *
 *  Created on: Feb 25, 2014
 *      Author: zubair
 */

#ifndef POINT_H_
#define POINT_H_

#include <iostream>
#include <string>
#include <sstream>				// for ostream

using namespace std;

class point {
public:
	double x;
	double y;
	double z;
	int type;
	//Default Constructor
	point();

	//overload Constructor
	point(double, double, double, int);

	void print(ostream &out);

	//Destructor
	~point();



};


#endif /* POINT_H_ */
