/**
 *	 Project: Coordination Number and Bond Valence for XAS
 *
 *
 *   Copyright (C) 2014 SESAME | SYNCHROTRON-LIGHT FOR EXPERIMENTAL SCIENCE AND
 *   											APPLICATIONS IN THE MIDDLE EAST
 *                           Allan, jordan
 *
 *   Principal authors: Zubair Nawaz (zubair.nawaz@sesame.org.jo)
 *   Last revision: 08/05/2014
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as published
 *   by the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   and the GNU Lesser General Public License  along with this program.
 *   If not, see <http://www.gnu.org/licenses/>.
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
