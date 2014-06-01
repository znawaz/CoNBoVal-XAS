/**
 *   Project: Coordination Number and Bond Valence for XAS
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

#include <iostream>
#include <vector>
#include <map>
#include <math.h> 				// sqrt
#include <string>
#include <fstream>				// ifstream
#include <sstream>				// for ostream
#include <string>
#include <cstring>				// string manipulation
#include <limits>
#include <stdlib.h>
#include <ctype.h>				// tolower
#include <algorithm>			// for find
#include <iomanip>  			// for setw
#include <set>
#include <ANN/ANN.h>			// ANN declarations

#include "point.h"

using namespace std;

map<int, string> atom_list;			// map to store atom type and atom notation from input file

int				k;					// number of nearest neighbors
int				dim	= 3;			// dimension
double			eps	= 0;			// error bound

//----------------------------------------------------------------------
//	Parameters that are set in getArgs()
//----------------------------------------------------------------------
char inputFilename[40];
char outputFilenameDist[40];	// prefix name for the file storing distances
char outputFilenameVal[40];		// suffix file name for storing bond valence and coordination number

void getArgs(int argc, char **argv);			// get command-line arguments

/**
 * find and replace 'toReplace' in 's' with 'replaceWith'
 * @param s
 */
string myreplace(string &s, string toReplace,string replaceWith)
{
	size_t found = s.find(toReplace);
	if ( found!= string::npos)
		return(s.replace(found, toReplace.length(), replaceWith));
	return s;
}
/**
 * takes a number as string and convert to double.
 * Furthermore, in the input file, the exponent is written with D, so it replaces D with E.
 */
double convertToDouble(string str)	{
	myreplace(str,"D","E");
	double val = atof(str.c_str());
	return val;

}

/**
 * print header to every output distance file
 */
void printHeader(ostream &out)	{
	out << "* This feff input file was generated personal program based on Molecular Dynamic model" << endl;
	out << "* -----------------------------------------------------------------" << endl << endl << endl << endl;
	out << "HOLE 4   1.0 " << endl << endl;
	out << "*         pot    xsph  fms   paths genfmt ff2chi" << endl;
	out << "CONTROL   1      1     1     1     1      1" << endl;
	out << "PRINT     1      0     0     0     0      3"  << endl << endl << endl;
	out << "EXCHANGE  0" << endl;
	out << "SCF       4.0" << endl;
	out << "* XANES     4.0" << endl;
	out << "* FMS        3.94559  0" << endl;
	out << "RPATH     8" << endl << endl;
	out << "POTENTIALS" << endl;
	out << "*    ipot\tZ\telement" << endl;

}

/**
 * reads the periodic table from a file and populate in atom_list vector
 */
vector<string> readPeriodicTable()	{
	ifstream infile("periodic-table.txt");
	string atom;
	vector<string> atom_list;
	if (infile.is_open()){
		while( getline( infile, atom ) )	{
			atom_list.push_back(atom);
		}
	}
	else
		cout << "could not open the file";
	return atom_list;
}
/**
 * loads the R values from "R-values.txt" and populate a map
 */
map<string,vector<float> > loadRValues()	{
	ifstream infile("R-values.txt");
	string line;

	map< string,vector<float> > Rmap;

	if (infile.is_open()){
		while( getline( infile, line ) )	{
			std::istringstream iss (line);

			string atomNotation;
			float val;
			vector<float> Rvals;
			iss >> atomNotation;

			while ( !iss.eof() )	{
				iss >> val;
				//cout << " " << val;
				Rvals.push_back(val);
			}

			Rmap[atomNotation] = Rvals;
		}
	}
	else
		cout << "could not open the file";

	return Rmap;
}

/**
 * returns the R value for the input atom.
 * if the atom has more than one R values, it gives user the option to
 * choose the appropriate value.
 *
 * @param Rmap [in]		map which stores the name of atom and corresponding R value(s).
 * @param atom [in]		name of the atom ("Mn") to look up the R value
 * @return				R value for the given atom
 *
 */
float getRValues(map<string,vector<float> > Rmap, string atom)	{
	vector<float> vals;
	float RVal;
	map<string,vector<float> > :: iterator Rmap_iter;
	Rmap_iter = Rmap.find(atom);
	if ( Rmap_iter != Rmap.end() )	{
		vals = Rmap_iter->second;
	}
	if (vals.size() > 1)	{
		cout << atom << " has more than one value. Choose one from below " << endl;
		cout << "index \t Atom \t R-val" << endl;
		int i = 1;
		vector<float>::iterator val_iter;
		vector<float> RVals;
		for (val_iter = vals.begin();val_iter != vals.end(); val_iter++)	{
			cout << i++ << "\t";
			cout << atom << *val_iter++ << "\t";
			RVals.push_back(*val_iter);
			cout << *val_iter << endl;

		}
		int idx;
		cout << "Enter the index : ";
		cin >> idx;
	
		RVal = RVals[idx-1];

	} else
		RVal = *vals.begin();

	return RVal;
}

/**
 * reads the input points from a file and return as a vector of points
 */
vector<point> readFile(char*  fp)	{

	ifstream infile(fp);
	string line;
	vector<point> vp;

	int lineno = 1;
	string no_of_atoms;


	if (infile.is_open()){
		while( getline( infile, line ) )	{
			if (lineno == 8)
				no_of_atoms = line;
			if (lineno > 8) {
				std::istringstream iss (line);

				string x_coord, y_coord, z_coord;

				int type;
				string atom;

				iss >> x_coord;
				iss >> y_coord;
				iss >> z_coord;

				iss >> type;
				iss >> atom;

				if (atom_list.find(type) == atom_list.end())	// add unique atom values to atom_list map
					atom_list[type]=atom;

				point p(convertToDouble(x_coord),convertToDouble(y_coord),convertToDouble(z_coord), type);
				vp.push_back(p);
			}
			lineno++;
		}
		cout << "atom_list size: " << atom_list.size() << endl;
		infile.close();
	}
	else
		cout << "could not open the file";

	cout << "no of atoms = " << no_of_atoms << endl;

	return vp;
}

/**
 * prints the list of points to screen or file
 *
 * @param out [in]		output stream to write to
 * @param vp [in]		vector of points to be printed
 */
void printPointlist(ostream &out, const vector<point>& vp)	{
	unsigned int size = vp.size();
	out.precision(numeric_limits<double>::digits10 + 1);
	for (unsigned int i=0; i<size;i++)	{
		out << vp[i].x << "  " << vp[i].y << "  " << vp[i].z << "  " << vp[i].type << endl;

	}
}

/**
 * structure to store the min, max and their respective positions
 */
struct MinMax	{
	double xmin;
	int xmin_idx;
	double xmax;
	int xmax_idx;
	double ymin;
	int ymin_idx;
	double ymax;
	int ymax_idx;
	double zmin;
	int zmin_idx;
	double zmax;
	int zmax_idx;
};

/**
 * find the min max from a vector of points
 */
void find_MinMax(const vector<point>& vp, MinMax& m)	{
	unsigned int size = vp.size();
	point p;

	p = vp[0];

	m.xmin = p.x; m.xmin_idx = 0;
	m.xmax = p.x; m.xmax_idx = 0;
	m.ymin = p.y; m.ymin_idx = 0;
	m.ymax = p.y; m.ymax_idx = 0;
	m.zmin = p.z; m.zmin_idx = 0;
	m.zmax = p.z; m.zmax_idx = 0;

	for (unsigned int i=1; i<size;i++)	{
		p = vp[i];

		if (p.x < m.xmin)	{
			m.xmin = p.x;
			m.xmin_idx = i;
		}
		if (p.y < m.ymin)	{
			m.ymin = p.y;
			m.ymin_idx = i;
		}
		if (p.z < m.zmin)	{
			m.zmin = p.z;
			m.zmin_idx = i;
		}
		if (p.x > m.xmax)	{
			m.xmax = p.x;
			m.xmax_idx = i;
		}
		if (p.y > m.ymax)	{
			m.ymax = p.y;
			m.ymax_idx = i;
		}
		if (p.z > m.zmax)	{
			m.zmax = p.z;
			m.zmax_idx = i;
		}
	}
}

/**
 * splits the box into 8 small boxes.
 *
 * @param vp [in]			vector of points in a box
 * @param gp_box [out] 		vector of vector of points, vector(group) of small vector points which are classified by the procedure
 * @param mm [in]			object of MinMax
 * @param dx [in]			double, span of box in x direction
 * @param dy [in]			double, span of box in y direction
 * @param dz [in]			double, span of box in z direction
 */

void Find_Small_boxes(const vector<point>& vp, vector<vector<point> >& gp_box, const MinMax& mm, const double& dx, const double& dy, const double& dz)	{
	unsigned int size = vp.size();
	cout << "Initial Vector size: " << size << endl;
	point p;
	// 8 vectors for small boxes of point
	vector <point> b0;
	vector <point> b1;
	vector <point> b2;
	vector <point> b3;
	vector <point> b4;
	vector <point> b5;
	vector <point> b6;
	vector <point> b7;

	for (unsigned int i=0; i<size;i++)	{
		p = vp[i];

		if ( (p.x >= mm.xmin ) && (p.y >= mm.ymin) && (p.z >= mm.zmin) && (p.x <= mm.xmin+dx/2) && (p.y <= mm.ymin+dy/2) && (p.z <= mm.zmin+dz/2) )
			b0.push_back(p);
		else if ( (p.x > mm.xmin+dx/2 ) && (p.y >= mm.ymin) && (p.z >= mm.zmin) && (p.x <= mm.xmax) && (p.y <= mm.ymin+dy/2) && (p.z <=mm.zmin+dz/2) )
			b1.push_back(p);
		else if ( (p.x >= mm.xmin ) && (p.y > mm.ymin+dy/2) && (p.z >= mm.zmin) && (p.x <= mm.xmin+dx/2) && (p.y <= mm.ymax) && (p.z <=mm.zmin+dz/2) )
			b2.push_back(p);
		else if ( (p.x > mm.xmin+dx/2 ) && (p.y > mm.ymin+dy/2) && (p.z >= mm.zmin) && (p.x <= mm.xmax) && (p.y <= mm.ymax) && (p.z <=mm.zmin+dz/2) )
			b3.push_back(p);
		else if ( (p.x >= mm.xmin ) && (p.y >= mm.ymin) && (p.z > mm.zmin+dz/2) && (p.x <= mm.xmin+dx/2) && (p.y <= mm.ymin+dy/2) && (p.z <= mm.zmax) )
			b4.push_back(p);
		else if ( (p.x > mm.xmin+dx/2 ) && (p.y >= mm.ymin) && (p.z > mm.zmin+dz/2) && (p.x <= mm.xmax) && (p.y <= mm.ymin+dy/2) && (p.z <=mm.zmax) )
			b5.push_back(p);
		else if ( (p.x >= mm.xmin ) && (p.y > mm.ymin+dy/2) && (p.z > mm.zmin+dz/2) && (p.x <= mm.xmin+dx/2) && (p.y <= mm.ymax) && (p.z <=mm.zmax) )
			b6.push_back(p);
		else if ( (p.x > mm.xmin+dx/2 ) && (p.y > mm.ymin+dy/2) && (p.z > mm.zmin+dz/2) && (p.x <= mm.xmax) && (p.y <= mm.ymax) && (p.z <=mm.zmax) )
			b7.push_back(p);

	}

	gp_box.push_back(b0);
	gp_box.push_back(b1);
	gp_box.push_back(b2);
	gp_box.push_back(b3);
	gp_box.push_back(b4);
	gp_box.push_back(b5);
	gp_box.push_back(b6);
	gp_box.push_back(b7);

}
/**
 * prints point object to output stream either file or screen
 */
void printPoint(ostream &out, point p)	{

	out << p.x << " " << p.y << " " << p.z << endl;
}

/**
 * print ANNpoint object to output stream either file or screen
 *
 * @param out [in]	object of outputstream
 * @param p [in]	object of ANNpoint
 */
void printPt(ostream &out, ANNpoint p)			// print ANNpoint
{
	out << "(" << p[0];
	for (int i = 1; i < dim; i++) {
		out << ", " << p[i];
	}
	out << ")\n";
}
/**
 * method used in debugging, prints each point in each small box
 *
 * @param gp_box [in] 		vector of vector of points, vector(group) of small vector points (small box)
 */
void printGPbox(vector<vector<point> >& gp_box)	{

	ofstream out;			// for debugging
	out.open("group-box.txt");	// for debugging

	vector<vector<point> > :: iterator box_iter;
	vector<point> :: iterator point_iter;

	for(box_iter=gp_box.begin(); box_iter != gp_box.end(); box_iter++)	{
		out << "-----------box size--------" << box_iter->size() << endl;


		//box_iter=gp_box.begin()+1;
		for (point_iter=box_iter->begin(); point_iter != box_iter->end(); point_iter++)	{
				out << point_iter->x << " " << point_iter->y << " " << point_iter->z << endl;
		}
	}
	out.close();
}

/**
 * This method extrapolates the points to make a big Box
 *
 * @param gp_box [in] 		vector of vector of points, vector(group) of small vector points (small box)
 * @param bigBox [out] 		vector of points, big box after extrapolation
 * @param dx [in]			double, span of box in x direction
 * @param dy [in]			double, span of box in y direction
 * @param dz [in]			double, span of box in z direction
 */
void fillBigBox(vector<vector<point> >& gp_box, vector<point>& bigBox, const double dx, const double dy, const double dz)	{

	//filling the layers

	double x_incr = dx;
	double y_incr = dy;
	double z_incr = dz;

	double i_x_add = -dx;	// first block adder
	double i_y_add = -dy;
	double i_z_add = -dz;

	double ii_x_add = 0;	// Second block adder
	double ii_y_add = -dy;
	double ii_z_add = -dz;

	double iii_x_add = -dx;	// third block adder
	double iii_y_add = 0;
	double iii_z_add = -dz;

	double iv_x_add = 0;	// fourth block adder
	double iv_y_add = 0;
	double iv_z_add = -dz;

	double v_x_add = -dx;	// fifth block adder
	double v_y_add = -dy;
	double v_z_add = 0;

	double vi_x_add = 0;	// sixth block adder
	double vi_y_add = -dy;
	double vi_z_add = 0;

	double vii_x_add = -dx;	// seventh block adder
	double vii_y_add = 0;
	double vii_z_add = 0;

	double viii_x_add = 0;	// eighth block adder
	double viii_y_add = 0;
	double viii_z_add = 0;

	vector< vector<point> >::iterator box_iter;
	vector<point>::iterator point_iter;
	point p;

	for (unsigned int i=0; i<2;i++)	{				// iterates in z-direction

		for (unsigned int j=0; j < 2; j++)	{		// iterates the first two rows of blocks twice

			for (unsigned int k=0; k < 2; k++)	{	// iterates twice block 1 and 2 in x-direction

				box_iter = gp_box.begin() + 7; // get the 7th box
				// for each point in the 7th box add the x_add, y_add and z_add respectively
				for (point_iter = box_iter->begin(); point_iter != box_iter->end(); point_iter++)	{
					p.x = point_iter->x + i_x_add;
					p.y = point_iter->y + i_y_add;
					p.z = point_iter->z + i_z_add;
					p.type = point_iter->type;

					bigBox.push_back(p);

				}

				box_iter = gp_box.begin() + 6; // get the 6th box
				// for each point in the 6th box add the x_add, y_add and z_add respectively
				for (point_iter = box_iter->begin(); point_iter != box_iter->end(); point_iter++)	{
					p.x = point_iter->x + ii_x_add;
					p.y = point_iter->y + ii_y_add;
					p.z = point_iter->z + ii_z_add;
					p.type = point_iter->type;

					bigBox.push_back(p);

				}
				i_x_add += x_incr;	// add x_incr to first block
				ii_x_add += x_incr;	// add x_incr to second block
			}

			i_x_add -= 2*x_incr;	// get x_add for 1st block back
			ii_x_add -= 2*x_incr;	// get x_add for 1st block back

			for (unsigned int k=0; k < 2; k++)	{	// iterates twice block 3 and 4 in x-direction
				//i_x_add += x_incr;

				box_iter = gp_box.begin() + 5; 		// get the 5th box
				// for each point in the 5th box add the x_add, y_add and z_add respectively
				for (point_iter = box_iter->begin(); point_iter != box_iter->end(); point_iter++)	{
					p.x = point_iter->x + iii_x_add;
					p.y = point_iter->y + iii_y_add;
					p.z = point_iter->z + iii_z_add;
					p.type = point_iter->type;

					bigBox.push_back(p);

				}

				box_iter = gp_box.begin() + 4; // get the 4th box
				// for each point in the 4th box add the x_add, y_add and z_add respectively
				for (point_iter = box_iter->begin(); point_iter != box_iter->end(); point_iter++)	{
					p.x = point_iter->x + iv_x_add;
					p.y = point_iter->y + iv_y_add;
					p.z = point_iter->z + iv_z_add;
					p.type = point_iter->type;

					bigBox.push_back(p);

				}
				iii_x_add += x_incr;	// add x_incr to third block
				iv_x_add += x_incr;		// add x_incr to fourth block
			}
			iii_x_add -= 2*x_incr;		// get x_add for 3rd block back
			iv_x_add -= 2*x_incr;		// get x_add for 4th block back

			i_y_add += y_incr;
			ii_y_add += y_incr;
			iii_y_add += y_incr;
			iv_y_add += y_incr;
		}

		i_y_add -= 2*y_incr;		// get y_add for 1st block back
		ii_y_add -= 2*y_incr;		// get y_add for 2nd block back
		iii_y_add -= 2*y_incr;		// get y_add for 3rd block back
		iv_y_add -= 2*y_incr;		// get y_add for 4th block back


		for (unsigned int j=0; j < 2; j++)	{		// iterates the first two rows of blocks twice

			for (unsigned int k=0; k < 2; k++)	{	// iterates twice block 5 and 6 in x-direction

				box_iter = gp_box.begin() + 3; // get the 3rd box
				// for each point in the 3rd box add the x_add, y_add and z_add respectively
				for (point_iter = box_iter->begin(); point_iter != box_iter->end(); point_iter++)	{
					p.x = point_iter->x + v_x_add;
					p.y = point_iter->y + v_y_add;
					p.z = point_iter->z + v_z_add;
					p.type = point_iter->type;

					bigBox.push_back(p);

				}

				box_iter = gp_box.begin() + 2; // get the 2nd box
				// for each point in the 2nd box add the x_add, y_add and z_add respectively
				for (point_iter = box_iter->begin(); point_iter != box_iter->end(); point_iter++)	{
					p.x = point_iter->x + vi_x_add;
					p.y = point_iter->y + vi_y_add;
					p.z = point_iter->z + vi_z_add;
					p.type = point_iter->type;

					bigBox.push_back(p);

				}
				v_x_add += x_incr;	// add x_incr to first block
				vi_x_add += x_incr;	// add x_incr to second block
			}

			v_x_add -= 2*x_incr;	// get x_add for 1st block back
			vi_x_add -= 2*x_incr;	// get x_add for 1st block back

			for (unsigned int k=0; k < 2; k++)	{	// iterates twice block 7 and 8 in x-direction

				box_iter = gp_box.begin() + 1; 		// get the 1st box
				// for each point in the 1st box add the x_add, y_add and z_add respectively
				for (point_iter = box_iter->begin(); point_iter != box_iter->end(); point_iter++)	{
					p.x = point_iter->x + vii_x_add;
					p.y = point_iter->y + vii_y_add;
					p.z = point_iter->z + vii_z_add;
					p.type = point_iter->type;

					bigBox.push_back(p);

				}

				box_iter = gp_box.begin() + 0; // get the 0th box
				// for each point in the 0th box add the x_add, y_add and z_add respectively
				for (point_iter = box_iter->begin(); point_iter != box_iter->end(); point_iter++)	{
					p.x = point_iter->x + viii_x_add;
					p.y = point_iter->y + viii_y_add;
					p.z = point_iter->z + viii_z_add;
					p.type = point_iter->type;

					bigBox.push_back(p);

				}
				vii_x_add += x_incr;	// add x_incr to third block
				viii_x_add += x_incr;	// add x_incr to fourth block
			}
		vii_x_add -= 2*x_incr;			// get x_add for 7th block back
			viii_x_add -= 2*x_incr;		// get x_add for 8th block back

			v_y_add += y_incr;			// add y_incr to 5th block
			vi_y_add += y_incr;			// add y_incr to 6th block
			vii_y_add += y_incr;		// add y_incr to 7th block
			viii_y_add += y_incr;		// add y_incr to 8th block
		}

		v_y_add -= 2*y_incr;			// get y_add for 5th block back
		vi_y_add -= 2*y_incr;			// get y_add for 6th block back
		vii_y_add -= 2*y_incr;			// get y_add for 7th block back
		viii_y_add -= 2*y_incr;			// get y_add for 8th block back

		i_z_add += z_incr;
		ii_z_add += z_incr;
		iii_z_add += z_incr;
		iv_z_add += z_incr;
		v_z_add += z_incr;
		vi_z_add += z_incr;
		vii_z_add += z_incr;
		viii_z_add += z_incr;
	}

}

/**
 * displays the list of atoms present in the input
 */
void displayAtomlist()	{
	map<int,string>::iterator atom_list_itr;
	cout << "List of available atoms" << endl;
	for(atom_list_itr = atom_list.begin();atom_list_itr != atom_list.end(); atom_list_itr ++)	{
		cout << atom_list_itr->first << "  " << atom_list_itr->second << endl;
	}
}

/**
 * structure to store the min, max and their respective positions
 */
struct queryPoint	{
	int lineNo;
	point qpoint;
};

/**
 * constructs querylist from choice
 *
 * @param point_list [in]		list of points
 * @param choice [in]			choice of atom for query
 * @param query_list [out]		constructed querylist from point_list using choice
 */
void constructQuerylist(vector<point>& point_list, int choice, vector<queryPoint> & query_list )	{
	vector<point>:: iterator point_list_itr;

	int i = 1;
	for(point_list_itr=point_list.begin(); point_list_itr != point_list.end(); point_list_itr++)	{

		if (point_list_itr->type == choice)	{
			queryPoint qp;
			qp.lineNo = point_list_itr - point_list.begin()+1;
			qp.qpoint = *point_list_itr;
			query_list.push_back(qp);
		}
		i++;
	}

}

/**
 * copies point from point object to ANNpoint object
 */
void readPt(point pt, ANNpoint p)
{
	p[0] = pt.x;
	p[1] = pt.y;
	p[2] = pt.z;

}

/**
 * returns the atomic number of an atom
 *
 * @param notation [in]			Atomic notation of the atom
 * @param table [in]			periodic table
 * @return						atomic number
 */
int lookupAtomicNumber(string notation, vector<string> table)	{

	vector<string>::iterator table_iter = find(table.begin(), table.end(), notation);
	if ( table_iter != table.end() )
		return distance(table.begin(), table_iter) + 1;
	else
		return -1;
}

/**
 * construct a unique list of atoms
 *
 * @param aset [in,out]			unique set of atom type
 * @param s [in, out] 			string stream
 * @param type [in]				atom type
 * @param number [in]			atomic number of the atom
 * @param notation [in]			atomic notation of the atom
 */
void addToListOfAtoms(set<int>& aset, stringstream& s, int type, int number, string name )	{
	if ( aset.find(type) == aset.end()	)	{		// if type is not found, then insert it
		aset.insert(type);
		
		s << "\t" << type << "\t" << number << "\t" << name << endl;
	}

}

/**
 * copies bigBox points to ANNpointArray datapoints
 */
void copyBigBoxPoints(vector<point> &ext_box, ANNpointArray &dataPts)	{
	vector<point>::iterator bigBox_iter;

		// copy all the points in bigBox vector to dataPts
		
		int i=0;
		for (bigBox_iter = ext_box.begin(); bigBox_iter != ext_box.end(); bigBox_iter++ )	{
			readPt(*bigBox_iter,dataPts[i]);		// copying the pt to dataPts

			i++;
		}
}

/**
 * prints the make script for one input
 *
 * @param out [in]			output stream to write in
 * @param out_file [in]		name of the output file
 * @param lineNo [in]		line number of the query atom
 */
void printMakeScript(ostream &out, string out_file, int lineNo)	{
	out << "mkdir cluster" << lineNo << endl;
	out << "cp " << out_file << " cluster" << lineNo << "/feff.inp" << endl;
	out << "cd cluster" << lineNo << endl;
	out << "run_feff 4" << endl;
	out << "cp chi.data ../chi.dat" << lineNo << endl;
	out << "cd .." << endl << endl;
}

/**
 * prints the distances of the atoms within the given radius,
 * computes the coordination number and bond valence and
 * prints the script to run for feff
 *
 * @param querylist [in]	list that contains all the query atoms
 * @param kdTree [in]		built kd-tree from all the points form extended box
 * @param ext_box [in]		list of points from extended box
 * @param radius [in]		input radius
 * @param coordType [in]	type of the coordination atom
 * @param coordDist	[in]	maximum distance for the coordination atom
 * @param Rval [in]			R value to compute bond valence
 */
void printDistances(vector<queryPoint> &querylist, ANNkd_tree* kdTree, vector<point> &ext_box, double radius, int coordType, float coordDist, float Rval)	{

	vector<queryPoint>::iterator query_iter;

	ANNpoint			queryPt;				// query point
	ANNidxArray			nnIdx;					// near neighbor indices
	ANNdistArray		dists;					// near neighbor distances

	vector<string> periodicTable = readPeriodicTable();

	queryPt = annAllocPt(dim);					// allocate query point with 'dim' dimensions

	nnIdx = NULL;
	dists = NULL;

	ofstream out, outbv, outmk;	// outbv for writing bond Valence into file, outmk for make feffscript
	stringstream st;		// string stream to output all the bond valence and coordination number for all the query atoms
	st.precision(5);		// set the precision of stream t to 5

	// creates a file for feffscript
	string str(outputFilenameDist);
	string outmkFilename = str + "_feffscript.make";
	char outmkFile[40];
	strcpy(outmkFile, outmkFilename.c_str());
	outmk.open(outmkFile);

	if (Rval > 0 )	{		// when bond valence needs to be computed
		string str(outputFilenameDist);
		string str2(outputFilenameVal);
		string outbvFilename = str +"_" +str2;
		char outbvfile[40];							// char array variable to store output file name
		strcpy(outbvfile, outbvFilename.c_str() ); 	// convert the string to char array
		outbv.open(outbvfile);						// opens the bond valence output file

		st << "row#\tAtom\tCN\tBV" << "\n";
	}
	// get each point from the query list and finds the neighborhood
	for (query_iter = querylist.begin(); query_iter != querylist.end(); query_iter++ )	{
		readPt(query_iter->qpoint,queryPt);	// reads the query point from query_iter and copy to queryPt
		string out_file = outputFilenameDist + static_cast<ostringstream*>( &(ostringstream() << query_iter->lineNo ) )->str() + ".inp"; // convert lineno to string and concatenate

		

		char outfile[30];					// char array variable to store output file name
		strcpy(outfile, out_file.c_str());	// convert the string to char array
		out.open(outfile);					// opens the output file

		
		//initially k is set to zero, as this function will return the number of points "NeighborCount" present in the radius*radius
		int NeighborCount = kdTree->annkFRSearch(						// search
										queryPt,						// query point
								radius * radius,						// square of the radius for search
											k=0,						// number of near neighbors
											nnIdx,						// nearest neighbors (returned)
											dists,						// distance (returned)
											eps);						// error bound

		k = NeighborCount;
		cout  << "within the given radius : " << NeighborCount << endl;
		
		nnIdx = new ANNidx[k];						// allocate near neigh indices
		dists = new ANNdist[k];						// allocate near neighbor dists

		// once again call the annkFRSearch to get the k values in nnIdx and dists
		NeighborCount = kdTree->annkFRSearch(						// search
												queryPt,						// query point
										radius * radius,						// radius
													k,							// number of near neighbors
													nnIdx,						// nearest neighbors (returned)
													dists,						// distance (returned)
													eps);						// error bound

		set<int> atom_set;			// set which contains the unique atoms in the neighborhood

		printHeader(out);			// prints the header information in the outputfile
		stringstream s;				// string stream to output the unique atoms in the neighborhood
		stringstream t;				// string stream to output all the atoms with distances in the neighborhood
		t.precision(5);				// set the precision of stream t to 5

		float bondVal = 0.0;		// bond valence of coordination type
		int coordNum = 0;			// coordination number of coordination type

		string queryAtomName="";	// name of the query atom

		for (int i = 0; i < k; i++) {			// for each atom in the querylist, prints the neighbors and its distances
			dists[i] = sqrt(dists[i]);			// unsquare distance
			point pt = ext_box[nnIdx[i]];		// nnIdx[i] points to neighbor in ext_box

			// populate the list of unique atoms, present in the neighborhood
			string atom_name = "";
			map <int,string>:: iterator atom_iter = atom_list.find(pt.type);
			if ( atom_iter != atom_list.end() )
				atom_name = atom_iter->second;

			if (Rval > 0)	{	// compute coordination number and valence
				if ( (pt.type == coordType) && (dists[i] <= coordDist) )	{	// we have got the coordination atom
					if ( ( atom_name == "Ti" ) && (dists[i] >= 1.84 )	)	{	// second case for 'Ti' otherwise the R value of Ti = 1.89
						Rval = 1.815;
					}
					float bv = exp( ( Rval - dists[i] ) / 0.37 );
					bondVal += bv;
					coordNum++ ;
				}
			}

			// lookup the atomic number from the periodic table
			int atomic_number = lookupAtomicNumber(atom_name, periodicTable);

			if (i==0)	{	// the query atom
				queryAtomName = atom_name;
				s << "\t0" << "\t" << atomic_number << "\t" << atom_name << endl;
				t << " \t" << fixed << setw(15) << pt.x - query_iter->qpoint.x;
				t << fixed << setw(15) << pt.y - query_iter->qpoint.y;
				t << fixed << setw(15)<< pt.z - query_iter->qpoint.z;
				t << setw(15) << "0";
				t << setw(4) << atom_name ;
				t << setw(15) << dists[i] << "\n";
			}
			else	{	// neighbors of the query atom

				addToListOfAtoms(atom_set,s,pt.type,atomic_number,atom_name);	// add unique atoms to a list
				t << " \t" << fixed << setw(15) << pt.x - query_iter->qpoint.x;
				t << fixed << setw(15) << pt.y - query_iter->qpoint.y;
				t << fixed << setw(15) << pt.z - query_iter->qpoint.z;
				t << setw(15) << pt.type;
				t << setw(4) << atom_name;
				t << setw(15) << dists[i] << "\n";
			}
		}

		if (Rval > 0)	{
			st << query_iter->lineNo << "\t" << queryAtomName << "\t" << coordNum << "\t" <<bondVal << "\n";

		}

		atom_set.clear();							// empty the set
		out << s.str() << endl << endl << endl;		// prints the unique atom list
		out << "ATOMS" << endl;
		out << "*\t" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << setw(15) << "ipot" << setw(4) << "tag" << setw(15)<< "Distance" << endl;
		out << t.str();

		out.flush();
		delete [] nnIdx;							// clean things up
		delete [] dists;
		out.close();								// close the output file

		// writes the make script
		printMakeScript(outmk, out_file, query_iter->lineNo);
	}

	if (Rval > 0)	{	// when bond valence needs to be computed
		outbv << st.str();
		outbv.flush();
		outbv.close();
	}
	outmk.flush();
	outmk.close();
	cout << "Read query Points " << endl;
}

/**
 * takes the input from the user
 */
void takesInput (int &atomChoice, double &radius, int &coordAtomChoice, float &coordDist, float &RVal)	{
	char continueValence; 			// 'y/n' choice to continue calculating valence
	bool correctEntry = false;		// check whether entered the right choice

	do {		// iterates till get the right entry
		cout << "Enter the atom number to search for neighborhood : ";
			cin >> atomChoice;

		map <int,string>:: iterator atom_iter = atom_list.find(atomChoice);
		if ( atom_iter != atom_list.end() )
			correctEntry = true;
		else
			cout << "Wrong Entry" << endl << "Again ";
	} while (!correctEntry);

	cout << "Enter the radius to search: ";
	cin >> radius;

	cout << "Would you like to calculate bond valence (y)es/n(o)? ";
	cin >> continueValence;

	if (tolower(continueValence) == 'y')	{

		map<int,string>::iterator atom_list_itr = atom_list.find(atomChoice);
		string atomNotation = atom_list_itr->second;

		map<string, vector <float> > RMap = loadRValues();		// loads the R values from the file

		correctEntry = false;									// check whether entered the right choice
		do {		// iterates till get the right entry
			cout << "Enter the coordination atom from the index: ";
			cin >> coordAtomChoice;

			map <int,string>:: iterator atom_iter = atom_list.find(coordAtomChoice);
			if ( atom_iter != atom_list.end() )
				correctEntry = true;
			else
				cout << "Wrong Entry" << endl << "Again ";
		} while (!correctEntry);

		cout << "Enter the coordination distance: ";
		cin >> coordDist;

		RVal = getRValues(RMap,atomNotation);
		cout << "Rval = " << RVal << endl;

	}
}

int main(int argc, char **argv)	{

	vector<point> box;							// vector of all points calling it a box

	ANNpointArray		dataPts;				// data points
	ANNkd_tree*			kdTree;					// search structure

	MinMax mm;									// object contain Min Max information of points

	vector<vector<point> > group_box;			// to store group of small boxes

	getArgs(argc, argv);						// read command-line arguments

	box = readFile(inputFilename);				// read the input file and populate the box


	find_MinMax(box, mm);	// find the min & max of the points in x,y & z.

	ofstream mout;			// for debugging
	mout.open("min-max.txt");	// for debugging

		mout << "xmax : " << mm.xmax << " idx = "<< mm.xmax_idx << endl;
		printPoint(mout,box[mm.xmax_idx]);
		mout << "xmin : " << mm.xmin << " idx = "<< mm.xmin_idx << endl;
		printPoint(mout,box[mm.xmin_idx]);
		mout << "ymax : " << mm.ymax << " idx = "<< mm.ymax_idx << endl;
		printPoint(mout,box[mm.ymax_idx]);
		mout << "ymin : " << mm.ymin << " idx = "<< mm.ymin_idx << endl;
		printPoint(mout,box[mm.ymin_idx]);
		mout << "zmax : " << mm.zmax << " idx = "<< mm.zmax_idx << endl;
		printPoint(mout,box[mm.zmax_idx]);
		mout << "zmin : " << mm.zmin << " idx = "<< mm.zmin_idx << endl;
		printPoint(mout,box[mm.zmin_idx]);
	mout.close();

	// computes the dx, dy and dz
	double dx = mm.xmax - mm.xmin;
	double dy = mm.ymax - mm.ymin;
	double dz = mm.zmax - mm.zmin;

	// partitions the box into 8 small boxes
	Find_Small_boxes(box,group_box,mm,dx,dy,dz);

	//After wrapping the box all around symmetrically, the size of points increases to 8 times the original
	vector<point> ext_box;	// larger box which also contains the points of neighboring boxes as well


	fillBigBox(group_box, ext_box, dx, dy, dz);	// fills the larger box with values


	int bigBoxSize = ext_box.size();			// size of the big box

	displayAtomlist();		// displays the list of atoms present in the input

	int atomChoice;			// atom number for neighborhood search
	double radius;			// radius in which search is made

	int coordAtomChoice = 0;	// coordination atom choice
	float coordDist = 0.0;		// coordination distance

	float RVal = 0.0;			// R Value

	takesInput (atomChoice, radius, coordAtomChoice, coordDist, RVal);		// takes the input from the user

	cout << "RVal = " << RVal << endl;

	vector<queryPoint> querylist;

	constructQuerylist(box,atomChoice,querylist);	// construct querylist from atomChoice

	dataPts = annAllocPts(bigBoxSize, dim);			// allocate array of data points

//	copy all the points in bigBox vector to dataPts
	copyBigBoxPoints(ext_box, dataPts);

	cout << "Read AnnPoints " << endl;

	kdTree = new ANNkd_tree(					// build search structure
						dataPts,				// the data points
						bigBoxSize,				// number of points
						dim);					// dimension of space
	cout << "kdTree built " << endl;

	printDistances(querylist, kdTree, ext_box, radius, coordAtomChoice, coordDist, RVal);	// prints the distances of the atoms within the given radius

	delete kdTree;								// delete the kdTree
	annClose();									// done with ANN


	return 0;
}
/**
 * parses the command line arguments and initializes the variable
 */
void getArgs(int argc, char **argv)
{
	strcpy(inputFilename,"MD-data");
	strcpy(outputFilenameDist,"feff");
	strcpy(outputFilenameVal,"bondVal.txt");

	int i = 1;
	while (i < argc) {							// read arguments
		if (!strcmp(argv[i], "-i")) {			// -i option
			strcpy(inputFilename,argv[++i]);
			
		}
		else if (!strcmp(argv[i], "-d")) {		// -d option
			strcpy(outputFilenameDist,argv[++i]);
			
		}
		else if (!strcmp(argv[i], "-b")) {		// -b option
			strcpy(outputFilenameVal,argv[++i]);
			
		}
		else if (!strcmp(argv[i], "--help")) {	// --help option
			cerr << "Usage:\n\n"
					<< "  CoNBoVal-XAS [-i inputfile] [-d prefix] [-b bondValence]"
					   "\n\n"
					<< "  where:\n"
					<< "    inputfile		name of the input file containing position of atoms (default = \"MD-data\")\n"
					<< "    prefix     		prefix name of the output files that contain distances (default = \"feff\")\n"
					<< "    bondValence     suffix name of the output file that contain bond valence (default = \"_bondVal.txt\")\n"
					<< " Results are stored in prefix and bondValence files.\n"
					<< "\n"
					<< " To run this demo use:\n"
					<< "    CoNBoVal-XAS -i MD-data -d feff -b bond-valence.txt \n"
					<< "\n"
					<< "or \n"
					<< "CoNBoVal-XAS \n";
					exit(0);
		}

		else {									// illegal syntax
			cerr << "Unrecognized option.\n";
			exit(1);
		}
		i++;
	}

}



