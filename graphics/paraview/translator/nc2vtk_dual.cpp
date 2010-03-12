////////////////////////////////////////////////////////////////////////////////
//
// This program translates a newpop netCDF data file to dual-grid legacy, ascii VTK format
// The variables that have time dim are automatically written.
//
// Assume all variables are of interest.
// Assume variable data type is double.
// Assume variable dims are (Time, nCells|nVertices, nVertLevels|nVertLevelsP1, [nTracers])
// Assume no more than 100 vars each for cell and point data
// Does not deal with edge data.
// Doubles converted to floats in .vtk files follow these rules:
//      if a NaN, then become -FLT_MAX.
//      if positive infinity, then become FLT_MAX
//      if negative infinity, then become -FLT_MAX
//      if smaller than a float, then become 0.
//      if not zero, but between -FLT_MIN and FLT_MIN make -FLT_MIN or FLT_MIN
//      if outside of -FLT_MAX and FLT_MAX make -FLT_MAX or FLT_MAX
//
//
// Version 1.3
// Christine Ahrens
// 3/8/2010
//
////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <sstream>
#include "stdlib.h"
#include "vtkCellType.h"
#include "netcdfcpp.h"
#include <string>
#include <cmath>
#include <cfloat>

using namespace std;

#define CHECK_MALLOC(ptr) \
	if (ptr == NULL) { \
		cerr << "malloc failed!\n"; \
		exit(1); \
	} 
#define MAX_VARS 100

int my_isinf_f(float x){
      	if (isinf(x)) {
		return (x < 0.0 ? -1 : 1);
	} else return 0;
}

float convertDouble2ValidFloat(double inputData) {

        // check for NaN
        if (inputData != inputData) {
                cerr << "found NaN!" << endl;
                return -FLT_MAX;
        }

        // check for infinity
        int retval = my_isinf_f((float)inputData);
        if (retval < 0) {
                return -FLT_MAX;
        } else if (retval > 0) {
                return FLT_MAX;
        }

        // check number too small for float
        if (abs(inputData) < 1e-126) return 0.0;

        if ((float)inputData == 0) return 0.0;

        if ((abs(inputData) > 0) && (abs(inputData) < FLT_MIN)) {
                if (inputData < 0) return -FLT_MIN; else return FLT_MIN;
        }

        if (abs(inputData) > FLT_MAX) {
                if (inputData < 0) return -FLT_MAX; else return FLT_MAX;
        }

        return (float)inputData;
}

int main(int argc, char* argv[])
{
	if ((argc < 3) || (argc > 4))  {
		cerr << "Usage: transdual infile.nc outfile.vtk [verticalLevel]" << endl;
		cerr << "Variables with time and vertical level are written out." << endl;
		cerr << "Tracer vars are named tracer1, tracer2, etc." << endl;
		cerr << "If vertical level is not specified, default is 0." << endl;
		cerr << "A series of vtk files will be created, one file for each time step." << endl;
		cerr << "with time prepended to the file name (e.g. 0outfile.vtk)." << endl;
		cerr << "vtk datafile Version 2.0, transdual Version 1.2" << endl;
		exit(1);
	}

		
	NcFile ncFile(argv[1]);

	if (!ncFile.is_valid()) 
	{
		cerr << "Couldn't open file: " << argv[1] << endl;
		exit(1);
	}

	NcDim* nCells = ncFile.get_dim("nCells");
	NcDim* nVertices = ncFile.get_dim("nVertices");
	NcDim* vertexDegree = ncFile.get_dim("vertexDegree");
	NcDim* Time = ncFile.get_dim("Time");
	NcDim* nVertLevels = ncFile.get_dim("nVertLevels");
	NcDim* nTracers = ncFile.get_dim("nTracers");

	// Can't check for this, b/c if not there it crashes program	
	//NcDim* nVertLevelsP1 = ncFile.get_dim("nVertLevelsP1");
	int maxNVertLevels = nVertLevels->size() + 1;
	/*if (nVertLevelsP1 != NULL) {
		maxNVertLevels = nVertLevelsP1->size();
	}
*/
 
	int outputVertLevel = 0;

	if (argc == 4) {
		outputVertLevel = atoi(argv[3]);
        	// cout << "outputVertLevel: " << outputVertLevel << endl;
		if (outputVertLevel > (maxNVertLevels-1)) {
			cerr << "Specified vertical level " << outputVertLevel;
			cerr << " doesn't exist.  The highest level is level ";
			cerr << maxNVertLevels-1 << "." << endl;
			exit(1);
		}
	}


	// figure out what variables to visualize
	NcVar* dualCellVars[MAX_VARS];
	NcVar* dualPointVars[MAX_VARS];
	int dualCellVarIndex = -1;
	int dualPointVarIndex = -1;
	int numDualCells = nVertices->size();
	int numDualPoints = nCells->size()+1;
	
	int numVars = ncFile.num_vars();

	bool tracersExist = false;

	for (int i = 0; i < numVars; i++) {
		NcVar* aVar = ncFile.get_var(i);

		// must have 3 dims 
		// (Time, nCells | nVertices, nVertLevels | nVertLevelsP1)

		int numDims = aVar->num_dims();
		//cout << "Num Dims of var: " << aVar->name() << " is " << numDims << endl;
		
		if ((numDims != 3) && (strcmp(aVar->name(), "tracers"))) {
			continue;
		} else {
			// TODO, check if it is a double
			// assume a double for now

			// check for Time dim 0
			NcToken dim0Name = aVar->get_dim(0)->name();
			if (strcmp(dim0Name, "Time")) 
				continue;

			// check for dim 1 being num vertices or cells 
			bool isVertexData = false;
			bool isCellData = false;
			NcToken dim1Name = aVar->get_dim(1)->name();
			if (!strcmp(dim1Name, "nVertices")) 
				isVertexData = true;
			else if (!strcmp(dim1Name, "nCells")) 
				isCellData = true; 
			else continue;

			// check if dim 2 is nVertLevels or nVertLevelsP1, too
			NcToken dim2Name = aVar->get_dim(2)->name();
			if ((strcmp(dim2Name, "nVertLevels")) 
				&& (strcmp(dim2Name, "nVertLevelsP1"))) {
				continue;
			}

			// check if we have data for the selected level
			if (aVar->get_dim(2)->size()-1 < outputVertLevel) {
				cout << "No data found for level ";
				cout << outputVertLevel << " for variable ";
				cout << aVar->name() << endl;
				continue;
			}

			// Add to cell or point var array
			if (isVertexData) {  // means it is dual cell data
				dualCellVarIndex++;
				if (dualCellVarIndex > MAX_VARS-1) {
					cerr << "Exceeded number of cell vars." << endl;
					exit(1);
				}
				dualCellVars[dualCellVarIndex] = aVar;
				//cout << "Adding var " << aVar->name() << " to dualCellVars" << endl;
			} else if (isCellData) { // means it is dual vertex data
				if (strcmp(aVar->name(), "tracers")) {
					dualPointVarIndex++;
					if (dualPointVarIndex > MAX_VARS-1) {
						cerr << "Exceeded number of point vars." << endl;
						exit(1);
					}
					dualPointVars[dualPointVarIndex] = aVar;
					//cout << "Adding var " << aVar->name() << " to dualPointVars" << endl;
				} else { // case of tracers, add each as "tracer0", "tracer1", etc.
					tracersExist = true;
					int numTracers = aVar->get_dim(3)->size();
					for (int t = 0; t < numTracers; t++) {
						dualPointVarIndex++;
						if (dualPointVarIndex > MAX_VARS-1) {
							cerr << "Exceeded number of point vars." << endl;
							exit(1);
						}
						dualPointVars[dualPointVarIndex] = aVar;
						//cout << "Adding var " << aVar->name() << " to dualPointVars" << endl;
					}
				}	
			}
		}
	}

	// TODO
 	// prompt the user to find out which fields are of interest?
	// for now, assume all are of interest
	
	// get points  (centers of primal-mesh cells)

	// TO DO check malloc return vals.

	double *xCellData = (double*)malloc(nCells->size() 
		* sizeof(double));
	CHECK_MALLOC(xCellData);
	NcVar *xCellVar = ncFile.get_var("xCell");
	xCellVar->get(xCellData, nCells->size());

	double *yCellData = (double*)malloc(nCells->size() 
		* sizeof(double));
	CHECK_MALLOC(yCellData);
	NcVar *yCellVar = ncFile.get_var("yCell");
	yCellVar->get(yCellData, nCells->size());

	double *zCellData = (double*)malloc(nCells->size() 
		* sizeof(double));
	//cout << "ptr for zCellData"  << zCellData << endl;
	CHECK_MALLOC(zCellData);
	NcVar *zCellVar = ncFile.get_var("zCell");
	zCellVar->get(zCellData, nCells->size());

	// get dual-mesh cells

	int *cellsOnVertex = (int *) malloc((nVertices->size()) * vertexDegree->size() * 
		sizeof(int));
	//cout << "ptr for cellsOnVertex"  << cellsOnVertex << endl;
	CHECK_MALLOC(cellsOnVertex);
	NcVar *cellsOnVertexVar = ncFile.get_var("cellsOnVertex");
	//cout << "getting cellsOnVertexVar\n";
	cellsOnVertexVar->get(cellsOnVertex, nVertices->size(), vertexDegree->size());

	// decls for data storage
	double* dualCellVarData;
	double* dualPointVarData;

	// for each variable, allocate space for variables

	//cout << "dualCellVarIndex: " << dualCellVarIndex << endl;

	int varVertLevels = 0;

	dualCellVarData = (double*)malloc((sizeof(double))  * numDualCells);
	CHECK_MALLOC(dualCellVarData);

	dualPointVarData = (double*)malloc((sizeof(double)) * nCells->size());
	CHECK_MALLOC(dualPointVarData);
	

	// for each time step, write a file with the time prepended to filename
	for (int i = 0; i < Time->size(); i++) 
	{

		ostringstream vtkFileName;

		vtkFileName << i << argv[2]; 	

		ofstream vtkFile(vtkFileName.str().c_str(), ios::out);

		if (!vtkFile)
		{
			cerr << "vtk output file could not be opened" <<endl;
			exit(1);
		}


		// write header

		vtkFile << "# vtk DataFile Version 2.0, transdual Version 1.2" << endl;
		vtkFile << "Translated newpop data to dual grid for timestep " << i << " from netCDF by Christine Ahrens" 
			<< endl;
		vtkFile << "ASCII" << endl;
		vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;


		// write points  (the points are the primal-grid cell centers)

		vtkFile << "POINTS " << numDualPoints << " float" << endl;

		// first write a dummy point, because the climate code
		// starts their cell numbering at 1 and VTK starts it at
		// 0
		vtkFile.precision(16);
		vtkFile << (float)0.0 << "\t" << (float)0.0 << "\t" << (float)0.0 
			<< endl;

		for (int j = 0; j < nCells->size(); j++ )
		{
			vtkFile.precision(16);
			if (abs(xCellData[j]) < 1e-126) xCellData[j] = 0;
			if (abs(yCellData[j]) < 1e-126) yCellData[j] = 0;
			if (abs(zCellData[j]) < 1e-126) zCellData[j] = 0;
			vtkFile << (float)xCellData[j] << "\t" << (float)yCellData[j] << "\t"
				<< (float)zCellData[j] << endl;
		}	
		vtkFile << endl;


		// Write dual-mesh cells
		// Dual-mesh cells are triangles with primal-mesh cell 
		// centers as the vertices.
		// The number of dual-mesh cells is the number of vertices in the
		// primal mesh.

		vtkFile << "CELLS " << numDualCells << " " 
			<< numDualCells * (vertexDegree->size() + 1) << endl;	

		// for each dual-mesh cell, write number of points for each
		// and then list the points by number

		for (int j = 0; j < numDualCells ; j++) {

			vtkFile << vertexDegree->size() << "\t" ;

			// since primal vertex(pt) numbers  == dual cell numbers
			// we go through the primal vertices, find the cells around
			// them, and since those primal cell numbers are dual 
			// point numbers,  
			// we can write the cell numbers for the cellsOnVertex
			// and those will be the numbers of the dual vertices (pts).
			
			int* dualCells = cellsOnVertex + (j * vertexDegree->size());

			for (int k = 0; k < vertexDegree->size(); k++) 
			{		
				vtkFile << dualCells[k] << "\t";
			}

			vtkFile << endl;
		}	
		vtkFile << endl;


		// write cell types 
		int cellType;
		if (vertexDegree->size() == 3) 
			cellType = VTK_TRIANGLE;
		else if (vertexDegree->size() == 4)
			cellType = VTK_QUAD;
		else cellType = VTK_POLYGON;

		vtkFile << "CELL_TYPES " << numDualCells << endl;

		for (int j = 0; j < numDualCells; j++)
		{
			vtkFile << cellType << endl;
		}
		vtkFile << endl;

		// Write attributes of dual-mesh cell data (attributes of primal-mesh
		// vertex data)

		// for each var, figure out if it is by cell or by point


		vtkFile.precision(16);

		// If by point, write out point data

		if (dualPointVarIndex >= 0) vtkFile << "POINT_DATA " << numDualPoints << endl;

		int printstep = -1;

		if (i == printstep) cout << "TIME STEP: " << i << endl;

		int tracerNum = 0;

		for (int v = 0; v <= dualPointVarIndex; v++) {
	
			// Read variable number v data for that timestep

			varVertLevels = dualPointVars[v]->get_dim(2)->size();

			bool isTracer = false;

			if (!strcmp(dualPointVars[v]->name(), "tracers")) {
				isTracer = true;				
				dualPointVars[v]->set_cur(i, 0, outputVertLevel, tracerNum);
				dualPointVars[v]->get(dualPointVarData, 1, nCells->size(), 1, 1);
			} else {	
				dualPointVars[v]->set_cur(i, 0, outputVertLevel);
				dualPointVars[v]->get(dualPointVarData, 1, nCells->size(), 1);
			}

			// Write variable number v data for that timestep
			vtkFile << "SCALARS " << dualPointVars[v]->name();
			if (isTracer) vtkFile << tracerNum+1;
			vtkFile << " float 1" <<  endl;
			vtkFile << "LOOKUP_TABLE default" << endl;

			//debugging
			if (i == printstep) cout << "SCALARS " << dualPointVars[v]->name() << " float 1" <<  endl;
			if (i==printstep) {
				for (int z = 0; z < varVertLevels; z++) {
					cout << z << ":" << *(dualPointVarData + z) << endl;
				}
			}

			// get starting point
			double *var_target = dualPointVarData;

			float validData;

			validData = convertDouble2ValidFloat (*var_target);
			
			// write dummy
			vtkFile << validData  << endl;
		
			// write data	
			for (int j = 0; j < nCells->size(); j++) 
			{
				validData = convertDouble2ValidFloat (*var_target);
				vtkFile << validData << endl;
				var_target++;
			}

			if (isTracer) tracerNum++;
		}

		// if by cell, then write out cell data

		if (dualCellVarIndex >= 0) vtkFile << "CELL_DATA " << numDualCells << endl;

		for (int v = 0; v <= dualCellVarIndex; v++) {

			// Read variable number v data for that timestep
			varVertLevels = dualCellVars[v]->get_dim(2)->size();
			dualCellVars[v]->set_cur(i, 0, outputVertLevel);
			dualCellVars[v]->get(dualCellVarData, 1, numDualCells, 1);

			// Write variable number v data for that timestep
			vtkFile << "SCALARS " << dualCellVars[v]->name() << " float 1" <<  endl;
			vtkFile << "LOOKUP_TABLE default" << endl;

			// debugging
			if (i==printstep) cout << "SCALARS " << dualCellVars[v]->name() << " float 1" <<  endl;
			if (i==printstep) {
				for (int z = 0; z < varVertLevels; z++) {
					cout << z << ":" << *(dualCellVarData + z) << endl;
				}
			}

			double *var_target = dualCellVarData; 
			float validData;

			for (int j = 0; j < numDualCells; j++) 
			{
				validData = convertDouble2ValidFloat (*var_target);
				vtkFile << validData << endl;
				var_target++;
			}
		}

	}

	return(0);

}	

