////////////////////////////////////////////////////////////////////////////////
//
// This program translates a newpop netCDF data file to dual-grid legacy, ascii VTK format
// The variables that have time dim are automatically written.
// Version Beta
//
// Assume all variables are of interest.
// Assume variable data type is double.
// Assume variable dims are (Time, nCells|nVertices, nVertLevels|nVertLevelsP1, [nTracers])
// Assume no more than 100 vars each for cell and point data
// Does not deal with edge data.
// Only vis up to nVertLevels, not nVertLevelsP1.
// Doubles converted to floats in .vtk files follow these rules:
//	if a NaN, then become -FLT_MAX.
//	if positive infinity, then become FLT_MAX
//	if negative infinity, then become -FLT_MAX
//	if smaller than a float, then become 0.
//	if not zero, but between -FLT_MIN and FLT_MIN make -FLT_MIN or FLT_MIN
//	if outside of -FLT_MAX and FLT_MAX make -FLT_MAX or FLT_MAX
//
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
#define DEFAULT_LAYER_THICKNESS 100000

extern "C" {
	void CartesianToSpherical(float x, float y, float z, float* rho, float* phi, float* theta);
	void SphericalToCartesian(float rho, float phi, float theta, float* x, float* y, float* z);
}	

// START
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

#define BUFF_SIZE 2048

// overwrite outFile if it already exists
int myCopyFile(string *inFile, string *outFile)
{
   char buff[BUFF_SIZE];
   int readBytes = 1;

   ifstream inFileStream(inFile->c_str(), ios::in|ios::binary);
   if(!inFileStream)
   {
     return -1;
   }

   ofstream outFileStream(outFile->c_str(), ios::out|ios::binary);
   if(!outFileStream)
   {
     return -1;
   }

   while(readBytes != 0)
   {
     inFileStream.read((char*)buff, BUFF_SIZE);
     readBytes = inFileStream.gcount();
     outFileStream.write((char*)buff, readBytes);
   }
   return 0;
}

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
		cerr << "Usage: proj3dsphere infile.nc infile.vtk [layer_thickness]" << endl;
		cerr << "Note: layer_thickness defaults to 100000." << endl;
		exit(1);
	}

	//cerr << "FLT_MAX: " << FLT_MAX << " FLT_MIN: " << FLT_MIN << endl;

	NcFile ncFile(argv[1]);

	if (!ncFile.is_valid()) 
	{
		cerr << "Couldn't open file: " << argv[1] << endl;
		exit(1);
	}

	int layerThickness =  DEFAULT_LAYER_THICKNESS;
	if (argc == 4) {
		layerThickness = atoi(argv[3]);
	}

	NcDim* nCells = ncFile.get_dim("nCells");
	NcDim* nVertices = ncFile.get_dim("nVertices");
	NcDim* vertexDegree = ncFile.get_dim("vertexDegree");
	NcDim* Time = ncFile.get_dim("Time");
	NcDim* nVertLevels = ncFile.get_dim("nVertLevels");

	// Can't check for this, b/c if not there it crashes program	
	//NcDim* nVertLevelsP1 = ncFile.get_dim("nVertLevelsP1");
	int maxNVertLevels = nVertLevels->size();

	//cout << "maxNVertLevels: " << maxNVertLevels << endl;

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
			continue; // try the next var
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


	// write a file with the geometry.

	ostringstream geoFileName;

	geoFileName << "geo_" << argv[2]; 	

	ofstream geoFile(geoFileName.str().c_str(), ios::out);

	if (!geoFile)
	{
		cerr << "vtk output file could not be opened" <<endl;
		exit(1);
	}


	// write header

	geoFile << "# vtk DataFile Version 2.0" << endl;
	geoFile << "Translated newpop geometry to dual grid on 3D sphere from netCDF by Christine Ahrens" 
		<< endl;
	geoFile << "ASCII" << endl;
	geoFile << "DATASET UNSTRUCTURED_GRID" << endl;


	// write points  (the points are the primal-grid cell centers)

	geoFile << "POINTS " << numDualPoints*(maxNVertLevels+1) << " float" << endl;

	// write the point at each vertical level, plus one level for last layer

	//cout << "Writing dummy points" << endl;
	for (int levelNum = 0; levelNum < maxNVertLevels+1; levelNum++) {

		// first write a dummy point, because the climate code
		// starts their cell numbering at 1 and VTK starts it at
		// 0
		geoFile.precision(16);
		geoFile << (float)0.0 << "\t" << (float)0.0 << "\t" << (float)0.0 
			<< endl;
	}

	//cout << "Writing points at each level" << endl;
	for (int j = 0; j < nCells->size(); j++ )
	{
		geoFile.precision(16);
		if (abs(xCellData[j]) < 1e-126) xCellData[j] = 0;
		if (abs(yCellData[j]) < 1e-126) yCellData[j] = 0;
		if (abs(zCellData[j]) < 1e-126) zCellData[j] = 0;

		float rho, rholevel, theta, phi, x, y, z;

		CartesianToSpherical((float)xCellData[j], (float)yCellData[j], (float)zCellData[j], 
				&rho, &phi, &theta);

		for (int levelNum = 0; levelNum < maxNVertLevels+1; levelNum++) {

			// decrease rho, so we go down (inwards towards the center of the sphere)
			rholevel = rho - (layerThickness * levelNum);

			SphericalToCartesian(rholevel, phi, theta, &x, &y, &z);

			geoFile << x << "\t" << y << "\t" << z << endl;
		}
	}	

	geoFile << endl;


	// Write dual-mesh cells
	// Dual-mesh cells are triangles with primal-mesh cell 
	// centers as the vertices.
	// The number of dual-mesh cells is the number of vertices in the
	// primal mesh.

	int newDegree = 6;

	geoFile << "CELLS " << numDualCells*maxNVertLevels << " " 
		<< numDualCells * (newDegree + 1) * maxNVertLevels << endl;	

	// for each dual-mesh cell, write number of points for each
	// and then list the points by number

	//cout << "Writing Cells" << endl;

	for (int j = 0; j < numDualCells ; j++) {

		// since primal vertex(pt) numbers  == dual cell numbers
		// we go through the primal vertices, find the cells around
		// them, and since those primal cell numbers are dual 
		// point numbers,  
		// we can write the cell numbers for the cellsOnVertex
		// and those will be the numbers of the dual vertices (pts).

		int* dualCells = cellsOnVertex + (j * vertexDegree->size());

		// for each level, write the prism
		for (int levelNum = 0; levelNum < maxNVertLevels; levelNum++) {
			geoFile << newDegree << "\t" ;
			for (int k = 0; k < vertexDegree->size(); k++) 
			{		
				geoFile << (dualCells[k]*(maxNVertLevels+1)) + levelNum << "\t";
			}
			for (int k = 0; k < vertexDegree->size(); k++) 
			{		
				geoFile << (dualCells[k]*(maxNVertLevels+1)) + levelNum+1 << "\t";
			}
			geoFile << endl;
		}
	}	
	geoFile << endl;


	// write cell types 
	int cellType;

	if (vertexDegree->size() == 3) 
		cellType = VTK_WEDGE;

	geoFile << "CELL_TYPES " << numDualCells*maxNVertLevels << endl;

	//multiply by number of levels
	for (int j = 0; j < numDualCells*maxNVertLevels; j++)
	{
		geoFile << cellType << endl;
	}
	geoFile << endl;


	// release resources
	geoFile.close();
	free(xCellData);
	free(yCellData);
	free(zCellData);
	free(cellsOnVertex);

	// For each timestep, write data for each level

	dualCellVarData = (double*)malloc((sizeof(double))  * numDualCells * maxNVertLevels);
	CHECK_MALLOC(dualCellVarData);

	dualPointVarData = (double*)malloc((sizeof(double)) * nCells->size() * maxNVertLevels);
	CHECK_MALLOC(dualPointVarData);

	// for each timestep, copy the geometry file, since we don't want to
	// have to recompute the points
	for (int i = 0; i < Time->size(); i++) {
		ostringstream vtkFileName;
		vtkFileName << i << argv[2];

		string geoFileNameString = geoFileName.str();
		string vtkFileNameString = vtkFileName.str();
		int copyRetVal = myCopyFile(&geoFileNameString, &vtkFileNameString);
		//cout << "myCopyFile returned: " << copyRetVal << endl;

		ofstream vtkFile(vtkFileName.str().c_str(), ios::out|ios::app);
		if (!vtkFile)
		{
			cerr << "vtk output file could not be opened" <<endl;
			exit(1);
		}

		if (!vtkFile)
		{
			cerr << "vtk output file could not be opened" <<endl;
			exit(1);
		}

		vtkFile.precision(16);

		// If by point, write out point data

		if (dualPointVarIndex >= 0) vtkFile << "POINT_DATA " << numDualPoints*(maxNVertLevels+1) << endl;

		int printstep = -1;

		if (i == printstep) cout << "TIME STEP: " << i << endl;

		int tracerNum = 0;

		for (int v = 0; v <= dualPointVarIndex; v++) {

			// Read variable number v data for that timestep

			varVertLevels = dualPointVars[v]->get_dim(2)->size();

			bool isTracer = false;

			if (!strcmp(dualPointVars[v]->name(), "tracers")) {
				isTracer = true;
			// Uncomment if want to exclude tracers.  
			//	continue;
			}
			// Write variable number v data for that timestep
			vtkFile << "SCALARS " << dualPointVars[v]->name();
			if (isTracer) vtkFile << tracerNum+1;
			vtkFile << " float 1" <<  endl;
			vtkFile << "LOOKUP_TABLE default" << endl;


			if (isTracer) {
				dualPointVars[v]->set_cur(i, 0, 0, tracerNum);
				dualPointVars[v]->get(dualPointVarData, 1, nCells->size(), maxNVertLevels, 1);
			} else {
				dualPointVars[v]->set_cur(i, 0, 0);
				dualPointVars[v]->get(dualPointVarData, 1, nCells->size(), maxNVertLevels);
			}


			float defaultPointVal = 0.0;

			//write dummy

			double *var_target = dualPointVarData;
			float validData;

			for (int levelNum = 0; levelNum < maxNVertLevels; levelNum++) {

				validData = convertDouble2ValidFloat (*var_target);

				// write dummy
				vtkFile << validData << endl;
	
				var_target++;
			}

			// write highest level dummy point
			vtkFile << validData << endl;

			var_target = dualPointVarData;

			for (int j = 0; j < nCells->size(); j++) {

				// write data for one point lowest level to highest
				for (int levelNum = 0; levelNum < maxNVertLevels; levelNum++) {

					validData = convertDouble2ValidFloat (*var_target);
					vtkFile << validData << endl;
					var_target++;
				}
			
				// for last layer of dual points, repeat last level's values
				// Need Mark's input on this one
				vtkFile << validData << endl;
			}

			if (isTracer) tracerNum++;
		}

		// if by cell, then write out cell data

		if (dualCellVarIndex >= 0) vtkFile << "CELL_DATA " << numDualCells*maxNVertLevels << endl;

		for (int v = 0; v <= dualCellVarIndex; v++) {

			// Write variable number v data for that timestep
			vtkFile << "SCALARS " << dualCellVars[v]->name() << " float 1" <<  endl;
			vtkFile << "LOOKUP_TABLE default" << endl;

			// Read variable number v data for that timestep and level
			dualCellVars[v]->set_cur(i, 0, 0);
			dualCellVars[v]->get(dualCellVarData, 1, numDualCells, maxNVertLevels);

			for (int j = 0; j < numDualCells; j++) {

				double *var_target = dualCellVarData;

				for (int levelNum = 0; levelNum < maxNVertLevels; levelNum++)
				{
					float validData = convertDouble2ValidFloat (*var_target);
					vtkFile << validData << endl;

					var_target++;
				}
			}
		}
	}
}
