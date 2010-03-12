
////////////////////////////////////////////////////////////////////////////////
//
// This program reads a newpop netCDF data file to allow paraview to
// display a dual-grid.
// The variables that have time dim are available to paraview.
//
// Assume all variables are of interest if they have dims (Time, nCells|nVertices, nVertLevels, [nTracers])
// Converts variable data type from double to float.
// Assume no more than 100 vars each for cell and point data
// Displays tracer vars as tracer1, tracer2, etc.
// Does not deal with edge data.
//
// Christine Ahrens
// 3/12/2010
// Version 1.1
//
////////////////////////////////////////////////////////////////////////////////

#include "vtkNewPopReader.h"
#include "vtkSmartPointer.h"
#include "vtkCallbackCommand.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkDataArraySelection.h"
#include "vtkDataObject.h"
#include "vtkErrorCode.h"
#include "vtkFloatArray.h"
#include "vtkStringArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkInformationDoubleVectorKey.h"
#include "vtkIntArray.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTableExtentTranslator.h"
#include "vtkToolkits.h"
#include "vtkUnstructuredGrid.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include "stdlib.h"
#include "netcdfcpp.h"
#include <string>
#include <cmath>
#include <cfloat>

using namespace std;

#define CHECK_MALLOC(ptr) \
	if (ptr == NULL) { \
		cerr << "malloc failed!\n"; \
		return(0); \
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

vtkCxxRevisionMacro(vtkNewPopReader, "$Revision$");
vtkStandardNewMacro(vtkNewPopReader);

//----------------------------------------------------------------------------
// Constructor for vtkNewPopReader
//----------------------------------------------------------------------------

vtkNewPopReader::vtkNewPopReader()
{
	// Debugging
	this->DebugOn();
	vtkDebugMacro(<< "Starting to create vtkNewPopReader..." << endl);

	this->SetNumberOfInputPorts(0);
	this->SetNumberOfOutputPorts(1);

	this->FileName = NULL;
	this->ncFile = NULL;

	this->infoRequested = false;

	this->PointDataArraySelection = vtkDataArraySelection::New();
	this->CellDataArraySelection = vtkDataArraySelection::New();

	this->NumberOfDualPoints = 0;
	this->NumberOfDualCells = 0;
	this->NumberOfVariables = 0;

	this->primalPointVarData = NULL;
	this->primalCellVarData = NULL;

	this->dualCellVarData = NULL;
	this->dualPointVarData = NULL;

	// put in defaults
	this->VerticalLevelRange[0] = 0;
	this->VerticalLevelRange[1] = 1;
	this->VerticalLevelSelected = 0;

	// Setup selection callback to modify this object when array selection changes
	this->SelectionObserver = vtkCallbackCommand::New();
	this->SelectionObserver->SetCallback(&vtkNewPopReader::SelectionCallback);
	this->SelectionObserver->SetClientData(this);
	this->CellDataArraySelection->AddObserver(vtkCommand::ModifiedEvent,
			this->SelectionObserver);
	this->PointDataArraySelection->AddObserver(vtkCommand::ModifiedEvent,
			this->SelectionObserver);

	vtkDebugMacro(<< "MAX_VARS:" << MAX_VARS << endl);

	vtkDebugMacro(<< "Created vtkNewPopReader" << endl);
	
}  

//----------------------------------------------------------------------------
// Destructor for NewPop Reader
//----------------------------------------------------------------------------

vtkNewPopReader::~vtkNewPopReader()
{
	vtkDebugMacro(<< "Destructing vtkNewPopReader..." << endl);
	if (this->FileName)
	{
		delete [] this->FileName;
	}

	if (this->ncFile)
	{
		delete ncFile;
	}

	vtkDebugMacro(<< "Destructing cell var array..." << endl);
	if (this->dualCellVarData) {
		for (int i = 0; i < this->numDualCellVars; i++) {
			if (this->dualCellVarData[i] != NULL) {
				this->dualCellVarData[i]->Delete();
			}
		}
		delete [] this->dualCellVarData;
	}

	vtkDebugMacro(<< "Destructing point var array..." << endl);
	if (this->dualPointVarData) {
		for (int i = 0; i < this->numDualPointVars; i++) {
			if (this->dualPointVarData[i] != NULL) {
				this->dualPointVarData[i]->Delete();
			}
		}
		delete [] this->dualPointVarData;
	}

	vtkDebugMacro(<< "Destructing other stuff..." << endl);
	if (this->primalCellVarData) free (this->primalCellVarData);	
	if (this->primalPointVarData) free (this->primalPointVarData);	
	if (this->TimeSteps) delete [] this->TimeSteps;
	if (this->PointDataArraySelection) this->PointDataArraySelection->Delete();
	if (this->CellDataArraySelection) this->CellDataArraySelection->Delete();
	if (this->SelectionObserver) this->SelectionObserver->Delete();

	vtkDebugMacro(<< "Destructed vtkNewPopReader" << endl);
}

//----------------------------------------------------------------------------
// Verify that the file exists, get dimension sizes and variables
//----------------------------------------------------------------------------

int vtkNewPopReader::RequestInformation(
		vtkInformation *reqInfo,
		vtkInformationVector **inVector,
		vtkInformationVector *outVector)
{ 
	vtkDebugMacro(<< "In vtkNewPopReader::RequestInformation" << endl);

	if (!this->Superclass::RequestInformation(reqInfo, inVector, outVector))
		return 0;

	// Verify that file exists
	if ( !this->FileName ) {
		vtkErrorMacro("No filename specified");
		return 0;
	}

	vtkDebugMacro(<< "In vtkNewPopReader::RequestInformation read filename okay" << endl);

	// Get ParaView information and output pointers
	vtkInformation* outInfo = outVector->GetInformationObject(0);
	vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
			outInfo->Get(vtkDataObject::DATA_OBJECT()));

	// RequestInformation() is called for every Modified() event which means
	// when more variable data is selected it will be called again

	if (!this->infoRequested) {
		this->infoRequested = true;

		vtkDebugMacro(<< "FileName: " << this->FileName << endl);
		ncFile = new NcFile(this->FileName);

		if (!ncFile->is_valid()) 
		{
			vtkErrorMacro( << "Couldn't open file: " << this->FileName << endl);
			return 0;
		}

		vtkDebugMacro(<< "In vtkNewPopReader::RequestInformation read file okay" << endl);

		this->nCells = ncFile->get_dim("nCells");
		this->nVertices = ncFile->get_dim("nVertices");
		this->vertexDegree = ncFile->get_dim("vertexDegree");
		this->Time = ncFile->get_dim("Time");
		this->nVertLevels = ncFile->get_dim("nVertLevels");

		vtkDebugMacro(<< "In vtkNewPopReader::RequestInformation setting VerticalLevelRange" << endl);

		this->VerticalLevelRange[0] = 0;
		if (this->nVertLevels != NULL) {
			this->VerticalLevelRange[1] = this->nVertLevels->size();
		} else {
			this->VerticalLevelRange[1] = 1;
		}

		if (!BuildVarArrays()) return 0;

		// Allocate the data arrays which will hold the NetCDF var data

		this->primalCellVarData = (double*)malloc(sizeof(double)* 
			this->nCells->size());
		CHECK_MALLOC(this->primalCellVarData);
		this->primalPointVarData = (double*) malloc(sizeof(double)* 
			this->nVertices->size());
		CHECK_MALLOC(this->primalPointVarData);

		// Allocate the ParaView data arrays which will hold the variables
		this->dualPointVarData = new vtkFloatArray*[this->numDualPointVars];
		for (int i = 0; i < this->numDualPointVars; i++) {
			this->dualPointVarData[i] = NULL;
		}
		this->dualCellVarData = new vtkFloatArray*[this->numDualCellVars];
		for (int i = 0; i < this->numDualCellVars; i++) {
			this->dualCellVarData[i] = NULL;
		}

		// Start with no data loaded into ParaView
		DisableAllPointArrays();
		DisableAllCellArrays();

		this->NumberOfDualCells = this->nVertices->size();
		this->NumberOfDualPoints = this->nCells->size() + 1;

		// Collect temporal information
		this->NumberOfTimeSteps = this->Time->size();
		this->TimeSteps = NULL;

		// At this time, NewPop doesn't have fine-grained time value, just
		// the number of the step, so that is what I store here for TimeSteps.
		this->TimeSteps = new double[this->NumberOfTimeSteps];
		for (int step = 0; step < this->NumberOfTimeSteps; step++)
			this->TimeSteps[step] = (double) step;

		// Tell the pipeline what steps are available
		outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
				this->TimeSteps, this->NumberOfTimeSteps);

		double tRange[2];
		tRange[0] = this->TimeSteps[0];
		tRange[1] = this->TimeSteps[this->NumberOfTimeSteps-1];
		outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
				tRange,
				2);

	}

	return 1;
}

int vtkNewPopReader::BuildVarArrays() 
{
	vtkDebugMacro(<< "In vtkNewPopReader::BuildVarArrays" << endl);

	// figure out what variables to visualize -
	int dualCellVarIndex = -1;
	int dualPointVarIndex = -1;

	int numVars = this->ncFile->num_vars();

	bool tracersExist = false;

	for (int i = 0; i < numVars; i++) {
		NcVar* aVar = this->ncFile->get_var(i);

		// must have 3 dims (Time, nCells|nVertices, nVertLevels)

		int numDims = aVar->num_dims();
		vtkDebugMacro( << "Num Dims of var: " << aVar->name() << " is " << numDims << endl);
		if ((numDims != 3) && (strcmp(aVar->name(), "tracers"))) {
			continue; // try the next var
		} else {
			// TODO, check if it is a double
			// assume a double for now

			// check for Time dim
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
                        if ((aVar->get_dim(2)->size()-1 < this->VerticalLevelSelected)
				|| (this->VerticalLevelSelected < 0)) {
                                //cout << "No data found for level ";
                                //cout << outputVertLevel << " for variable ";
                                //cout << aVar->name() << endl;
                                continue;
                        }

                        // Add to cell or point var array
                        if (isVertexData) {  // means it is dual cell data
                                dualCellVarIndex++;
				if (dualCellVarIndex > (MAX_VARS-1)) {
					vtkDebugMacro(<< "too many dual cell vars!" << endl);
					exit(0);
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
							vtkDebugMacro(<< "too many point vars!" << endl);
                                                        exit(1);
                                                }
                                                dualPointVars[dualPointVarIndex] = aVar;
                                                //cout << "Adding var " << aVar->name() << " to dualPointVars" << endl;
						ostringstream tracerName;
						tracerName << "tracer" << t+1;
						strcpy(tracerNames[t], tracerName.str().c_str());
                                        }
                                }
                        }
		}
	}

	this->numDualCellVars = dualCellVarIndex+1;
	this->numDualPointVars = dualPointVarIndex+1;

	vtkDebugMacro( << "numDualCellVars: " << this->numDualCellVars << " numDualPointVars: " << this->numDualPointVars << endl);

	int tracerNum = 0;

	for (int var = 0; var < this->numDualPointVars; var++) {
		if (!strcmp(this->dualPointVars[var]->name(), "tracers") ) {
			this->PointDataArraySelection->
				EnableArray((const char*)(this->tracerNames[tracerNum]));
			vtkDebugMacro(<< "Adding point var: " << tracerNames[tracerNum] << endl);
			tracerNum++;
		} else {
			this->PointDataArraySelection->
				EnableArray((const char*)(this->dualPointVars[var]->name()));
				vtkDebugMacro(<< "Adding point var: " << this->dualPointVars[var]->name() << endl);
		}
	}

	for (int var = 0; var < this->numDualCellVars; var++) {
		vtkDebugMacro(<< "Adding cell var: " << this->dualCellVars[var]->name() << endl);
		this->CellDataArraySelection->
			EnableArray((const char*)(this->dualCellVars[var]->name()));
	}

	vtkDebugMacro(<< "Leaving vtkNewPopReader::BuildVarArrays" << endl);

	return(1);
}


//----------------------------------------------------------------------------
// Data is read into a vtkUnstructuredGrid
//----------------------------------------------------------------------------

int vtkNewPopReader::RequestData(
			vtkInformation *vtkNotUsed(reqInfo),
			vtkInformationVector **vtkNotUsed(inVector),
			vtkInformationVector *outVector)
{
	vtkDebugMacro(<< "In vtkNewPopReader::RequestData" << endl);

	// get the info object
	vtkInformation *outInfo = outVector->GetInformationObject(0);

	// Output will be an ImageData
	vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
			outInfo->Get(vtkDataObject::DATA_OBJECT()));

	// Output the unstructured grid from the netCDF file
	if (!ReadAndOutputDualGrid()) return 0;

	// Collect the time step requested
	double* requestedTimeSteps = NULL;
	int numRequestedTimeSteps = 0;
	vtkInformationDoubleVectorKey* timeKey =
		static_cast<vtkInformationDoubleVectorKey*>
		(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS());
	if (outInfo->Has(timeKey)) {
		numRequestedTimeSteps = outInfo->Length(timeKey);
		requestedTimeSteps = outInfo->Get(timeKey);
	}

	// print out how many steps are requested, just for my information
	vtkDebugMacro( << "Num Time steps requested: " << numRequestedTimeSteps << endl);

	// At this time, it seems to only get one timestep of info, why?

	this->dTime = requestedTimeSteps[0];
	double dTimeTemp = this->dTime;
	output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEPS(), &dTimeTemp, 1);
	this->dTime = dTimeTemp;

	// Examine each variable to see if it is selected
	for (int var = 0; var < this->numDualPointVars; var++) {

		// Is this variable requested
		if (this->PointDataArraySelection->GetArraySetting(var)) {
			vtkDebugMacro( << "Loading Point Variable: " << var << endl);
			LoadPointVarData(var, dTime);
			output->GetPointData()->AddArray(this->dualPointVarData[var]);

		}
	}

	for (int var = 0; var < this->numDualCellVars; var++) {
		if (this->CellDataArraySelection->GetArraySetting(var)) {
			vtkDebugMacro( << "Loading Cell Variable: " << this->dualCellVars[var]->name() << endl);
			LoadCellVarData(var, dTime);
			output->GetCellData()->AddArray(this->dualCellVarData[var]);

		}
	}

	vtkDebugMacro( << "Returning from RequestData" << endl);
	return 1;
}

int vtkNewPopReader::ReadAndOutputDualGrid() 
{

	vtkDebugMacro(<< "In vtkNewPopReader::ReadAndOutputDualGrid" << endl);

	// read points  (centers of primal-mesh cells)

	double *xCellData = (double*)malloc(nCells->size()
			* sizeof(double));
	CHECK_MALLOC(xCellData);
	NcVar *xCellVar = ncFile->get_var("xCell");
	xCellVar->get(xCellData, nCells->size());

	double *yCellData = (double*)malloc(nCells->size()
			* sizeof(double));
	CHECK_MALLOC(yCellData);
	NcVar *yCellVar = ncFile->get_var("yCell");
	yCellVar->get(yCellData, nCells->size());

	double *zCellData = (double*)malloc(nCells->size()
			* sizeof(double));
	//cout << "ptr for zCellData"  << zCellData << endl;
	CHECK_MALLOC(zCellData);
	NcVar *zCellVar = ncFile->get_var("zCell");
	zCellVar->get(zCellData, nCells->size());

	// read dual-mesh cells  (triangles formed by primal-mesh cell centers)

	// cellsOnVertex is a 2D array of cell numbers
	// of dimensions (nVertices X vertexDegree)
	int *cellsOnVertex = (int *) malloc((nVertices->size()) * 
		this->vertexDegree->size() * sizeof(int));
	//cout << "ptr for cellsOnVertex"  << cellsOnVertex << endl;
	CHECK_MALLOC(cellsOnVertex);
	NcVar *cellsOnVertexVar = ncFile->get_var("cellsOnVertex");
	//cout << "getting cellsOnVertexVar\n";
	cellsOnVertexVar->get(cellsOnVertex, nVertices->size(), 
		this->vertexDegree->size());


	vtkUnstructuredGrid* output = GetOutput();

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->Allocate(this->NumberOfDualPoints, this->NumberOfDualPoints);

	// first write a dummy point, because the climate code
	// starts their cell numbering at 1 and VTK starts it at
	// 0
	points->InsertNextPoint(0, 0, 0);

	for (int i = 0; i < nCells->size(); i++) {	
		points->InsertNextPoint(xCellData[i], yCellData[i], zCellData[i]);
	}	

	output->SetPoints(points);
	
	free(xCellData);
	free(yCellData);
	free(zCellData);

	// Write dual-mesh cells
	// Dual-mesh cells are triangles with primal-mesh cell 
	// centers as the vertices.
	// The number of dual-mesh cells is the number of vertices in the
	// primal mesh.

	// for each dual-mesh cell

	// write cell types
	int cellType;
	if (vertexDegree->size() == 3)
		cellType = VTK_TRIANGLE;
	else if (vertexDegree->size() == 4)
		cellType = VTK_QUAD;
	else cellType = VTK_POLYGON;

	output->Allocate(this->NumberOfDualCells, this->NumberOfDualCells);

	vtkIdType* polygon = new vtkIdType[this->vertexDegree->size()];
	for (int i = 0; i < nVertices->size(); i++) {

		// since primal vertex(pt) numbers  == dual cell numbers
		// we go through the primal vertices, find the cells around
		// them, and since those primal cell numbers are dual pt numbers
		// we can write the cell numbers for the cellsOnVertex
		// and those will be the numbers of the dual vertices (pts).

		int* dualMeshCells = cellsOnVertex + (i * this->vertexDegree->size());
		for (int j = 0; j < vertexDegree->size(); j++) {
			polygon[j] = dualMeshCells[j];
		}

		// InsertNextCell(type, number of points, array of points)
		output->InsertNextCell(cellType, vertexDegree->size(), polygon);
	}	
	free(cellsOnVertex);
	free(polygon);

	vtkDebugMacro(<< "Leaving vtkNewPopReader::ReadAndOutputDualGrid" << endl);

	return(1);
}

int vtkNewPopReader::LoadPointVarData(int variableIndex, double dTimeStep)
{	
	
	vtkDebugMacro(<< "In vtkNewPopReader::LoadPointVarData" << endl);

	NcVar* ncVar = this->dualPointVars[variableIndex];

	vtkDebugMacro(<< "got ncVar in vtkNewPopReader::LoadPointVarData" << endl);
	if (ncVar == NULL) {
		cerr << "Can't find data for variable " << variableIndex << endl;
		return 0;
	}

	// Allocate data array for this variable

	bool isTracer = false;
	int tracerNum = 0;

	// if it is a tracer var
	if (!strcmp(dualPointVars[variableIndex]->name(), "tracers")) {
		isTracer = true;

		// find its tracer number
		int firstTracerIndex = 0;
		for (int v = 0; v < numDualPointVars; v++) {
			if (!strcmp(dualPointVars[v]->name(), "tracers"))  {
				firstTracerIndex = v;
				break;
			}
		}
		tracerNum = variableIndex - firstTracerIndex;
	}

	vtkDebugMacro(<< "isTracer: " << isTracer << " and tracerNum: " << tracerNum << endl);

	if (this->dualPointVarData[variableIndex] == NULL) {
		vtkDebugMacro(<< "allocating data array in vtkNewPopReader::LoadPointVarData" << endl);
		this->dualPointVarData[variableIndex] = vtkFloatArray::New();
		if (isTracer) {
			this->dualPointVarData[variableIndex]->SetName(tracerNames[tracerNum]);
			vtkDebugMacro(<< "set name to : " << (tracerNames[tracerNum]) << endl);
		} else {
			this->dualPointVarData[variableIndex]->SetName(dualPointVars[variableIndex]->name());
		}
		this->dualPointVarData[variableIndex]->SetNumberOfTuples(this->NumberOfDualPoints);
		this->dualPointVarData[variableIndex]->SetNumberOfComponents(1);
	}

	vtkDebugMacro(<< "getting pointer in vtkNewPopReader::LoadPointVarData" << endl);
	float* dataBlock = this->dualPointVarData[variableIndex]->GetPointer(0);

	// allocate for doubles, will convert to vtkFloatArray 
	int timestep = dTimeStep;

	if (isTracer) {
		dualPointVars[variableIndex]->set_cur(timestep, 0, this->VerticalLevelSelected, tracerNum);
		dualPointVars[variableIndex]->get(primalCellVarData, 1, nCells->size(), 1, 1);
	} else {
		dualPointVars[variableIndex]->set_cur(timestep, 0, this->VerticalLevelSelected);
		dualPointVars[variableIndex]->get(primalCellVarData, 1, nCells->size(), 1);
	}

	vtkDebugMacro(<< "got point data in vtkNewPopReader::LoadPointVarData" << endl);

	double *var_target = this->primalCellVarData;

	// write dummy, just make it the same as the first elt, so we stay in range.
	float validData;
        validData = convertDouble2ValidFloat (*var_target);

	dataBlock[0] = validData;

	for (int j = 0; j < nCells->size(); j++) 
	{
        	validData = convertDouble2ValidFloat (*var_target);
		dataBlock[j+1] =  validData;
		var_target ++;
	}

	vtkDebugMacro(<< "converted and stored point data in vtkNewPopReader::LoadPointVarData" << endl);
	return (1);
}


int vtkNewPopReader::LoadCellVarData(int variableIndex, double dTimeStep)
{	
	vtkDebugMacro(<< "In vtkNewPopReader::LoadCellVarData" << endl);

	NcVar* ncVar = dualCellVars[variableIndex];

	if (ncVar == NULL) {
		cerr << "Can't find data for variable index:" << variableIndex << endl;
		return 0;
	}

	// Allocate data array for this variable
	if (this->dualCellVarData[variableIndex] == NULL) {
		this->dualCellVarData[variableIndex] = vtkFloatArray::New();
		vtkDebugMacro( << "Allocated cell var index: " << dualCellVars[variableIndex]->name() << endl);
		this->dualCellVarData[variableIndex]->SetName(dualCellVars[variableIndex]->name());
		this->dualCellVarData[variableIndex]->SetNumberOfTuples(this->NumberOfDualCells);
		this->dualCellVarData[variableIndex]->SetNumberOfComponents(1);
	}
	float* dataBlock = this->dualCellVarData[variableIndex]->GetPointer(0);

	
	int timestep = dTimeStep;

	ncVar->set_cur(timestep, 0, this->VerticalLevelSelected);

	ncVar->get(this->primalPointVarData, 1, nVertices->size(), 1);

	vtkDebugMacro( << "Got data for cell var: " << dualCellVars[variableIndex]->name() << endl);

	double *var_target = this->primalPointVarData;

	for (int j = 0; j < nVertices->size(); j++) 
	{
       		float validData = convertDouble2ValidFloat (*var_target);
		dataBlock[j] =  validData;
		var_target++;
	}

	vtkDebugMacro( << "Converted and stored data for cell var: " << dualCellVars[variableIndex]->name() << endl);

	return(1);
}


//----------------------------------------------------------------------------
void vtkNewPopReader::SelectionCallback(vtkObject*, unsigned long eventid,
		void* clientdata, void* calldata)
{
	static_cast<vtkNewPopReader*>(clientdata)->Modified();
}

//----------------------------------------------------------------------------
vtkUnstructuredGrid* vtkNewPopReader::GetOutput()
{
	return this->GetOutput(0);
}

//----------------------------------------------------------------------------
vtkUnstructuredGrid* vtkNewPopReader::GetOutput(int idx)
{
	if (idx)
	{
		return NULL;
	}
	else
	{
		return vtkUnstructuredGrid::SafeDownCast( this->GetOutputDataObject(idx) );
	}
}

//----------------------------------------------------------------------------
int vtkNewPopReader::GetNumberOfPointArrays()
{   
	return this->PointDataArraySelection->GetNumberOfArrays();
}   

//----------------------------------------------------------------------------
int vtkNewPopReader::GetNumberOfCellArrays()
{   
	return this->CellDataArraySelection->GetNumberOfArrays();
}
   
//----------------------------------------------------------------------------
void vtkNewPopReader::EnableAllPointArrays()
{     
	this->PointDataArraySelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
void vtkNewPopReader::DisableAllPointArrays()
{   
	this->PointDataArraySelection->DisableAllArrays();
}   

//----------------------------------------------------------------------------
void vtkNewPopReader::EnableAllCellArrays()
{
	this->CellDataArraySelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
void vtkNewPopReader::DisableAllCellArrays()
{   
	this->CellDataArraySelection->DisableAllArrays();
}   


//----------------------------------------------------------------------------
const char* vtkNewPopReader::GetPointArrayName(int index)
{   

        // if it is a tracer var
        if (!strcmp(dualPointVars[index]->name(), "tracers")) {
                // find its tracer number
                int firstTracerIndex = 0;
                for (int v = 0; v < numDualPointVars; v++) {
                        if (!strcmp(dualPointVars[v]->name(), "tracers"))  {
                                firstTracerIndex = v;
                                break;
                        }
                }
                int tracerNum = index - firstTracerIndex;
		vtkDebugMacro( << "GetPointArrayName: " << tracerNames[tracerNum] << endl);
		return (const char*) tracerNames[tracerNum];
	} else {
		vtkDebugMacro( << "GetPointArrayName: " << this->dualPointVars[index]->name() << endl);
		return (const char*)(this->dualPointVars[index]->name());
	}

} 

//----------------------------------------------------------------------------
int vtkNewPopReader::GetPointArrayStatus(const char* name)
{ 
	return this->PointDataArraySelection->ArrayIsEnabled(name);
} 

//----------------------------------------------------------------------------
void vtkNewPopReader::SetPointArrayStatus(const char* name, int status)
{           
	if (status) 
		this->PointDataArraySelection->EnableArray(name);
	else
		this->PointDataArraySelection->DisableArray(name);
}   

//----------------------------------------------------------------------------
const char* vtkNewPopReader::GetCellArrayName(int index)
{   
	return (const char*)(this->dualCellVars[index]->name());
} 

//----------------------------------------------------------------------------
int vtkNewPopReader::GetCellArrayStatus(const char* name)
{ 
	return this->CellDataArraySelection->ArrayIsEnabled(name);
} 

//----------------------------------------------------------------------------
void vtkNewPopReader::SetCellArrayStatus(const char* name, int status)
{           
	if (status) 
		this->CellDataArraySelection->EnableArray(name);
	else
		this->CellDataArraySelection->DisableArray(name);
}   


void vtkNewPopReader::SetVerticalLevel(int level)
{
        this->VerticalLevelSelected = level;
        vtkDebugMacro( << "Set VerticalLevelSelected to: " << level << endl);

        vtkDebugMacro( << "infoRequested?: " << infoRequested << endl);

	if (!this->infoRequested) { return; }
	
	// Examine each variable to see if it is selected
	for (int var = 0; var < this->numDualPointVars; var++) {
		// Is this variable requested
		if (this->PointDataArraySelection->GetArraySetting(var)) {
			vtkDebugMacro( << "Loading Point Variable: " << this->dualPointVars[var]->name() << endl);
			LoadPointVarData(var, dTime);
		}
	}

	for (int var = 0; var < this->numDualCellVars; var++) {
		if (this->CellDataArraySelection->GetArraySetting(var)) {
			vtkDebugMacro( << "Loading Cell Variable: " << this->dualCellVars[var]->name() << endl);
			LoadCellVarData(var, dTime);
		}
	}

	this->PointDataArraySelection->Modified();
	this->CellDataArraySelection->Modified();
}

void vtkNewPopReader::PrintSelf(ostream& os, vtkIndent indent)
{
	os << "NewPop Reader PrintSelf" << endl;
}

