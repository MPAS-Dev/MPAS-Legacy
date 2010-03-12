#ifndef __vtkNewPopReader_h
#define __vtkNewPopReader_h

#define MAX_VARS 100
#define MAX_VAR_NAME 100
 
#include "vtkUnstructuredGridAlgorithm.h"
#include "netcdfcpp.h"
 
class vtkCallbackCommand;
class vtkDataArraySelection;
class vtkFloatArray;
class vtkStdString;
class vtkStringArray;

 
class VTK_IO_EXPORT vtkNewPopReader : public vtkUnstructuredGridAlgorithm 
{
public:
  static vtkNewPopReader *New();
  vtkTypeRevisionMacro(vtkNewPopReader,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify file name of NewPop data file to read.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Get the number of data cells
  vtkGetMacro(NumberOfDualCells, int);

  // Description:
  // Get the number of points
  vtkGetMacro(NumberOfDualPoints, int);

  // Description:
  // Get the number of data variables at the cell centers and points
  vtkGetMacro(NumberOfVariables, int);

  // Description:
  // Get the reader's output
  vtkUnstructuredGrid *GetOutput();
  vtkUnstructuredGrid *GetOutput(int index);

  // Description:
  // The following methods allow selective reading of solutions fields.
  // By default, ALL data fields on the nodes are read, but this can
  // be modified.
  int GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int GetPointArrayStatus(const char* name); 
  void SetPointArrayStatus(const char* name, int status);
  void DisableAllPointArrays();
  void EnableAllPointArrays();

  int GetNumberOfCellArrays();
  const char* GetCellArrayName(int index);
  int GetCellArrayStatus(const char* name); 
  void SetCellArrayStatus(const char* name, int status);
  void DisableAllCellArrays();
  void EnableAllCellArrays();
  void SetVerticalLevel(int level);
  vtkGetVector2Macro(VerticalLevelRange, int);

protected:
  vtkNewPopReader();
  ~vtkNewPopReader();

  char *FileName;			// First field part file giving path
/*  
  int Rank;				// Number of this processor
  int TotalRank;			// Number of processors
*/

//  int NumberOfPieces;			// Number of files in dataset
  vtkIdType NumberOfDualPoints;		// Number of points in grid
  vtkIdType NumberOfDualCells;		// Number of cells in grid
// vtkIdType NumberOfTuples;		// Number of tuples in sub extent

  int NumberOfVariables;		// Number of variables to display
  vtkStdString* VariableName;		// Names of each variable
  int* VariableType;			// Scalar, vector or tensor

  int NumberOfTimeSteps;		// Temporal domain
  double* TimeSteps;			// Times available for request
  double dTime;

  vtkFloatArray** dualCellVarData;	// Actual data arrays
  vtkFloatArray** dualPointVarData;	// Actual data arrays

  // Selected field of interest
  vtkDataArraySelection* PointDataArraySelection;
  vtkDataArraySelection* CellDataArraySelection;

  int VerticalLevelSelected;
  int VerticalLevelRange[2];

  // Observer to modify this object when array selections are modified
  vtkCallbackCommand* SelectionObserver;

  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *);
  int RequestInformation(vtkInformation *, vtkInformationVector **,
                         vtkInformationVector *);

  void LoadGeometryData(int var, double dTime);
  void LoadPointData(int var);
  void LoadCellData(int var);

  static void SelectionCallback(vtkObject* caller, unsigned long eid,
                                void* clientdata, void* calldata);

  bool infoRequested;
  NcDim* nCells;
  NcDim* nVertices;
  NcDim* vertexDegree;
  NcDim* Time;
  NcDim* nVertLevels;
  NcFile* ncFile;
  NcVar* dualCellVars[MAX_VARS];
  NcVar* dualPointVars[MAX_VARS];
  char tracerNames[MAX_VAR_NAME][MAX_VARS];
  int numDualCellVars;
  int numDualPointVars;
  double* primalPointVarData;
  double* primalCellVarData;
  int ReadAndOutputDualGrid();
  int ReadAndOutputVariableData();
  int LoadPointVarData(int variable, double dTime);
  int LoadCellVarData(int variable, double dTime);
  int BuildVarArrays();

private:
  vtkNewPopReader(const vtkNewPopReader&);	// Not implemented.
  void operator=(const vtkNewPopReader&);	// Not implemented.

};

#endif


