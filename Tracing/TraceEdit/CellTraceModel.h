/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/
#ifndef CELLTRACEMODEL_H
#define CELLTRACEMODEL_H

#include <vector>
#include <string>
#include <set>
#include <fstream>
//QT INCLUDES
#include <QtCore>
#include <QtGui>
//vtk table includes
#include "vtkTable.h"
#include <vtkTableToGraph.h>
#include <vtkViewTheme.h>
#include <vtkStringToCategory.h>
#include <vtkGraphLayout.h>
#include <vtkGraphLayoutView.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkGraphToGlyphs.h>
#include <vtkRenderer.h>
#include <vtkFast2DLayoutStrategy.h>
#include <vtkArcParallelEdgeStrategy.h>
#include <vtkPolyDataMapper.h>
#include <vtkEdgeLayout.h>
#include <vtkGraphToPolyData.h>
#include <vtkActor.h>
#include <vtkProperty.h>


#include <vtkAbstractArray.h>
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkAbstractArray.h"
#include "vtkVariantArray.h"

//ftk includes
#include <ftkGUI/ObjectSelection.h>
#include <ftkGUI/GraphWindow.h>
#include "ftkGraphs/kNearestObjects.h"

class QStandardItemModel;
class QItemSelectionModel;
class TraceLine;
class CellTrace;
class CellTraceModel : public QObject
{
	Q_OBJECT

public:
	CellTraceModel();
	CellTraceModel(std::vector<CellTrace*> Cells);
	~CellTraceModel();
	void setCells(std::vector<CellTrace*> Cells);
	vtkSmartPointer<vtkTable> getDataTable();
	void setDataTable(vtkSmartPointer<vtkTable> table);
	vtkSmartPointer<vtkTable> getCellBoundsTable();
	ObjectSelection * GetObjectSelection();
	////////////////////////////////////////////
	ObjectSelection * GetObjectSelectionColumn();
	///////////////////////////////////////////////
	void SelectByRootTrace(std::vector<TraceLine*> roots);
	void SelectByIDs(std::vector<int> IDs);
	std::set<long int> GetSelectedIDs();
	std::vector<CellTrace*> GetSelectedCells();
	unsigned int getCellCount();
	CellTrace * GetCellAt( int i);
	CellTrace * GetCellAtNoSelection( int i);
	void WriteCellCoordsToFile(const char* fileName);
	std::vector<CellTrace*> getCells(std::vector<long> IDs);
	void createCellToCellGraph();
	double average(std::vector< std::pair<unsigned int, double> > ID);
	int AddNewFeatureHeader(std::string NewHeader);
signals:
	void selectionChanged(void);
private:
	std::vector<CellTrace*> Cells;
	std::vector<QString> headers;
	std::vector<QString> AdditionalHeaders;
	void SetupHeaders();
	void SyncModel();
	vtkSmartPointer<vtkTable> DataTable;
	ObjectSelection * Selection;
	//////////////////////////////////////////////
	ObjectSelection * ColumnSelection;
	////////////////////////////////////////////
	GraphWindow * graphVisualize;
};

#endif