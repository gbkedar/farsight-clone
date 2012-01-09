#ifndef HEATMAPWINDOW_H
#define HEATMAPWINDOW_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <strstream>
#include <string>
#include <math.h>
 
#include <QVTKWidget.h>
#include <QtGui/QAction>
#include <QtGui/QMainWindow>
#include <QtGui/QApplication>
#include <QtGui/QDesktopWidget>
#include <QtGui/QWidget>
#include <QApplication>
#include <QFileDialog>
#include <QFile>
#include <QCoreApplication>
#include <QTextStream>

#include <vtkTable.h>
#include <vtkLookupTable.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkViewTheme.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkGraphLayout.h>
#include <vtkGraphLayoutView.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderedGraphRepresentation.h>
#include <vtkGraphToGlyphs.h>
#include <vtkPolyDataMapper.h>
#include <vtkGraphToPolyData.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkLine.h>
#include <vtkLineWidget2.h>
#include <vtkLineSource.h>
#include <vtkLineRepresentation.h>
#include <vtkIdTypeArray.h>
#include <vtkAbstractArray.h>
#include <vtkAnnotationLink.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPoints.h>
#include <vtkCallbackCommand.h>
#include <vtkFloatArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkRendererCollection.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdTypeArray.h>
#include <vtkCommand.h>
#include <vtkCellArray.h>
#include <vtkPlaneSource.h>
#include <vtkCellPicker.h>
#include <vtkPicker.h>
#include <vtkPointPicker.h>
#include <vtkInteractorStyleImage.h>
#include <vtkExtractSelection.h>
#include <vtkObjectFactory.h>
#include <vtkStringArray.h>
#include <vtkPointSetToLabelHierarchy.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkLabelPlacementMapper.h>
#include <vtkQtLabelRenderStrategy.h>
#include <vtkPointData.h>
#include <vtkScalarBarActor.h>


#include <boost/math/distributions/normal.hpp>

#include <ftkCommon/ftkUtils.h>
#include "ClusClus/clusclus.h"
#include "ObjectSelection.h"
#include "ftkGUI/ColorMap.h"

using namespace std;

class Heatmap : public QMainWindow
{
    Q_OBJECT;

public:
	Heatmap(QWidget * parent = 0);
	~Heatmap();
	void setDataForHeatmap(double** features, int* optimalleaforder1, int* optimalleaforder2,int num_samples, int num_features);
	void setDataForDendrograms(double** treedata1, double** treedata2 = NULL);
	void creatDataForHeatmap(double powCof);
	void setModels(vtkSmartPointer<vtkTable> table = NULL, ObjectSelection * sels = NULL, ObjectSelection * sels2 = NULL);
	void runClusclus();
	void runClus();
	void showGraph();
	void GetSelRowCol(int &r1, int &c1, int &r2, int &c2);
	void SetSelRowCol(int r1, int c1, int r2, int c2);
	void SetInteractStyle();
	void showDendrogram1();
	void showDendrogram2();	

	int              num_samples;
	int              num_features;
	double**         mapdata;
	double**         connect_Data_Tree1;
	double**         connect_Data_Tree2;
	int*             Optimal_Leaf_Order1;
	int*             Optimal_Leaf_Order2;
	
	vector<vector<double > > Processed_Coordinate_Data_Tree1;
	vector<vector<double > > Processed_Coordinate_Data_Tree2;	

	ObjectSelection * Selection;
	ObjectSelection * Selection2;

signals:
	void SelChanged();

protected slots:
	void SetdenSelectedIds1(std::set<long int>& IDs, bool bfirst);
	void GetSelecectedIDs();
	static void SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void SelectionCallbackFunction1(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void SelectionCallbackFunction2(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void SelectionCallbackFunction3(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void HandleKeyPress(vtkObject* caller, long unsigned eventId, void* clientData, void* callData );

private:
	QVTKWidget mainQTRenderWidget;
	vtkSmartPointer<vtkTable > table;
	vtkSmartPointer<vtkPlaneSource> aPlane;
	vtkSmartPointer<vtkFloatArray> cellData;
	vtkSmartPointer<vtkLookupTable> celllut;
	vtkSmartPointer<vtkPolyDataMapper> mapper;
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkGraphLayoutView> view;
	vtkSmartPointer<vtkViewTheme> theme;
	vtkSmartPointer<vtkIdTypeArray>ids1;
	vtkSmartPointer<vtkIdTypeArray>ids2;

	vtkSmartPointer<vtkIdTypeArray> v;
	vtkSmartPointer<vtkMutableUndirectedGraph> graph_Layout;
	vtkSmartPointer<vtkPoints> points;
	vtkSmartPointer<vtkIntArray> vertexColors;
	vtkSmartPointer<vtkLookupTable> vetexlut;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback1;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback2;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback3;
	vtkSmartPointer<vtkCallbackCommand> keyPress;
	vtkSmartPointer<vtkCellPicker> myCellPicker;
	vtkSmartPointer<vtkIdTypeArray> ids;

	vtkSmartPointer<vtkPoints> denpoints1;
	vtkSmartPointer<vtkCellArray> denlines1;
	vtkSmartPointer<vtkPolyData> denlinesPolyData1;
	vtkSmartPointer<vtkPolyDataMapper> denmapper1;
	vtkSmartPointer<vtkActor> denactor1;
	vtkSmartPointer<vtkUnsignedCharArray> dencolors1;

	vtkSmartPointer<vtkPoints> denpoints2;
	vtkSmartPointer<vtkCellArray> denlines2;
	vtkSmartPointer<vtkPolyData> denlinesPolyData2;
	vtkSmartPointer<vtkPolyDataMapper> denmapper2;
	vtkSmartPointer<vtkActor> denactor2;
	vtkSmartPointer<vtkUnsignedCharArray> dencolors2;

	vtkSmartPointer<vtkVariantArray> featureName;

	rgb GetRGBValue(double val);
	void readmustd(double** mustd);
	void scaleData(double** mustd);
	void scaleData();
	void drawPoints1();
	void drawPoints3();
	void setselectedCellIds();
	void computeselectedcells();
	void createDataForDendogram1(double powCof);
	void createDataForDendogram2(double powCof);
	void createDataForDendogram2();
	void reselectIds1(std::set<long int>& selectedIDs, long int id);
	void reselectIds2(std::set<long int>& selectedIDs2, long int id);

	std::map<int, int> indMapFromVertexToInd;
	std::vector<int> indMapFromIndToVertex;

	int     r1;
	int     r2;
	int     c1;
	int     c2;
	int     removeActorflag;
	int     denResetflag1;
	int     denResetflag2;
	int     continueselectnum;
	bool	clusflag;
	bool    continueselect;
	bool    intersectionselect;
	vtkIdType id1;
	vtkIdType id2;
	std::set<long int> interselectedIDs;

	clusclus *cc1;
	clusclus *cc2;
};

#endif