/*
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/

#ifndef NUCLEUS_EDITOR_H
#define NUCLEUS_EDITOR_H

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//QT Includes:
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#include <QtGui/QToolBar>
#include <QtGui/QProgressBar>
#include <QtGui/QDialog>
#include <QtGui/QVBoxLayout>
#include <QtGui/QRadioButton>
#include <QtCore/QFileInfo>
#include <QtCore/QThread>
#include <QtCore/QSettings>
#include <QtCore/QSignalMapper>
#include <QtGui/QGridLayout>
#include <QtGui/QLineEdit>
#include <QtGui/QCheckBox>
#include <QtGui/QLabel>
#include <QtGui/QComboBox>
#include <QtGui/QPushButton>

#include "ProjectFilenamesDialog.h"
#include "ExclusionDialog.h"
#include "PreferencesDialog.h"
#include "ftkProjectProcessor.h"
#include "ftkProjectFiles.h"
#include "ActiveLearningDialog.h"
#include "GalleryDialog.h"

#ifdef USE_TRACKING
#include "KymoGraphView.h" 
#include "TrackingDialog.h"
#include "Image3DView.h"
#include "CellTrackerLib/MultiFrameCellTracker.h"
#endif

#ifdef USE_Clusclus
#include "ClusClus/HeatmapWindow.h"
#endif

//Farsight Includes:
#include "NuclearSegmentation/ftkNuclearSegmentation.h"
#include "NuclearSegmentation/CytoplasmSegmentation/CytoplasmSegmentation.h"
#include "NuclearSegmentation/Nuclear_Association/ftkNuclearAssociationRules.h"
#include "ftkFeatures/ftkLabelImageToFeatures.h"
#include "ftkCommon/ftkUtils.h"
#include "ftkGUI/TrainingDialog.h"
#include "ftkGUI/PatternAnalysisWizard.h"
#include "ftkGUI/TableWindow.h"
#include "ftkGUI/PlotWindow.h"
#include "ftkGUI/ImageBrowser5D.h"
#include "ftkGUI/HistoWindow.h"
#include "ftkGUI/LabelImageViewQT.h"
#include "ftkGUI/PreprocessDialog.h"
#include "ftkGraphs/kNearestObjects.h"
#include "ftkSpectralUnmixing/ftkSpectralUnmixing.h"

//VTK includes:
#include "vtkQtTableView.h"

//ITK Preprocessing includes:
#include "ftkPreprocessDialog.h"
#include "ftkPreprocess.h"

class ParamsFileDialog;
class ProcessThread;


class NucleusEditor : public QMainWindow
{
    Q_OBJECT;

public:

	NucleusEditor(QWidget * parent = 0, Qt::WindowFlags flags = 0);
	~NucleusEditor();
	void loadModelFromFile(std::string file_name);
	void getColorInfo(int numChann,std::vector<std::string> *channelNames,std::vector<unsigned char> *channelColors);
	void getColor(int numChann,std::vector<unsigned char> *channelColors);
	void loadTableSeries(QString fileName);

protected:
	void closeEvent(QCloseEvent *event);

protected slots:
	void setMouseStatus(int,int,int,int,list<int>);
	void EnableModels() {modelsMenu->setEnabled(true);};
	//Loading:
	void askLoadImage(void);
	void loadImage(QString fileName);
	void askload5DImage(void);
	void load5DImage(std::vector<QStringList> filesChannTimeList, int numChann);

	void askload5DLabelImage(void);
	void load5DLabelImage(QStringList labelfilesChannTimeList);
	void update5DTable(void);

	void updateMultiLabels(void);
	


	void askLoadResult(void);
	void loadResult(QString fileName);
	void askLoadTable(void);
	void loadTable(QString fileName);
	void loadAdjTables(QString fileName);
	void loadProject(void);
	void saveDisplayImageToFile();
	void saveCompositeImageToFile();
	//Processing:
	void processProject(void);
	void startProcess(void);
	void process(void);
	void abortProcess(void);
	void continueProcess(void);
	void deleteProcess(void);

	//Saving:
	bool saveProject(void);
	bool saveSomaImage(void);
	bool askSaveImage(void);
	bool saveImage(void);
	bool askSaveResult(void);
	bool saveResult(void);
	bool askSaveTable(void);
	bool saveTable(void);
	bool askSaveAdjTables(void);
	bool saveAdjTables(void);
	void createDefaultLogName(void);
	
	//Views:
	void setPreferences();
	void toggleBounds();
	void toggleIDs();
	void toggleCentroids();
	void toggleCrosshairs();
	void toggleNucAdjacency();
	void toggleCellAdjacency();
	void toggleChannel(int chNum );
	void DisplayChannelsMenu();
#ifdef USE_OPENSLIDE
	void DisplayLevelsMenu();
	void toggleLevel(int levNum);
	void toggleSaveSlide(int levNum);
#endif
	void CreateNewPlotWindow();
	void CreateNewTableWindow();
	void CreateNewHistoWindow();
	void CreateNewNucRAG();
	void CreateNewCellRAG();
	void updateViews();
	void viewClosing(QWidget * view);
	void closeViews();

	//For Editing Menu
	void setCommonEnabled(bool val);
	void setEditsEnabled(bool val);
	void clearSelections(void);
	void addCell(int x1, int y1, int x2, int y2, int z);//Once box is drawn we call this to add a cell
	void mergeCells(void);
	void deleteCells(void);
	void splitCellAlongZ(void);		//Split single cell along current z
	void splitCells(void);
	void splitCell(int x1, int y1, int z1, int x2, int y2, int z2);
	void applyExclusionMargin(void);
	void changeClass(void);
	void markVisited(void);
	void fillCells(void);

	//***************************************************
	// Preprocessing Menu
	void CropToRegion(void);
	void InvertIntensities(void);
	void AnisotropicDiffusion(void);
	void MedianFilter(void);
	void SigmoidFilter(void);
	void GrayscaleErode(void);
	void GrayscaleDilate(void);
	void GrayscaleOpen(void);
	void GrayscaleClose(void);
	void CurvAnisotropicDiffusion(void);
	//void Resample(void);
	void preprocess(QString id);

	//***************************************************
	// Models Menu
	void createTrainer(void);
	void appendTrainer(void);

	//Queries Menu
	void queryKNearest(void);
	void queryInRadius(void);
	void queryViewsOff(void);
	double average(std::vector< std::pair<unsigned int, double> > ID);

	//*****************************************************

	//For Tools menu
	void getCentroids(void);
	void startROI(void);
	void endROI(void);
	void updateROIinTable(void);
	void loadROI(void);
	void saveROI(void);
	void clearROI(void);
	void roiStatistics(void);
	void preprocessImage(void);
	void unmixChannels(void);
	void segmentNuclei(void);
	void startEditing(void);
	void stopEditing(void);
	void startSVM();
	void updateDatabase();
	void startTraining();
	void startKPLS();
	void startActiveLearningwithFeat();
	void startActiveValidation();
	//void BuildGallery();
	void SaveActiveLearningResults(void);
	void SaveActiveLearningModel();
	void StartTraining(vtkSmartPointer<vtkTable> pTable);
	void ALDialogPopUP(bool first_pop, std::vector<std::pair<int,int> > query_labels);
	//void CreateActiveLearningModel(MCLR* mclr_alm,  vtkSmartPointer<vtkTable> pWizard_table);
	void classifyFromActiveLearningModel();
	void Start_Classification(bool create_model = false);
	std::vector< vtkSmartPointer<vtkTable> > Perform_Classification(std::vector< vtkSmartPointer<vtkTable> > table_vector);
	// Clus
	void runClus();
	
	//******************************************************
	//5D Views Menu
#ifdef USE_TRACKING
	void startTracking(void);
	void displayKymoGraph(void);
#endif
	void about(void);
	void menusEnabled(bool val);

signals:

protected:
	void createMenus();
	void createStatusBar();
	void createProcessToolBar();
	void createPreprocessingMenu();

	void writeSettings();
	void readSettings();

	void updateNucSeg(bool ask = false);

	bool askSaveChanges(QString text);
	int requestChannel(ftk::Image::Pointer img);	//Request a channel from this image
	QVector<QString> getChannelStrings(void);

	LabelImageViewQT *segView;
	std::vector<PlotWindow *> pltWin;
	std::vector<TableWindow *> tblWin;
	std::vector<HistoWindow *> hisWin;
	PatternAnalysisWizard *pWizard;

	QMenu *fileMenu;
	QAction *loadImageAction;
	QAction *load5DImageAction;
	QAction *load5DLabelImageAction;

	QAction *loadLabelAction;
	QAction *loadTableAction;
	QAction *loadProjectAction;
	QAction *processProjectAction;
	QAction *saveProjectAction;
	QAction *saveSomaImageAction;
	QAction *saveImageAction;
	QAction *saveLabelAction;
	QAction *saveTableAction;
	QAction *saveDisplayAction;
	QAction *saveCompositeAction;
	QAction *exitAction;
	vtkSmartPointer<vtkTable> model_table;

	QMenu *viewMenu;
	QAction *setPreferencesAction;
	QAction *showBoundsAction;
	QAction *showIDsAction;
	QAction *showCentroidsAction;
	QAction *showCrosshairsAction;
	QMenu *adjacencyMenu;
	QAction *showNucAdjAction;
	QAction *showCellAdjAction;
	QMenu *zoomMenu;
	QAction *zoomInAction;
	QAction *zoomOutAction;
	QMenu *displayChannelMenu;
#ifdef USE_OPENSLIDE
	QMenu *displayLevelMenu;
	QVector<QAction *> displayLevelAction;
	QSignalMapper *levSignalMapper;
	QSignalMapper *levSignalMapper1;
#endif
	QSignalMapper *chSignalMapper;
	QVector<QAction *> displayChannelAction;
	QAction *newTableAction;
	QAction *newScatterAction;
	QAction *newHistoAction;
	QMenu *ragMenu;
	QAction *nucRagAction;
	QAction *cellRagAction;
	QAction *imageIntensityAction;

	QMenu *toolMenu;
	QAction *getCentroidAction;
	QMenu *roiMenu;
	QAction *drawROIAction;
	QAction *loadROIAction;
	QAction *saveROIAction;
	QAction *clearROIAction;
	QAction *roiStatsAction;
	QAction *preprocessAction;
	QAction *unmixChannelsAction;
	QAction *segmentNucleiAction;
	QAction *editNucleiAction;
	QAction *svmAction;		//Start the One-Class SVM outlier detecter
	QAction *databaseAction;
	QMenu *activeMenu;
	QAction *activeAction; // Active Learning 
	QAction *activeVAction; // Active Learning 
	QAction *showGalleryAction;
	QAction *saveActiveResultsAction;
	QAction *saveActiveLearningModelAction;
	QAction *classifyFromActiveLearningModelAction;
	QMenu *classifyMenu;
	QAction *trainAction;	//Train the KPLS Classifier
	QAction *kplsAction;	//Start the KPLS Classifier
	QAction *runClusAction;

	//For Editing Menu
	QMenu *editMenu;
	QAction *clearSelectAction;
	QAction *addAction;
	QAction *mergeAction;
	QAction *deleteAction;
	QAction *fillAction;
	QAction *splitZAction;			//for split along z direction
	QAction *splitAction;			//for split along x-y direction
	QAction *exclusionAction;
	QAction *classAction;
	QAction *visitAction;			//Mark an object as visited

	QMenu *helpMenu;
	QAction *aboutAction;

    //For Models Menu
	QMenu *modelsMenu;
	QAction *createTrainingAction;
	QAction *appendTrainingAction;

	//For Queries Menu
	QMenu *queriesMenu;
	QAction *kNearestNeighborsAction;
	QAction *inRadiusNeighborsAction;
	QAction *queryViewsOffAction;
#ifdef USE_TRACKING
	//For 5D Image Menu
	QMenu * fiveDMenu;
	QAction * trackingAction;
	QAction * kymoViewAction;
#endif
	//************************************************************************
	//Preprocess menu
	QMenu *PreprocessMenu;
	QAction *cropAction;
	QAction *blankAction;
	QAction *invertPixAction;
	QAction *AnisotropicAction;
	QAction *MedianAction;
	QAction *SigmoidAction;
	QAction *GSErodeAction;
	QAction *GSDilateAction;
	QAction *GSOpenAction;
	QAction *GSCloseAction;
	QAction *CurvAnisotropicAction;
	//QAction *ResampleAction;

	//*********************************************************************
	QLabel *statusLabel;						//Shown at bottom of main window
	QString standardImageTypes;					//Image types that can be opened

	//Settings:
	QString lastPath;							//Last path that has been navigated to
	QMap<QString, QColor> colorItemsMap;		//Colors for specific objects

	ftk::NuclearSegmentation *nucSeg;			//Used for editing a nuclear segmentation
	ftk::ProjectProcessor *pProc;				//My project processor
	ftk::Image::Pointer myImg;					//My currently visible image
	ftk::Image::Pointer labImg;					//Currently visible label image
	//Amin
	ftk::Image::PtrMode mode;
	std::vector<QStringList> * filesChannTimeList;
	ftk::SpectralUnmixing * SpecUnmix;

#ifdef USE_TRACKING
	Image3DView * myview;
	TrackingKymoView * kymoView;
	MultiFrameCellTracker * mfcellTracker;
#endif

	ObjectSelection * selection;				//object selection list
	vtkSmartPointer<vtkTable> table;			//table
	vtkSmartPointer<vtkTable> pawTable;			//Pattern Analaysis Wizard (paw) table. Contains user selected features only
	vtkSmartPointer<vtkTable> NucAdjTable;		//Nuclear Adjacency Table
	vtkSmartPointer<vtkTable> CellAdjTable;		//Cellular Adjacency Table
	ftk::ProjectFiles projectFiles;				//files in the currently visible project
	ftk::ProjectDefinition projectDefinition;	//the project definition currently being used.
	unsigned int kplsRun;
	
	///////////////////////////////////////////////////////////////////////////////
	// Active Learning Variables
	///////////////////////////////////////////////////////////////////////////////
	unsigned int activeRun;
	MCLR *mclr;
	ActiveLearningDialog *dialog;
	bool from_model;
	std::vector< std::pair<int,int> > id_time;	
	std::vector<int> active_queries;
	std::vector<QImage> snapshots;
	vtkSmartPointer<vtkTable> featureTable;
	vnl_matrix<double> act_learn_matrix;
	double confidence_thresh;
	std::string classification_name;
	double sample_number;
	vnl_vector<double> std_dev_vec;
	vnl_vector<double> mean_vec; 
	std::vector< std::pair< std::string, vnl_vector<double> > >active_model;
	std::string prediction_col_name;
	std::string confidence_col_name;
	std::vector< std::vector< std::pair<int,int> > >  validation_samples; 
	// Gallery contains both the image of the query nuclei and their class values
	std::vector<std::pair<QImage,std::vector<int> > > gallery; 
	////////////////////////////////////////////////////////////////////////////////

	//This does not belong here, but is a temporary fix:
	void CreateDefaultAssociationRules();

	//Processing toolbar and thread pointers:
	bool abortProcessFlag;
	bool continueProcessFlag;
	QAction * processAbort;
	QAction * processContinue;
	QLabel * processTaskLabel;
	QProgressBar * processProgress;
	QToolBar * processToolbar;
	ProcessThread *processThread;

	//KPLS training and prediction
	std::vector<std::string> training_names;
	std::vector<std::string> prediction_names;
	std::string p_name;
	std::vector<std::string> class_names;
	int trainName;
	int predictName;

	//Clustering
	#ifdef	USE_Clusclus
		Heatmap *HeatmapWin;
	#endif

	// Loads all the tables of the time series
	std::vector< vtkSmartPointer<vtkTable> > tableVector;
	std::vector<std::string> imageNames;

	std::vector<std::pair<int,int> > ground_truth;


};

class ParamsFileDialog : public QDialog
{
	Q_OBJECT
public:
	ParamsFileDialog(QString lastPth, QVector<QString> chs, QWidget *parent = 0);
	QString getFileName();
	int getChannelNumber();
private slots:
	void ParamBrowse(QString);
private:
	QLabel *channelLabel;
	QComboBox *channelCombo;
	QRadioButton *autoButton;
	QRadioButton *fileButton;
	QComboBox *fileCombo;
	QString lastPath;
	QPushButton *okButton;
};

class ClassName_Confidence_Dialog : public QDialog
{
	Q_OBJECT
public:
	ClassName_Confidence_Dialog(QWidget *parent = 0);
	std::string getClassName();
	double getConfThresh();

private:
	QLabel *classNameLabel;
	QLineEdit *class_name;
	QHBoxLayout *classNameLayout;
	QLabel *confidenceLabel;
	QLineEdit *conf_thresh;
	QHBoxLayout *confLayout;
	QPushButton *okButton;
	QHBoxLayout *bLayout;
	QVBoxLayout *layout;
};



class SamplePercentDialog : public QDialog
{
	Q_OBJECT
public:
	SamplePercentDialog(int no_of_samples,int no_of_classes,QWidget *parent = 0);
	int samples;
	QSpinBox * vSpinNumber;
	int class_number;
	public slots:	
	
	void setNumber(double x);
	void setPercent(int x);

private:
	QLabel *sampleNumberLabel;
	QLabel *numberLabel;
	QLabel *sampleLabel;
	QLabel *dummyLabel;
	QLineEdit *sample_percent;
	QPushButton *okButton;
	QHBoxLayout *bLayout;
	QVBoxLayout *layout;
	QDoubleSpinBox * vSpinPercent;
	

};



class QueryDialog : public QDialog
{
	Q_OBJECT
public:
	QueryDialog(int QueryType, QVector<QString> classes, QWidget *parent = 0);
	std::vector<unsigned int> parseIDs(void);
	unsigned int parseK(void);
	double parseRad(void);
	unsigned short getSourceClass(void);
	unsigned short getDestClass(void);
	bool getKMutual(void);

private:
	QLabel * idLabel;
	QLineEdit * IDs;
	QHBoxLayout *idLayout;
	QLabel * kLabel;
	QLineEdit * K;
	QHBoxLayout *kLayout;
	QLabel * radLabel;
	QLineEdit * RAD;
	QHBoxLayout *radLayout;
	QLabel * classLabel1;
	QComboBox * classCombo1;
	QHBoxLayout * classLayout1;
	QLabel * classLabel2;
	QComboBox * classCombo2;
	QHBoxLayout * classLayout2;
	QCheckBox * check;
	QPushButton * okButton;
	QHBoxLayout * bLayout;
	QVBoxLayout * layout;
};

class ProcessThread : public QThread
{
public:
	ProcessThread(ftk::ProjectProcessor *proc = NULL);
	void run();
private:
	ftk::ProjectProcessor *myProc;
};

class PredictionDialog : public QDialog
{
	Q_OBJECT
public:
	PredictionDialog(QVector<QString> training_fields, QWidget *parent = 0);
	
public:
	int getTrainNumber();

private:
	QLabel *channelLabel;
	QComboBox *channelCombo;
	
	QPushButton *okButton;
	QPushButton *cancelButton;
};

#endif