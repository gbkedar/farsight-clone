#ifndef MONTAGE_VIEW_H
#define MONTAGE_VIEW_H

#include <iostream>

#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QAction>
#include <QtGui/QWidget>
#include <QtGui/QToolBar>
#include <QtGui/QButtonGroup>
#include <QtCore/QSignalMapper>
#include <QApplication>
#include <QFileDialog>

#include "itkMaximumProjectionImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkAdaptiveHistogramEqualizationImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkResampleImageFilter.h"
//#include "itkKdTreeGenerator.h"

#include "MontageDiplayArea.h"
#include "MontageRegionSelection.h"

#include "ftkCommon/ftkUtils.h"
#include "NuclearSegmentation/NucleusEditor/ftkProjectProcessor.h"

#include <boost/bind.hpp>

class MontageView : public QMainWindow
{
  Q_OBJECT;
public:
  MontageView(QWidget * parent = 0);
  ~MontageView();
  void Initialize();
  MontageRegionSelection *RegionSelection;

  struct TableEntryList{
    itk::SizeValueType x, y, z, LabelImId, xDownSampled, yDownSampled;
  };

signals:
  void NewRegionSelected();
  void MontageViewClosed();

protected:
  virtual void closeEvent(QCloseEvent *event);
  MontageDiplayArea *imageViewer;
  CellTypingDialog *cellTypingDialog;

protected slots:
  void loadProject(void);
  void askLoadImage(void);
  void DisplayChannelsMenu(void);
  void toggleChannel(int chNum);
  bool loadImage(QString fileName);
  void resetSubsampledImageAndDisplayImage(void);
  void SetChannelImage(void);
  void enableRegionSelButton(bool);
  void cropRegion(void);
  void processProject(void);
  void LaunchCellTypingWindow(void);
  void cellTypingDialogClosing(void);

private:
  typedef itk::Image<unsigned short, 3> NucleusEditorLabelType;
  typedef itk::Image<unsigned int,   3> LabelUnsignedIntType;
  typedef itk::LabelStatisticsImageFilter
    < LabelUnsignedIntType, LabelUnsignedIntType > LabelStatisticsImageFilterType;
  /*typedef itk::Vector< float, 2 > MeasurementVectorType;
  typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
  typedef itk::Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
  typedef TreeGeneratorType::KdTreeType TreeType;
  TreeType::Pointer tree;		//Should switch Faster indexing of the table*/

  ftk::Image::Pointer Image;
  ftk::Image::Pointer SubsampledImage;
  ftk::Image::Pointer LabelImage;
  vtkSmartPointer<vtkTable> Table;
  ftk::ProjectFiles projectFiles;			//files in the currently visible project
  ftk::ProjectDefinition projectDefinition;		//the project definition currently being used.
  double scaleFactor;
  itk::SizeValueType NumberOfLabelsFound;
  int nucChannel;
  CellTypingDialog *cellTypeDialog;

  template<typename pixelType> NucleusEditorLabelType::Pointer
  				RelabelImage(ftk::Image::Pointer InputImage);
  std::vector< TableEntryList > TableEntryVector;
  std::vector< LabelStatisticsImageFilterType::BoundingBoxType > BoundingBoxes;
  std::map< itk::SizeValueType, itk::SizeValueType > LabelToTableMap;
  std::map< itk::SizeValueType, itk::SizeValueType > LabelToRelabelMap;
  std::map<int, ftk::Object::Point> CroppedCentroidsMap;
  std::map<int, ftk::Object::Box> CroppedBoundBoxesMap;
  std::vector< std::vector <std::string> > ClassificationGroups;

  //Utility functions
  vtkSmartPointer<vtkTable> GetCroppedTable( itk::SizeValueType x1, itk::SizeValueType y1,
	itk::SizeValueType x2, itk::SizeValueType y2 );
  itk::SizeValueType InsertNewLabelToRelabelMap( itk::SizeValueType NewKey );
  itk::SizeValueType CheckBoundsAndSubtractMin( itk::SizeValueType CoOrd,
  					itk::SizeValueType Min, itk::SizeValueType Max );
  void Process(void);
  void computeBoundBoxesAndTableMapsForDisplay(void);

protected:
  void createMenus();
  void readSettings();
  void writeSettings();
  void connectSlotsToViewer();
  void createToolBar();

  QMenu *fileMenu;
  QAction *loadProjectAction;
  QAction *loadImageAction;
  QAction *processAction;
  QAction *exitAction;

  QMenu *viewMenu;
  QMenu *displayChannelMenu;
  QSignalMapper *chSignalMapper;
  QVector<QAction *> displayChannelAction;
  QString standardImageTypes;

  QMenu *toolsMenu;
  QAction *cellTypeDialogMenu;

  QPushButton *cropButton;
  QToolBar *toolbar;
  QString lastPath;
};

#include "MontageView.txx"
#endif //MONTAGE_VIEW_H
