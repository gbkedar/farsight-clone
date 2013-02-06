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

//Testing
#include "itkImageFileWriter.h"

#include "ftkCommon/ftkUtils.h"
#include "NuclearSegmentation/NucleusEditor/ftkProjectFiles.h"
#include "NuclearSegmentation/NucleusEditor/ftkProjectDefinition.h"

#include <boost/bind.hpp>

class MontageView : public QMainWindow
{
  Q_OBJECT;
public:
  MontageView(QWidget * parent = 0);
  ~MontageView();
  void Initialize();
  MontageRegionSelection *RegionSelection;

private:

public:
  struct TableEntryList{
    itk::SizeValueType x, y, z, LabelImId, TabInd;
  };
  struct TableEntryComparator
  {  bool operator()( const TableEntryList& t, itk::SizeValueType Value ) const
     { return t.x < Value; }
     bool operator()( itk::SizeValueType Value, const TableEntryList& t ) const
     { return Value < t.x; }
     //Following operator:Indexing table and image initiall-->Can afford to be inefficient
     bool operator()( const TableEntryList& t1, const TableEntryList& t2 ) const
     { if( t1.x == t2.x ) return t1.y < t2.y;
       else return t1.x < t2.x; }
  };

protected:
  virtual void closeEvent(QCloseEvent *event);
  MontageDiplayArea *imageViewer;

protected slots:
  void loadProject(void);
  void askLoadImage(void);
  void DisplayChannelsMenu(void);
  void toggleChannel(int chNum);
  bool loadImage(QString fileName);
  void resetSubsampledImageAndDisplayImage(void);
  void SetChannelImage(void);
  void enableRegionSelButton(bool);

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

  template<typename pixelType> NucleusEditorLabelType::Pointer
  				RelabelImage(ftk::Image::Pointer InputImage);
  std::vector< TableEntryList > TableEntryVector;
  std::vector< LabelStatisticsImageFilterType::BoundingBoxType > BoundingBoxes;
  std::map< itk::SizeValueType, itk::SizeValueType > LabelToTableMap;
  std::map< itk::SizeValueType, itk::SizeValueType > LabelToRelabelMap;
  std::map<int, ftk::Object::Point> CroppedCentroidsMap;
  std::map<int, ftk::Object::Box> CroppedBoundBoxesMap;

  //Utility functions
  void cropRegion(void);
  vtkSmartPointer<vtkTable> GetCroppedTable( itk::SizeValueType x1, itk::SizeValueType y1,
	itk::SizeValueType x2, itk::SizeValueType y2 );
  itk::SizeValueType InsertNewLabelToRelabelMap( itk::SizeValueType NewKey );
  itk::SizeValueType CheckBoundsAndSubtractMin( itk::SizeValueType CoOrd,
  					itk::SizeValueType Min, itk::SizeValueType Max );

signals:
protected:
  void createMenus();
  void readSettings();
  void writeSettings();
  void connectSlotsToViewer();
  void createToolBar();

  QMenu *fileMenu;
  QAction *loadProjectAction;
  QAction *loadImageAction;
  QAction *exitAction;

  QMenu *viewMenu;
  QMenu *displayChannelMenu;
  QSignalMapper *chSignalMapper;
  QVector<QAction *> displayChannelAction;
  QString standardImageTypes;

  QPushButton *cropButton;
  QToolBar *toolbar;
  QString lastPath;
};

#include "MontageView.txx"
#endif //MONTAGE_VIEW_H
