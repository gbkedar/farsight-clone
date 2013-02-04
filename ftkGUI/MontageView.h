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
#include "itkRelabelComponentImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkKdTreeGenerator.h"

#include "MontageDiplayArea.h"

//Testing
#include "ftkCommon/ftkUtils.h"
#include "itkImageFileWriter.h"
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

public:
  struct TableEntryList{
    itk::SizeValueType LabelImId, TabInd;
    LabelStatisticsImageFilterType::BoundingBoxType BoundBox;
  };
  struct TableEntryComparator
  {  bool operator()( const TableEntryList& t, itk::SizeValueType Value ) const
     { return t.BoundBox.at(0) < Value; }
     bool operator()( itk::SizeValueType Value, const TableEntryList& t ) const
     { return Value < t.BoundBox.at(0); }
     //Following operator:Indexing table and image initiall-->Can afford to be inefficient
     bool operator()( const TableEntryList& t1, const TableEntryList& t2 ) const
     { if( t1.BoundBox.at(0) == t2.BoundBox.at(0) ) return t1.BoundBox.at(1) < t2.BoundBox.at(2);
       return t1.BoundBox.at(0) < t2.BoundBox.at(0); }
  };

//private:

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

  //Utility functions
  void cropRegion(void);
  vtkSmartPointer<vtkTable> GetCroppedTable( itk::SizeValueType x1, itk::SizeValueType y1,
	itk::SizeValueType x2, itk::SizeValueType y2, NucleusEditorLabelType::Pointer CropLabel );

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
