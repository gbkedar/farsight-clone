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
#include "itkResampleImageFilter.h"

#include "MontageDiplayArea.h"

//Testing
#include "ftkCommon/ftkUtils.h"
#include "itkImageFileWriter.h"

class MontageView : public QMainWindow
{
  Q_OBJECT;
public:
  MontageView(QWidget * parent = 0);
  ~MontageView();
  void Initialize();

//private:

protected:
  virtual void closeEvent(QCloseEvent *event);
  MontageDiplayArea *imageViewer;

protected slots:
  void loadProject(void);
  void askLoadImage(void);
  void DisplayChannelsMenu(void);
  void cropRegion(void);
  void toggleChannel(int chNum);
  void loadImage(QString fileName);
  void resetSubsampledImageAndDisplayImage(void);
  void SetChannelImage(void);
  void enableRegionSelButton(bool);

private:
  ftk::Image::Pointer Image;
  ftk::Image::Pointer SubsampledImage;
  ftk::Image::Pointer LabelImage;
  vtkSmartPointer<vtkTable> Table;

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
#endif //MONTAGE_VIEW_H
