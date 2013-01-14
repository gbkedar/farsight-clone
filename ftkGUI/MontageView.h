#ifndef MONTAGE_VIEW_H
#define MONTAGE_VIEW_H

#include <iostream>

#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QAction>
#include <QtGui/QWidget>
#include <QApplication>
#include <QFileDialog>

class MontageView : public QMainWindow
{
  Q_OBJECT;
public:
  MontageView(QWidget * parent = 0);
  ~MontageView();
  void Initialize();

private:

protected:
  virtual void closeEvent(QCloseEvent *event);
protected slots:
  void loadProject(void);
  void askLoadImage(void);
  void DisplayChannelsMenu(void);
  void toggleChannel(int chNum);
  void loadImage(QString fileName);


private:
  QVTKWidget MainQTWidget;
  ftk::Image::Pointer Image;
  ftk::Image::Pointer SubsampledImage;
  ftk::Image::Pointer LabelImage;
  vtkSmartPointer<vtkTable> Table;

signals:
protected:
  void createMenus();
  void readSettings();
  QMenu *fileMenu;
  QAction *loadProjectAction;
  QAction *loadImageAction;
  QAction *exitAction;

  QMenu *viewMenu;
  QMenu *displayChannelMenu;
  QSignalMapper *chSignalMapper;
  QVector<QAction *> displayChannelAction;
  QString standardImageTypes;

  QString lastPath;
}
#endif //MONTAGE_VIEW_H
