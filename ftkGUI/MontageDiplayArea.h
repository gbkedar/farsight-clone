#ifndef MONTAGE_DISPLAY_H
#define MONTAGE_DISPLAY_H

#include <QtGui/QWidget>
#include <QtGui/QPainter>
#include <QtGui/QLabel>
#include <QtGui/QScrollBar>
#include <QtGui/QScrollArea>
#include <QtGui/QStylePainter>
#include <QtGui/QMouseEvent>
#include <QtCore/QRect>
#include <QtCore/QSize>
#include <QtCore/QPoint>
#include <QtCore/QMap>
#include <QtCore/QSettings>

#include <ftkImage/ftkImage.h>

#include <iostream>

#include "LabelImageViewQT.h"

class MontageDiplayArea  : public QWidget
{
  Q_OBJECT

public:
  MontageDiplayArea (QMap<QString, QColor> * new_colorItemsMap = NULL, QWidget *parent = 0);
  ~MontageDiplayArea();

  void SetChannelImage(ftk::Image::Pointer img);

public slots:


signals:
    void boxDrawn(int x1, int y1, int x2, int y2, int z);

protected slots:
  void refreshBaseImage(void);
  void updateVSlider(void);
  void updateHSlider(void);
  void initChannelFlags(void);

protected:
  void moveEvent (QMoveEvent * event);
  void mousePressEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *event);
  void paintEvent(QPaintEvent *event);
  void writeSettings();
  void readSettings();
  void setupUI(void);
  void CreateRubberBand(void);

  QMap<QString, QColor> * colorItemsMap;
  //UI Widgets:
  QScrollArea *scrollArea;
  QSlider *hSlider;
  QSlider *vSlider;

  ftk::Image::Pointer channelImg;

  std::vector<bool> channelFlags;	//is channel is visible or not

  QImage displayImage;			//Currently displayed image

  //For Getting a Box:
  QPoint origin;
  MyRubberBand *rubberBand;
}

#endif //MONTAGE_DISPLAY_H
