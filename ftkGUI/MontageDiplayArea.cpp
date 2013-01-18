#include "MontageDiplayArea.h"

MontageDiplayArea::MontageDiplayArea(QWidget *parent)
  : QWidget(parent)
{
  channelImg = NULL;
  this->setupUI();

//  regionSelection = NULL; //Pointer to the region selection class ******

  setMouseTracking(true);
  rubberBand = NULL;

  channelFlags.clear();

}


MontageDiplayArea::~MontageDiplayArea()
{

}

void MontageDiplayArea::setupUI(void)
{

  //Setup Sliders
  vSlider = new QSlider();
  vSlider->setOrientation(Qt::Vertical);
  vSlider->setDisabled(true);
  vSlider->setRange(0,0);
  vSlider->setValue(0);

  QGridLayout *vsliderLayout = new QGridLayout;
  vsliderLayout->addWidget(vSlider,1,0,1,1);

  hSlider = new QSlider();
  hSlider->setOrientation(Qt::Horizontal);
  hSlider->setDisabled(true);
  hSlider->setRange(0,0);
  hSlider->setValue(0);

  QHBoxLayout *hsliderLayout = new QHBoxLayout;
  hsliderLayout->addWidget(hSlider);

  //Setup image display area
  imageLabel = new QLabel();
  imageLabel->setBackgroundRole(QPalette::Base);
  imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
  imageLabel->setMouseTracking(true);
  imageLabel->setScaledContents(true);
  imageLabel->resize(0,0);
  scrollArea = new QScrollArea();
  scrollArea->setMouseTracking(true);
  scrollArea->setBackgroundRole(QPalette::Dark);
  scrollArea->setWidget(imageLabel);
  scrollArea->horizontalScrollBar()->setRange(0,0);
  scrollArea->verticalScrollBar()->setRange(0,0);
  scrollArea->horizontalScrollBar()->setValue(0);
  scrollArea->verticalScrollBar()->setValue(0);

  QGridLayout *viewerLayout = new QGridLayout();
  viewerLayout->addWidget(scrollArea, 0, 0);
  viewerLayout->addLayout(vsliderLayout, 0, 1);
  viewerLayout->addLayout(hsliderLayout, 1, 0);

  QGridLayout *allLayout = new QGridLayout;
  allLayout->addLayout(viewerLayout,1,0);

  setLayout(allLayout);

  setAttribute ( Qt::WA_DeleteOnClose );
  setWindowTitle(tr("Downsampled Montage Viewer"));
}

void MontageDiplayArea::CreateRubberBand(void)
{
  if(!rubberBand)
    rubberBand = new MyRubberBand(this);
}

void MontageDiplayArea::mousePressEvent(QMouseEvent *event)
{
  origin = event->pos();		//Do i need this
  if(rubberBand)
  {
    rubberBand->setGeometry(QRect(origin,QSize()));
    rubberBand->show();
    return;
  }
}

void MontageDiplayArea::mouseMoveEvent(QMouseEvent *event)
{
  QPoint pos = event->pos();
  if(rubberBand)
  {
    if(pos.x() <= origin.x() || pos.y() <= origin.y())
      rubberBand->setGeometry(QRect(origin,QSize()));
    else
      rubberBand->setGeometry(QRect(origin, pos).normalized());
  }
}

void MontageDiplayArea::SetChannelImage(ftk::Image::Pointer img)
{
  if(!img)
  {
    channelImg = NULL;
    initChannelFlags();
    refreshBaseImage();
    return;
  }
  channelImg = img;
  initChannelFlags();
  refreshBaseImage();
}

void MontageDiplayArea::initChannelFlags()
{
  channelFlags.clear();
  std::vector<std::string> channel_names = channelImg->GetChannelNames();
  for (unsigned ch=0; ch<channel_names.size(); ++ch)
    channelFlags.push_back(true);
}

void MontageDiplayArea::refreshBaseImage()
{
  if(!channelImg)
    return;
  const ftk::Image::Info *info;
  if(channelImg)    info = channelImg->GetImageInfo();
  else return;
  itk::SizeValueType totalWidth  = (*info).numColumns;
  itk::SizeValueType totalHeight = (*info).numRows;
  baseImage = QImage(totalWidth, totalHeight, QImage::Format_ARGB32_Premultiplied);
  baseImage.fill(qRgb(0,0,0));
  
  QPainter painter(&baseImage);
  painter.setCompositionMode(QPainter::CompositionMode_Plus);

  for (int i=0; i < (*info).numChannels; i++)
  {
    if (channelFlags[i])
    {
      QImage gray((*info).numColumns, (*info).numRows, QImage::Format_ARGB32_Premultiplied);
      std::vector<unsigned char> color = (*info).channelColors[i];
      gray.fill(qRgb(color[0],color[1],color[2]));
      unsigned char * p = channelImg->GetSlicePtr<unsigned char>( 0, i, 0 );
      if(p)
      {
	//Get the image:
	QImage img(p, (*info).numColumns, (*info).numRows, (*info).numColumns, QImage::Format_Indexed8);
	gray.setAlphaChannel(img);      //Set it to the alpha channel
      }
      painter.drawImage(0,0,gray);
    }
  }
  this->repaint();
}

void MontageDiplayArea::moveEvent ( QMoveEvent * event )
{
	QWidget::moveEvent ( event );
}


void MontageDiplayArea::mouseReleaseEvent(QMouseEvent *event)
{
  if(rubberBand)
  {
    QPoint corner = scrollArea->pos();
    QPoint pos = event->pos() - corner;
    QPoint org = origin - corner;
    int x1 = (org.x() + scrollArea->horizontalScrollBar()->value());
    int y1 = (org.y() + scrollArea->verticalScrollBar()->value());
    int x2 = (pos.x() + scrollArea->horizontalScrollBar()->value());
    int y2 = (pos.y() + scrollArea->verticalScrollBar()->value());

    delete rubberBand;

    CreateRubberBand();
  }
}

void MontageDiplayArea::paintEvent(QPaintEvent * event)
{	
  QWidget::paintEvent(event);

  if(baseImage.height() <= 0 || baseImage.width() <= 0)
		return;

  displayImage = baseImage;
  QPainter painter(&displayImage);

  //Do zooming:
  int oldX = scrollArea->horizontalScrollBar()->value();
  int oldY = scrollArea->verticalScrollBar()->value();

  imageLabel->setPixmap(QPixmap::fromImage(displayImage));
  imageLabel->adjustSize();

  scrollArea->horizontalScrollBar()->setValue(oldX);
  scrollArea->verticalScrollBar()->setValue(oldY);
}

void MontageDiplayArea::SetChannelFlags( std::vector<bool> chFlags )
{
  channelFlags = chFlags;
  refreshBaseImage();
}
