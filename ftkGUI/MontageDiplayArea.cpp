#include "MontageDiplayArea.h"

MontageDiplayArea::MontageDiplayArea(QWidget *parent)
  : QWidget(parent)
{
  channelImg = NULL;
  totalWidth = totalHeight = 0;
  channelFlags.clear();
  centroidPoints.clear();

  this->setupUI();

//  regionSelection = NULL; //Pointer to the region selection class ******

  setMouseTracking(true);
  rubberBand = NULL;
  mousePress = false;
  paintingCentroidClass = false;

  numPointsToPaint    = 0;
  paintCentroidsColor = Qt::green;		
}


MontageDiplayArea::~MontageDiplayArea()
{

}

void MontageDiplayArea::setupUI(void)
{
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

  QGridLayout *allLayout = new QGridLayout;
  allLayout->addLayout(viewerLayout,1,0);

  setLayout(allLayout);

  setAttribute ( Qt::WA_DeleteOnClose );
  setWindowTitle(tr("Downsampled Montage Viewer"));
}

void MontageDiplayArea::CreateRubberBand(void)
{
  if(rubberBand)
    delete rubberBand;
  rubberBand = new MyRubberBand(this);
}

void MontageDiplayArea::mousePressEvent(QMouseEvent *event)
{
  if(!channelImg)
    return;

  mousePress = true;
  this->CreateRubberBand();
  origin = event->pos();
  rubberBand->setGeometry(QRect(origin,QSize()));
  rubberBand->show();
  return;
}

void MontageDiplayArea::mouseMoveEvent(QMouseEvent *event)
{
  QPoint pos = event->pos();
  if(mousePress)
  {
    QPoint corner = scrollArea->pos();
    QPoint v_pos = pos - corner;  // This is a local position (in viewport coordinates)
    int xx = v_pos.x() + scrollArea->horizontalScrollBar()->value();
    int yy = v_pos.y() + scrollArea->verticalScrollBar()->value();
    if( xx>=0 && xx<totalWidth && yy>=0 && yy<totalHeight )
    if(pos.x() <= origin.x() || pos.y() <= origin.y())
    {
      rubberBand->setGeometry(QRect(origin,QSize()));
    }
    else
    {
      rubberBand->setGeometry(QRect(origin, pos).normalized());
    }
  }
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
    x1 = org.x() + scrollArea->horizontalScrollBar()->value();
    y1 = org.y() + scrollArea->verticalScrollBar()->value();
    x2 = pos.x() + scrollArea->horizontalScrollBar()->value();
    y2 = pos.y() + scrollArea->verticalScrollBar()->value();
    rubberBand->setMouseTracking(false);
    mousePress = false;
    if( x1<x2 && y1<y2 )
      emit selectionDrawn(true);
    else
      emit selectionDrawn(false);
  }
}

std::vector< int > MontageDiplayArea::GetCoordinates()
{
  std::vector< int > outputVector;
  if(rubberBand)
  {
    outputVector.push_back( x1 );
    outputVector.push_back( y1 );
    outputVector.push_back( x2 );
    outputVector.push_back( y2 );
  }
  else
  {
    outputVector.clear();
  }
    return outputVector;
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
  totalWidth  = (*info).numColumns;
  totalHeight = (*info).numRows;
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

void MontageDiplayArea::paintEvent(QPaintEvent * event)
{

  QWidget::paintEvent(event);

  if(baseImage.height() <= 0 || baseImage.width() <= 0)
		return;

  displayImage = baseImage; //May need to be replaced by deep copy
  if( paintingCentroidClass )
    this->repaintCentroids(0);
  QPainter painter(&displayImage);

  int oldX = scrollArea->horizontalScrollBar()->value();
  int oldY = scrollArea->verticalScrollBar()->value();

  imageLabel->setPixmap(QPixmap::fromImage(displayImage));
  imageLabel->adjustSize();

  scrollArea->horizontalScrollBar()->setValue(oldX);
  scrollArea->verticalScrollBar()->setValue(oldY);
}

void MontageDiplayArea::repaintCentroids( itk::SizeValueType numPointsToPaint = 0 )
{
  if( !numPointsToPaint )
    numPointsToPaint = centroidPoints.size();
    unsigned char *buf = displayImage.bits();
  const ftk::Image::Info *info = channelImg->GetImageInfo();
#ifdef _OPENMP
  #pragma omp parallel for
#if _OPENMP < 200805L
  for( itk::IndexValueType i=0; i<numPointsToPaint; ++i )
#else
  for(  itk::SizeValueType i=0; i<numPointsToPaint; ++i )
#endif
#else
  for(  itk::SizeValueType i=0; i<numPointsToPaint; ++i )
#endif
  {	
    QImage myImg = QImage(buf, (*info).numColumns, (*info).numRows, QImage::Format_ARGB32_Premultiplied);
    QPainter painter(&myImg);
    int x = centroidPoints.at( centroidPoints.size()-i-1 ).first;
    int y = centroidPoints.at( centroidPoints.size()-i-1 ).second;
    painter.setBrush(paintCentroidsColor);
    painter.drawEllipse(x-1, y-1, 2, 2);
  }
}

void MontageDiplayArea::SetChannelFlags( std::vector<bool> chFlags )
{
  channelFlags = chFlags;
  refreshBaseImage();
}

void MontageDiplayArea::respondToSlider(
	std::vector< std::pair<itk::SizeValueType,itk::SizeValueType> > &paintPoints,
	QColor &paintColor, unsigned &displayChannel )
{
  paintingCentroidClass = true;
  if( displayChannel || paintPoints.size()<centroidPoints.size() )
  //Display image needs to reset if channels are changed or some centroids need to be erased
  {
    numPointsToPaint = 0; //Need to repaint all in vector
    paintCentroidsColor = paintColor;
    centroidPoints = paintPoints;
    if( displayChannel )
    //Channel changing base image needs to be repainted
    {
      for( unsigned i=0; i<channelFlags.size(); ++i )
	if( (displayChannel-1)==i )
	  channelFlags.at(i) = true;
	else
	  channelFlags.at(i) = false;
      refreshBaseImage();
    }
    else
    //Just need to delete all centroids and repaint them
    {
      displayImage = baseImage; //May need to be replaced by deep copy
      this->repaintCentroids();
    }
  }
  else
  //Just adding points to base image
  {
    numPointsToPaint = paintPoints.size()-centroidPoints.size();
    centroidPoints = paintPoints;
    this->repaintCentroids( numPointsToPaint );
  }
}