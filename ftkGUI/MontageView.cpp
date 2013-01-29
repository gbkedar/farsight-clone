#include "MontageView.h"

MontageView::MontageView( QWidget * parent )
{
  setWindowTitle(tr("Montage View"));
  standardImageTypes = tr("XML Image Definition (*.xml)\n"
  			  "All Files (*.*)");
  Image = NULL;
  SubsampledImage = NULL;
  LabelImage = NULL;
  Table = NULL;
  tree = NULL;
  scaleFactor = DBL_MAX;

  chSignalMapper = NULL;

  imageViewer = new MontageDiplayArea();
  QPushButton *allButton = new QPushButton(tr("All"));
  this->setCentralWidget(imageViewer);

  this->createMenus();
  this->createToolBar();
  this->readSettings();
  this->connectSlotsToViewer();

  this->resize(300,100);
}

MontageView::~MontageView()
{

}

void MontageView::closeEvent(QCloseEvent *event)
{
	//Delete pointer from Nucleus Editor *****
}

void MontageView::readSettings()
{
  QSettings settings;
  lastPath = settings.value("lastPath", ".").toString();
}

void MontageView::writeSettings()
{
  QSettings settings;
  settings.setValue("lastPath", lastPath);
}

void MontageView::createMenus()
{
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->setObjectName("fileMenu");

  loadProjectAction = new QAction(tr("Load Project"), this);
  loadProjectAction->setObjectName("loadProjectAction");
  loadProjectAction->setStatusTip(tr("Load ftk project"));
  loadProjectAction->setShortcut(tr("Ctrl+L"));
  connect(loadProjectAction, SIGNAL(triggered()), this, SLOT(loadProject()));
  fileMenu->addAction(loadProjectAction);

  loadImageAction = new QAction(tr("Load Xml Image"), this);
  loadImageAction->setObjectName("loadImageAction");
  loadImageAction->setStatusTip(tr("Load an xml image"));
  loadImageAction->setShortcut(tr("Ctrl+O"));
  connect(loadImageAction, SIGNAL(triggered()), this, SLOT(askLoadImage()));
  fileMenu->addAction(loadImageAction);

  fileMenu->addSeparator();

  exitAction = new QAction(tr("Exit"), this);
  exitAction->setObjectName("exitAction");
  exitAction->setShortcut(tr("Ctrl+Q"));
  exitAction->setStatusTip(tr("Close the widget"));
  connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));
  fileMenu->addAction(exitAction);

  //VIEW MENU
  viewMenu = menuBar()->addMenu(tr("&View"));
  viewMenu->setObjectName("viewMenu");
  displayChannelMenu = viewMenu->addMenu(tr("Display Channel"));
  displayChannelMenu->setObjectName("displayChannelMenu");
  connect(displayChannelMenu, SIGNAL(aboutToShow()), this, SLOT(DisplayChannelsMenu()));

}

void MontageView::createToolBar()
{
  toolbar = this->addToolBar(tr("toolbar"));
  cropButton = new QPushButton(tr("Export Region"));
  connect(cropButton, SIGNAL(clicked()), this, SLOT(cropRegion()));
  cropButton->setDefault(false);
  cropButton->setAutoDefault(false);
  cropButton->setEnabled(false);
  toolbar->addWidget(cropButton);
}

void MontageView::connectSlotsToViewer()
{
  connect( imageViewer, SIGNAL(selectionDrawn(bool)), this, SLOT(enableRegionSelButton(bool)) );
}

void MontageView::enableRegionSelButton( bool enable )
{
  cropButton->setEnabled(enable);
}

void MontageView::DisplayChannelsMenu()
{
  if( !SubsampledImage )
    return;

  std::vector<std::string> channel_names = SubsampledImage->GetChannelNames();
  std::vector<bool> channel_status = imageViewer->GetChannelFlags();

  //remove all existing display channel actions
  for(unsigned i=0; i<displayChannelAction.size(); ++i)
    delete displayChannelAction.at(i);
  displayChannelMenu->clear();
  displayChannelAction.clear();

  if(chSignalMapper)
    delete chSignalMapper;

  chSignalMapper = new QSignalMapper(this);
  for(unsigned i=0; i<channel_names.size(); ++i)
  {
    QAction * action = new QAction( tr(channel_names.at(i).c_str()), this );
    action->setObjectName(tr(channel_names.at(i).c_str()));
    action->setCheckable(true);
    action->setChecked( channel_status.at(i) );
    action->setStatusTip( tr("Turn on/off this channel") );
    action->setShortcut( QString::number(i) );
    connect(action, SIGNAL(triggered()), chSignalMapper, SLOT(map()));
    chSignalMapper->setMapping( action, i );
    displayChannelMenu->addAction(action);
  }
  connect(chSignalMapper, SIGNAL(mapped(int)), this, SLOT(toggleChannel(int)));
}

void MontageView::askLoadImage()
{
  QString fileName = QFileDialog::getOpenFileName(this, "Open Image", lastPath, standardImageTypes);
  if(fileName == "")
    return; //No file returned
  if( this->loadImage(fileName) )
    resetSubsampledImageAndDisplayImage();
  else
    return;
}

bool MontageView::loadImage( QString fileName )
{
  lastPath = QFileInfo(fileName).absolutePath() + QDir::separator();
  QString name = QFileInfo(fileName).fileName();
  QString myExt = QFileInfo(fileName).suffix();
  if(myExt == "xml")
    Image = ftk::LoadXMLImage(fileName.toStdString());
  else
  {
    std::cerr<<"Can only load xml image files\n";
    return false;
  }
  //Clear image in nucleus editor *****
  //Clear label and table *****
  if(!Image)
  {
    std::cerr<<"Failed to load montage image\n";
    return false;
  }
  return true;
}

void MontageView::loadProject()
{
  QString filename = QFileDialog::getOpenFileName(this, "Open project...", lastPath, tr("XML Project File(*.xml)"));
  if(filename == "")
    return;

  QString path = QFileInfo(filename).absolutePath();
  lastPath = path;

  projectFiles.Read(filename.toStdString());

  if( this->loadImage( QString::fromStdString( projectFiles.GetFullInput()) ) )
    resetSubsampledImageAndDisplayImage();
  else return;

  projectDefinition.Load( projectFiles.GetFullDef() );

  QString labelName = QString::fromStdString( projectFiles.GetFullOutput() );
  QString myExt = QFileInfo( labelName ).suffix();
  if( myExt == "xml")
    LabelImage = ftk::LoadXMLImage( labelName.toStdString() );
  else
  {
    LabelImage = ftk::Image::New();
    if( !LabelImage->LoadFile( labelName.toStdString() ) )
      LabelImage = NULL;
  }

  if( !LabelImage )
  {
    std::cerr<<"Could not load the label image. Table will not be loaded\n";
    if( Table )
      Table = NULL;
    return;
  }

  Table = ftk::LoadTable( projectFiles.GetFullTable() );
  if(!Table)
  {
    std::cerr<<"Could not load the table\n";
    return;
  }

  this->IndexTable();

}

void MontageView::IndexTable()
{
  if(!Table)
  {
    std::cerr<<"Could not load the table\n";
    return;
  }
  //We only need the x-y as the whole z-stack is used when cropping
  SampleType::Pointer sample = SampleType::New();
  MeasurementVectorType mv;

  for( itk::SizeValueType i=0; i<Table->GetNumberOfRows(); ++i )
  {
    mv[0] = Table->GetValueByName(i,"centroid_x").ToDouble();
    mv[1] = Table->GetValueByName(i,"centroid_y").ToDouble();
    sample->PushBack( mv );
  }
  typedef itk::Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
  TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();

  treeGenerator->SetSample( sample );
  treeGenerator->SetBucketSize( 10 );
  treeGenerator->Update();
  tree = treeGenerator->GetOutput();
}

void MontageView::cropRegion()
{
  std::vector< int > coordinates = imageViewer->GetCoordinates();
  if( coordinates.size()!=4 )
  {
    std::cerr<<"Error in reading selection coordinates from montage display\n";
    return;
  }

  const ftk::Image::Info * imInfo = Image->GetImageInfo();

  //Translate the coordinates into the space of the original image
  itk::SizeValueType x1, y1, z1, x2, y2, z2;
  z1 = 0; z2 = imInfo->numZSlices-1;
  x1 = (((double)coordinates.at(0))*scaleFactor);
  y1 = (((double)coordinates.at(1))*scaleFactor);
  x2 = (((double)coordinates.at(2))*scaleFactor);
  y2 = (((double)coordinates.at(3))*scaleFactor);

  if( x2>=imInfo->numColumns )
    x2 = imInfo->numColumns-1;
  if( y2>imInfo->numRows )
    y2 = imInfo->numRows-1;

  if( x1>=x2 || y1>=y2 )
  {
    std::cerr<<"Error in reading selection coordinates from montage display\n";
    return;
  }

  //Crop channels
 // ftk::Image::Pointer cropChannelsImage = 

}

//Take the maximum intensity projection and subsample each channel
//Subsample the image setting the maximum x-y dimension to 1000pixels
//Redo the contrast, cast, smooth(for aliasing), then resample
void MontageView::resetSubsampledImageAndDisplayImage()
{
  typedef itk::Image< unsigned char, 3 > Uchar3DImageType;
  typedef itk::Image< float, 3 > Float3DImageType;
  typedef itk::MaximumProjectionImageFilter< Uchar3DImageType, Uchar3DImageType > MaxProjectType;
  typedef itk::AdaptiveHistogramEqualizationImageFilter< Uchar3DImageType > AdaptiveHistEqType;
  typedef itk::RescaleIntensityImageFilter< Uchar3DImageType >  RescaleIntensityType;
  typedef itk::CastImageFilter< Uchar3DImageType, Float3DImageType > CastFilterType;
  typedef itk::RecursiveGaussianImageFilter< Float3DImageType, Float3DImageType > GaussianFilterType;
  typedef itk::ResampleImageFilter< Float3DImageType, Uchar3DImageType > ResampleFilterType;
  typedef itk::IdentityTransform< double, 3 >  TransformType;
  typedef itk::LinearInterpolateImageFunction< Float3DImageType, double > InterpolatorType;

  if( SubsampledImage )
    delete SubsampledImage;

  SubsampledImage = ftk::Image::New();
  const ftk::Image::Info * imInfo = Image->GetImageInfo();
  std::vector< std::string > filenames = Image->GetFilenames();
  Uchar3DImageType::SizeType size;
  if( imInfo->numColumns<1000 && imInfo->numRows<1000 )
    SubsampledImage = Image; //Just in case...
  else
  {
   for( unsigned i=0; i<imInfo->channelNames.size(); ++i )
   {
    Uchar3DImageType::Pointer currentChannel = Image->GetItkPtr<unsigned char>( 0, i, ftk::Image::DEEP_COPY);
    Uchar3DImageType::Pointer currentChannelProjection;
    //Get projection if the image is 3D
    if( imInfo->numZSlices > 1 )
    {
      MaxProjectType::Pointer MaxProjectFilter = MaxProjectType::New();
      MaxProjectFilter->SetInput( currentChannel );
      MaxProjectFilter->SetProjectionDimension( 3 );
      try
      {
	MaxProjectFilter->Update();
      }
      catch( itk::ExceptionObject & excp )
      {
	std::cerr << "Exception in projection for subsampling:"
		  << excp << std::endl;
      }
      currentChannelProjection = MaxProjectFilter->GetOutput();
    }
    else
    {
      currentChannelProjection = currentChannel;
    }

    scaleFactor = imInfo->numColumns > imInfo->numRows ? imInfo->numColumns : imInfo->numRows ;
    scaleFactor = scaleFactor/1000.0;
    const Uchar3DImageType::SpacingType& inputSpacing = currentChannelProjection->GetSpacing();

    CastFilterType::Pointer caster = CastFilterType::New();
    caster->SetInput( currentChannelProjection );

    GaussianFilterType::Pointer smoother = GaussianFilterType::New();
    smoother->SetInput( caster->GetOutput() );
    smoother->SetSigma( inputSpacing[0]*(scaleFactor/2) );
    smoother->SetNormalizeAcrossScale( true );

    ResampleFilterType::Pointer resampler = ResampleFilterType::New();

    TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();
    resampler->SetTransform( transform );

    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    resampler->SetInterpolator( interpolator );

    resampler->SetDefaultPixelValue( 0 );

    Uchar3DImageType::SpacingType spacing;
    spacing[0] = inputSpacing[0] * scaleFactor;
    spacing[1] = inputSpacing[1] * scaleFactor;
    spacing[2] = inputSpacing[2];

    resampler->SetOutputSpacing( spacing );
    resampler->SetOutputOrigin( currentChannelProjection->GetOrigin() );
    resampler->SetOutputDirection( currentChannelProjection->GetDirection() );

    size[0] = ((double)imInfo->numColumns)/scaleFactor;
    size[1] = ((double)imInfo->numRows)/scaleFactor;
    size[2] = imInfo->numZSlices;
    resampler->SetSize( size );

    resampler->SetInput( smoother->GetOutput() );

/*  //Good idea to use but needs to be multithreaded.. Too slow..
    AdaptiveHistEqType::Pointer histeq = AdaptiveHistEqType::New();
    histeq->SetInput( resampler->GetOutput() );
    histeq->SetAlpha( 1 );
    histeq->SetBeta ( 1 );
    histeq->SetRadius( 100 );
*/
    RescaleIntensityType::Pointer rescale = RescaleIntensityType::New();
    rescale->SetInput( resampler->GetOutput() );
    rescale->SetOutputMaximum( itk::NumericTraits<Uchar3DImageType::PixelType>::max() );
    rescale->SetOutputMinimum( itk::NumericTraits<Uchar3DImageType::PixelType>::min() );

    std::string opstring = ftk::GetFilePath( filenames.at(i) )+ "/" 
    			   + imInfo->channelNames.at(i) + "_subsample.tif";
    typedef itk::ImageFileWriter< Uchar3DImageType > ImageFileWriterType;
    ImageFileWriterType::Pointer writer = ImageFileWriterType::New();
    writer->SetFileName( opstring.c_str() );
    writer->SetInput( rescale->GetOutput() );
    std::cout<<"Writing downsampled image:"<<opstring<<std::endl;

    try
    {
      std::cout<<"Resampling Channel "<< imInfo->channelNames.at(i) << std::endl;
      //rescale->Update();
      writer->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << "Exception in subsampling:\n"
		<< excep << std::endl;
    }
    SubsampledImage->AppendChannelFromData3D( (void *)resampler->GetOutput()->GetBufferPointer(),
	itk::ImageIOBase::UCHAR, sizeof(Uchar3DImageType::PixelType), size[0], size[1], size[2],
	imInfo->channelNames.at(i), imInfo->channelColors.at(i), true );
    std::cout	<<"Channel image name "<<SubsampledImage->GetChannelNames().at(i)
		<<"\tsize:\tz="<<size[2]<<"\ty="<<size[1]<<"\tx="<<size[0]<<std::endl;
   }
  }
  int resizeRows, resizeCols;
  resizeRows = (int)(size[1]+78);  //Always less than 1000
  resizeCols = (int)(size[0]+20);  //Always less than 1000
  this->resize(resizeCols,resizeRows); //Image size plus some for the menus, toolbar and bezel
  SetChannelImage();
}

void MontageView::SetChannelImage()
{
  imageViewer->SetChannelImage( SubsampledImage );
  //Clear nucleus editor views, label and table ***************
}

void MontageView::toggleChannel( int chNum )
{
	std::vector<bool> ch_stats = imageViewer->GetChannelFlags();
	ch_stats[chNum] = !ch_stats[chNum];
	imageViewer->SetChannelFlags( ch_stats );
}
