#include "MontageView.h"

#define DEBUG_MONTV 1

MontageView::MontageView( QWidget * parent )
{
  setWindowTitle(tr("Montage View"));
  standardImageTypes = tr("XML Image Definition (*.xml)\n"
  			  "All Files (*.*)");
  Image = NULL;
  SubsampledImage = NULL;
  LabelImage = NULL;
  Table = NULL;
  //tree = NULL;
  TableEntryVector.clear();
  BoundingBoxes.clear();
  LabelToRelabelMap.clear();
  LabelToTableMap.clear();
  CroppedCentroidsMap.clear();
  CroppedBoundBoxesMap.clear();
  scaleFactor = DBL_MAX;
  NumberOfLabelsFound = 0;

  chSignalMapper = NULL;
  RegionSelection = NULL;

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
  if( RegionSelection )
  {
    delete RegionSelection;
    RegionSelection = NULL;
  }
}

void MontageView::closeEvent(QCloseEvent *event)
{
  //Check if there are any edits made in the Nucleus Editor through RegionSelection*****
  emit MontageViewClosed();
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

  TableEntryVector.clear();
  BoundingBoxes.clear();
  LabelToTableMap.clear();

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

  std::cout<<"Computing label stats for the label image\n";

  LabelStatisticsImageFilterType::Pointer LabelStats = LabelStatisticsImageFilterType::New();
  LabelStats->SetInput( LabelImage->GetItkPtr<unsigned int>(0,0) );
  LabelStats->SetLabelInput( LabelImage->GetItkPtr<unsigned int>(0,0) );
  LabelStats->UseHistogramsOff();

  try
  {
    LabelStats->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Error in Label Image stats" << excp << std::endl;
  }
  for( LabelStatisticsImageFilterType::ValidLabelValuesContainerType::const_iterator
	  it = LabelStats->GetValidLabelValues().begin()+1;/*Since label 0 counted*/
	it != LabelStats->GetValidLabelValues().end(); ++it )
  {
    LabelStatisticsImageFilterType::BoundingBoxType bb = LabelStats->GetBoundingBox( *it );
    BoundingBoxes.push_back( bb );
  }

  Table = ftk::LoadTable( projectFiles.GetFullTable() );
  if( !Table || ( Table->GetNumberOfRows() != BoundingBoxes.size() ) )
  {
    if( Table->GetNumberOfRows() != BoundingBoxes.size() )
      std::cerr<< "The label image and the table do not have the same number of IDs\n"
		<< Table->GetNumberOfRows() << " rows, " << BoundingBoxes.size()
		<< " labels\n";
    else
      std::cerr<< "No Table Loaded\n";
    std::cerr<< "Setting Label Image and Table to NULL\n";
    Table = NULL;
    LabelImage =NULL;
    if( !BoundingBoxes.empty() )
      BoundingBoxes.clear();
    return;
  }

  for( itk::SizeValueType i=0; i<Table->GetNumberOfRows(); ++i )
  {
    TableEntryList CurrentEntry;
    CurrentEntry.x  = Table->GetValueByName(i,"centroid_x").ToTypeUInt64();
    CurrentEntry.y  = Table->GetValueByName(i,"centroid_y").ToTypeUInt64();
    CurrentEntry.z  = Table->GetValueByName(i,"centroid_z").ToTypeUInt64();
    CurrentEntry.LabelImId = Table->GetValueByName(i,"ID").ToTypeUInt64();
    if(!(CurrentEntry.x>=BoundingBoxes.at(i).at(0) && CurrentEntry.x<=BoundingBoxes.at(i).at(1)
      && CurrentEntry.y>=BoundingBoxes.at(i).at(2) && CurrentEntry.y<=BoundingBoxes.at(i).at(3)
      && CurrentEntry.z>=BoundingBoxes.at(i).at(4) && CurrentEntry.z<=BoundingBoxes.at(i).at(5) ) )
      std::cerr<< i << "The centroid in the table and label image don't match for id "
		<< CurrentEntry.LabelImId << " Nucleus editor may crash if you click on this cell\n";
    if( LabelStats->GetValidLabelValues().at(i+1)/*Since label 0 counted*/!=CurrentEntry.LabelImId )
      std::cerr<< LabelStats->GetValidLabelValues().at(i+1) << " The ids don't match for label " 
		<< CurrentEntry.LabelImId << " and table entry number " << i
		<<" Nucleus editor may crash if you click on this cell\n";
    CurrentEntry.TabInd = i;
    TableEntryVector.push_back( CurrentEntry );
    LabelToTableMap.insert(std::map<itk::SizeValueType, itk::SizeValueType>::value_type
    				( CurrentEntry.LabelImId, i ) );
  }
  std::sort( TableEntryVector.begin(), TableEntryVector.end(), MontageView::TableEntryComparator() );
}

void MontageView::cropRegion()
{
  this->enableRegionSelButton(false);
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

  int bytePerPix = 1;
  //Crop channels
  ftk::Image::Pointer CropChannel = Image->
    CropImage<unsigned char>( x1, y1, z1, x2, y2, z2, itk::ImageIOBase::UCHAR, bytePerPix );

  NucleusEditorLabelType::Pointer CropLabelItk;
  vtkSmartPointer<vtkTable> CropTable;
  LabelToRelabelMap.clear();

  if( RegionSelection )
  {
    delete RegionSelection;
    RegionSelection = NULL;
  }
  RegionSelection = new MontageRegionSelection;
  RegionSelection->SetChannelImage( CropChannel );

  //Check if there are any labels and if there are get their table entries
  if( LabelImage )
  {
    if( LabelImage->IsMatch<unsigned short>( LabelImage->GetImageInfo()->dataType ) )
    {
      bytePerPix = 2;
      CropLabelItk = this->RelabelImage<unsigned short>(
			LabelImage->CropImage<unsigned short>
			( x1, y1, z1, x2, y2, z2, itk::ImageIOBase::USHORT, bytePerPix ) );
    }
    else if( LabelImage->IsMatch<unsigned int>( LabelImage->GetImageInfo()->dataType ) )
    {
      bytePerPix = 4;
      CropLabelItk = this->RelabelImage<unsigned int>(
			LabelImage->CropImage<unsigned int>
			( x1, y1, z1, x2, y2, z2, itk::ImageIOBase::UINT, bytePerPix ) );
    }
    else
    {
      std::cerr<<"Unrecognized label image. unsigned and unsigned short are the valid types\n";
      return;
    }
    std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer CropLabel = ftk::Image::New();
    CropLabel->AppendChannelFromData3D( CropLabelItk->GetBufferPointer(), itk::ImageIOBase::USHORT,
	sizeof(unsigned short), (x2-x1+1), (y2-y1+1), (z2-z1+1), "nuc", color, true );
    RegionSelection->SetLabelImage( CropLabel );
    if( NumberOfLabelsFound && Table )
    {
      CropTable = this->GetCroppedTable( x1, y1, x2, y2 );
      RegionSelection->SetTable( CropTable );
      RegionSelection->SetBoundsMap( CroppedBoundBoxesMap );
      RegionSelection->SetCenterMap( CroppedCentroidsMap );
    }
  }
#if DEBUG_MONTV
  if( LabelImage )
  {
    typedef itk::ImageFileWriter< NucleusEditorLabelType > BinaryWriterType;
    BinaryWriterType::Pointer binwriter = BinaryWriterType::New();
    binwriter->SetInput( CropLabelItk );
    std::string binWriterStr = "label_crop.tif";
    binwriter->SetFileName( binWriterStr.c_str() );
    try{ binwriter->Update(); }
    catch( itk::ExceptionObject & excp ){ std::cerr << "Crop Label writer" << excp << std::endl; }
  }

  if( CropTable )
  {
    ftk::SaveTable( "crop_table.txt", CropTable );
  }

  typedef itk::ImageFileWriter< itk::Image<unsigned char, 3> > BinaryWriterType1;
  BinaryWriterType1::Pointer binwriter1 = BinaryWriterType1::New();
  binwriter1->SetInput( CropChannel->GetItkPtr<unsigned char>( 0, 0) );
  std::string binWriterStr1 = "channel0_crop.tif";
  binwriter1->SetFileName( binWriterStr1.c_str() );
  try{ binwriter1->Update(); }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Crop Image writer" << 
      CropChannel->GetImageInfo()->numColumns << "\t" <<
      CropChannel->GetImageInfo()->numRows << excp << std::endl;
  }
#endif
  emit NewRegionSelected();
}

vtkSmartPointer<vtkTable> MontageView::GetCroppedTable( itk::SizeValueType x1, itk::SizeValueType y1,
		itk::SizeValueType x2, itk::SizeValueType y2 )
{
  vtkSmartPointer<vtkTable> cropTable = vtkSmartPointer<vtkTable>::New();
  for( vtkIdType i = 0; i<Table->GetNumberOfColumns(); ++i )
  {
    vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
    column->SetName( Table->GetColumnName(i) );
    cropTable->AddColumn(column);
  }

  for( std::map< itk::SizeValueType , itk::SizeValueType >::iterator it = LabelToRelabelMap.begin();
	it != LabelToRelabelMap.end(); ++it )
  {
    vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
    std::map< itk::SizeValueType , itk::SizeValueType >::iterator it1;
    it1 = LabelToTableMap.find( it->first );
    if( it1!=LabelToTableMap.end() )
    {
#if _OPENMP >= 200805L
#pragma omp task
#endif
      {
	row->DeepCopy( Table->GetRow( it1->second ) );
      }
      cropTable->InsertNextRow( row );    
    }
    else
    {
      std::cerr << it->second << " has no table entry. Nucleus Editor will crash if you click on it.\n";
    }
  }

  std::map< itk::SizeValueType , itk::SizeValueType >::iterator it = LabelToRelabelMap.begin();
  CroppedCentroidsMap.clear();
  CroppedBoundBoxesMap.clear();
  
  //Should be replaced with iteration over LabelToRelabelMap in parallel
  for( vtkIdType i=0; i<cropTable->GetNumberOfRows(); ++i, ++it )
  {
#if _OPENMP >= 200805L
#pragma omp task
#endif
    {
      ftk::Object::Point CurrentCentroid, CurrentMin, CurrentMax;
      ftk::Object::Box CurrentBoundBox;
      cropTable->SetValue(i,0,i+1);
      std::map< itk::SizeValueType , itk::SizeValueType >::iterator it1;
      it1 = LabelToTableMap.find( it->first );
      itk::SizeValueType x = TableEntryVector.at( it1->second ).x;
      CurrentCentroid.x = CheckBoundsAndSubtractMin( x, x1, x2 );
      cropTable->SetValueByName( i, "centroid_x", x );
      itk::SizeValueType y = TableEntryVector.at( it1->second ).y;
      CurrentCentroid.y = CheckBoundsAndSubtractMin( y, y1, y2 );
      cropTable->SetValueByName( i, "centroid_y", y );
      CurrentCentroid.z = TableEntryVector.at( it1->second ).z;
      CroppedCentroidsMap[ it->second ] = CurrentCentroid;
      itk::SizeValueType XMin, YMin, ZMin, XMax, YMax, ZMax;
      XMin = BoundingBoxes.at(it1->second).at(0); XMax = BoundingBoxes.at(it1->second).at(1);
      YMin = BoundingBoxes.at(it1->second).at(2); YMax = BoundingBoxes.at(it1->second).at(3);
      CurrentMin.z = BoundingBoxes.at(it1->second).at(4);
      CurrentMax.z = BoundingBoxes.at(it1->second).at(5);
      CurrentMin.x = CheckBoundsAndSubtractMin( XMin, x1, x2 );
      CurrentMax.x = CheckBoundsAndSubtractMin( XMax, x1, x2 );
      CurrentMin.y = CheckBoundsAndSubtractMin( YMin, y1, y2 );
      CurrentMax.y = CheckBoundsAndSubtractMin( YMax, y1, y2 );
      CurrentBoundBox.min = CurrentMin; CurrentBoundBox.max = CurrentMax;
      CroppedBoundBoxesMap[ it->second ] = CurrentBoundBox;
    }
  }
#if DEBUG_MONTV
    ftk::SaveTable( "crop_table1.txt", cropTable );
#endif
  return cropTable;
}

itk::SizeValueType MontageView::CheckBoundsAndSubtractMin( itk::SizeValueType CoOrd,
					itk::SizeValueType Min, itk::SizeValueType Max )
{
  if( CoOrd<Min ) CoOrd = 0;
  else if( CoOrd>Max ) CoOrd = Max-Min-1;
  else CoOrd -= Min;
  return CoOrd;
}
//Take the maximum intensity projection and subsample each channel
//Subsample the image setting the maximum x-y dimension to 1000pixels
//Redo the contrast, cast, smooth(for aliasing), then resample
void MontageView::resetSubsampledImageAndDisplayImage()
{
  typedef itk::Image< unsigned char, 3 > Uchar3DImageType;
  typedef itk::Image< float, 3 > Float3DImageType;
  typedef itk::MaximumProjectionImageFilter< Uchar3DImageType, Uchar3DImageType > MaxProjectType;
//typedef itk::AdaptiveHistogramEqualizationImageFilter< Uchar3DImageType > AdaptiveHistEqType;
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
      currentChannelProjection = currentChannel;

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

itk::SizeValueType MontageView::InsertNewLabelToRelabelMap( itk::SizeValueType NewKey )
{
  if( LabelToRelabelMap.empty() )
  {
    LabelToRelabelMap.insert( std::map<itk::SizeValueType, itk::SizeValueType>::value_type
				( NewKey, 1 ) );
    return 1;
  }
  //Just in case another thread has inserted a key
#ifdef _OPENMP
  std::map< itk::SizeValueType , itk::SizeValueType >::iterator it;
  it = LabelToRelabelMap.find( NewKey );
  if( it!=LabelToRelabelMap.end() )
    return it->second;
#endif
  LabelToRelabelMap.insert( std::map<itk::SizeValueType, itk::SizeValueType>::value_type
    				( NewKey, (LabelToRelabelMap.size()+1) ) );
  return LabelToRelabelMap.size();
}
