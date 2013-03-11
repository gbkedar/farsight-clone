#include "CellTypingDialog.h"

CellTypingDialog::CellTypingDialog( std::vector<unsigned> &channelNums, vtkSmartPointer<vtkTable> tbl,
  	std::vector< std::vector<std::string> > &clusGroups, std::vector<std::string> &clusGroupsNames,
	std::vector< std::pair<itk::SizeValueType,itk::SizeValueType> > &DwnSmplCords, QWidget *parent = 0 )
	: QDialog(parent)
{
  this->setModal(false);

  table = tbl;
  clusteringGroups = clusGroups;
  DownSampledCoords = DwnSmplCords;
  clusteringGroupNames = clusGroupsNames;

  clusteringGroup = new QButtonGroup;
  clusteringGroup->setObjectName("clusteringGroup");
  clusteringGroup->setExclusive(true);
  masterLayout = new QVBoxLayout;
  
  classThresholdSliderVector.resize( clusteringGroupNames.size() );
  classThresholdSpinBoxVector.resize( clusteringGroupNames.size() );
  defaultThresholds.resize( clusteringGroupNames.size() );
  resetButtons.resize( clusteringGroupNames.size() );

  resetSignalMapper = new QSignalMapper(this);

  for( unsigned c=0; c<clusteringGroupNames.size(); ++c )
  {
    std::string currentGroupName = clusteringGroupNames.at(c);
    QGridLayout *channelLayout = new QGridLayout();

    QCheckBox *check = new QCheckBox( QString( currentGroupName.c_str() ) );
    check->setObjectName( QString( currentGroupName.c_str() ) );
    check->setChecked( false );
    clusteringGroup->addButton(check, c);
		
    channelLayout->setColumnMinimumWidth(0, 70);
    channelLayout->addWidget(check, 0, 0);

    QPushButton *resetButton = new QPushButton(tr("Reset"));
    resetButton->setDefault(false);
    resetButton->setAutoDefault(false);
    resetButtons.at(c) = resetButton;
    channelLayout->addWidget(resetButton, 0, 4);
    connect(resetButton, SIGNAL(clicked()), this, SLOT(map()));
    resetSignalMapper->setMapping( resetButton, c );

    classThresholdSliderVector.at(c).resize( clusteringGroups.at(c).size() );
    classThresholdSpinBoxVector.at(c).resize( clusteringGroups.at(c).size() );
    defaultThresholds.at(c).resize( clusteringGroups.at(c).size() );

    for( unsigned i=0; i<clusteringGroups.at(c).size(); ++i )
    {
      unsigned found = clusteringGroups.at(c).at(i).find_first_of("_");
      std::string assocName = clusteringGroups.at(c).at(i).substr(found+1);

      QLabel *assocLabel = new QLabel(QString(assocName.c_str()));
      assocLabel->setDisabled(false);
      channelLayout->setColumnMinimumWidth(1, 50);
      channelLayout->addWidget(assocLabel, i, 1);

      featureSlider *assocSlider = new featureSlider(c, i, this);
      assocSlider->setOrientation(Qt::Horizontal);
      assocSlider->setDisabled(true);
      assocSlider->setRange(0,0);
      assocSlider->setValue(0);
      classThresholdSliderVector.at(c).at(i) = assocSlider;
      channelLayout->addWidget(assocSlider, i, 2);
      channelLayout->setColumnMinimumWidth(2, 150);

      featureSpinBox *assocSpinBox = new featureSpinBox(c, i, this);
      assocSpinBox->setDisabled(true);
      assocSpinBox->setRange(0,0);
      assocSpinBox->setValue(0);
      classThresholdSpinBoxVector.at(c).at(i) = assocSpinBox;
      channelLayout->addWidget(assocSpinBox, i, 3);
      channelLayout->setColumnMinimumWidth(3, 60);
    }
    masterLayout->addLayout(channelLayout);
  }
  connect(resetSignalMapper, SIGNAL(mapped(int)), this, SLOT(ResetToDefault(int)));

  quitButton = new QPushButton(tr("Cancel"));
  connect(quitButton, SIGNAL(clicked()), this, SLOT(reject()));
  quitButton->setDefault(false);
  quitButton->setAutoDefault(false);

  doneButton = new QPushButton(tr("Done"));
  connect(doneButton, SIGNAL(clicked()), this, SLOT(accept()));
  doneButton->setDefault(false);
  doneButton->setAutoDefault(false);

  QHBoxLayout *endLayout = new QHBoxLayout;
  endLayout->addStretch(25);
  endLayout->addWidget(doneButton);
  endLayout->addWidget(quitButton);
  masterLayout->addLayout(endLayout);
  this->setLayout(masterLayout);

  this->setWindowTitle(tr("Cell Typing Dialog"));

  Qt::WindowFlags flags = this->windowFlags();
  flags &= ~Qt::WindowContextHelpButtonHint;
  this->setWindowFlags(flags);

  SetSlidersToDefault();

  for(unsigned i=0; i<classThresholdSliderVector.size(); ++i)
  {
    for(unsigned j=0; j<classThresholdSliderVector[i].size(); ++j)
    {
      connect( classThresholdSliderVector.at(i).at(j), SIGNAL(featureValueChanged(unsigned, unsigned, unsigned, bool)),
		this, SLOT(featureSliderChange(unsigned, unsigned, unsigned, bool)));
      connect( classThresholdSpinBoxVector.at(i).at(j), SIGNAL(featureValueChanged(unsigned, unsigned, unsigned, bool)),
		this, SLOT(featureSliderChange(unsigned, unsigned, unsigned, bool)));
    }
  }
}

void CellTypingDialog::SetSlidersToDefault(void)
{
  //######### Generate column names for the table for channel positivity values
  outputColNames.clear();
  for(unsigned i=0; i<clusteringGroupNames.size(); ++i)
  {
    std::string currentGroupName = clusteringGroupNames.at(i) + "_positive";
    outputColNames.push_back( currentGroupName );
  }

  //######### Compute the max and min values for each associative feature
  //######### Accessing the associative features of each channel
#ifdef _OPENMP
  #pragma omp parallel for
#if _OPENMP < 200805L
  for( itk::IndexValueType i=0; i<clusteringGroups.size(); ++i )
#else
  for(  itk::SizeValueType i=0; i<clusteringGroups.size(); ++i )
#endif
#else
  for(  itk::SizeValueType i=0; i<clusteringGroups.size(); ++i )
#endif
  {
    //######### Accessing the each associative feature of this channel
    for(unsigned j=0; j<clusteringGroups.at(i).size(); ++j)
    {
      double min, max, val;
      min = table->GetValueByName(0, clusteringGroups.at(i).at(j).c_str()).ToDouble();
      max = table->GetValueByName(0, clusteringGroups.at(i).at(j).c_str()).ToDouble();
      for(unsigned row=0; row<table->GetNumberOfRows(); ++row)
      {
	double val = table->GetValueByName(row, clusteringGroups.at(i).at(j).c_str()).ToDouble();
	if( val > max ) max = val;
	else if( val < min ) min = val;
      }
      classThresholdSliderVector.at(i).at(j)->setRange( floor(min), ceil(max) );
      classThresholdSpinBoxVector.at(i).at(j)->setRange( floor(min), ceil(max) );
    }
  }
  featurePositiveBools.resize(clusteringGroups.size());

  //######### Accessing the associative features of each channel
#ifdef _OPENMP
  #pragma omp parallel for
#if _OPENMP < 200805L
  for( itk::IndexValueType i=0; i<clusteringGroups.size(); ++i )
#else
  for(  itk::SizeValueType i=0; i<clusteringGroups.size(); ++i )
#endif
#else
  for(  itk::SizeValueType i=0; i<clusteringGroups.size(); ++i )
#endif
  {
    //######### Accessing each associative feature of this channel
    featurePositiveBools.at(i).resize(clusteringGroups.at(i).size());
    for(unsigned j=0; j<clusteringGroups.at(i).size(); ++j)
    {
      featurePositiveBools.at(i).at(j).resize(table->GetNumberOfRows());
      std::vector<std::vector<double > > featureMatrix;
      //######### Run k-means on each feature, get the cluster centers and assign the average value of the
      //######### cluster centers to the default slider value
      for(unsigned row=0; row<table->GetNumberOfRows(); ++row)
      {
	std::vector<double> dataPoint;
	dataPoint.push_back(table->GetValueByName(row, clusteringGroups.at(i).at(j).c_str()).ToDouble());
	featureMatrix.push_back(dataPoint);
      }
      
      Kmeans *kmeans = new Kmeans();
      kmeans->setDataToKmeans(featureMatrix);
      kmeans->setClusterNumber(2);
      kmeans->Clustering();
      std::vector< std::vector <double> > clusterCenters = kmeans->getClusterCenters();
      double sliderDefaultValue = ( clusterCenters.at(0).at(0) + clusterCenters.at(1).at(0) )/2;
      defaultThresholds.at(i).at(j) = sliderDefaultValue;
      classThresholdSliderVector.at(i).at(j)->setValue(sliderDefaultValue);
      classThresholdSliderVector.at(i).at(j)->setEnabled(true);
      classThresholdSpinBoxVector.at(i).at(j)->setValue(sliderDefaultValue);
      classThresholdSpinBoxVector.at(i).at(j)->setEnabled(true);
      
      //######### Determine wether the data point is positive for this associative feature w.r.t to the default slider value
      for(unsigned row=0; row<(unsigned)table->GetNumberOfRows(); ++row)
      {
	double val = table->GetValueByName(row, clusteringGroups.at(i).at(j).c_str()).ToDouble();
	featurePositiveBools.at(i).at(j).at(row) = (val >= sliderDefaultValue);
      }
    }
  }
}

void CellTypingDialog::ResetToDefault(unsigned groupNumber)
{
  for(unsigned i=0; i<defaultThresholds.size(); ++i)
  {
    for(unsigned j=0; j<defaultThresholds[i].size(); ++j)
    {
      classThresholdSliderVector[i][j]->setValue(defaultThresholds[i][j]);
      classThresholdSpinBoxVector[i][j]->setValue(defaultThresholds[i][j]);
    }
  }
}

void CellTypingDialog::featureSliderChange(unsigned val, unsigned ch, unsigned assoc, bool from_slider)
{
/*	std::vector< bool > channel_pos_bools;
	channel_pos_bools.resize((unsigned)table->GetNumberOfRows());
	double hreshold_value;
	if(from_slider)
	{
		threshold_value = classThresholdSliderVector[ch][assoc]->value();
		classThresholdSpinBoxVector[ch][assoc]->setValue(threshold_value);
	}
	else
	{
		threshold_value = classThresholdSpinBoxVector[ch][assoc]->value();
		classThresholdSliderVector[ch][assoc]->setValue(threshold_value);
	}

	for(unsigned row=0; row<(unsigned)table->GetNumberOfRows(); ++row)
	{
		double val = table->GetValueByName(row, outputColNames[ch][assoc].c_str()).ToDouble();
		featurePositiveBools[ch][assoc][row] = (val >= threshold_value);

		bool channel_positive = true;
		for(unsigned j=0; j<outputColNames[ch].size(); ++j)
		{
			channel_positive = (channel_positive && featurePositiveBools[ch][j][row]);
		}
		channel_pos_bools[row] = channel_positive;
		if(channel_positive)
			table->SetValueByName(row, channelNames[ch].c_str(), vtkVariant(1));
		else
			table->SetValueByName(row, channelNames[ch].c_str(), vtkVariant(0));
	}

	emit thresholdChanged(channel_pos_bools);
*/
}

void CellTypingDialog::closeEvent(QCloseEvent *event)
{
	emit closing();
}
