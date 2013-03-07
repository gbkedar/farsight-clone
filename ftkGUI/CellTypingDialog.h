#ifndef _CELL_TYPING_DIALOG_H_
#define _CELL_TYPING_DIALOG_H_

//QT INCLUDES
#include <QtGui/QDialog>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QButtonGroup>
#include <QtGui/QVBoxLayout>
#include <QtGui/QSlider>
#include <QtGui/QSpinBox>
#include <QtGui/QGroupBox>
#include <QtGui/QCheckBox>
#include <QtGui/QLineEdit>
#include <QSignalMapper>
#include <QLabel>
#include <QtCore/QMap> 

//VTK INCLUDES
#include <vtkTable.h>

//ITK Includes
#include <itkIntTypes.h> 

#include "ftkCommon/ftkUtils.h"

#include "ClusClus/Kmeans.h"

#include <iostream>
#include <vector>
#include <map>

class featureSpinBox;
class featureSlider;

class CellTypingDialog : public QDialog
{
  Q_OBJECT;

public:
  CellTypingDialog( std::vector<unsigned>&, vtkSmartPointer<vtkTable>,
  	std::vector< std::vector<std::string> >&, std::vector<std::string>&,
	std::vector< std::pair< itk::SizeValueType, itk::SizeValueType > >&, QWidget *parent ); 

signals:
  void thresholdChanged( std::vector< std::pair<itk::SizeValueType,itk::SizeValueType> >, QColor,
  			 unsigned );
  void closing(void);

private slots:
  void featureSliderChange(unsigned, unsigned, unsigned, bool);
  void ResetToDefault(unsigned);

protected:
  virtual void closeEvent(QCloseEvent *event);

private:
  void SetSlidersToDefault(void);

  QButtonGroup *clusteringGroup;
  QSignalMapper *resetSignalMapper;
  QVBoxLayout *masterLayout;
  QPushButton *quitButton;
  QPushButton *doneButton;
  std::vector< std::vector< featureSlider* > > classThresholdSliderVector;
  std::vector< std::vector< featureSpinBox* > > classThresholdSpinBoxVector;
  std::vector< QPushButton * > resetButtons;

  QVector<QColor> centroidColorOptions;
  vtkSmartPointer<vtkTable> table;
  std::vector< std::vector<std::string> > clusteringGroups;
  std::vector<std::string> clusteringGroupNames;
  std::vector<unsigned> displayChannelNumbers;
  std::vector< std::pair<itk::SizeValueType,itk::SizeValueType> > DownSampledCoords;
  std::vector< std::vector<double> > defaultThresholds;
  std::vector<std::string> outputColNames;
  std::vector< std::vector< std::vector< bool > > > featurePositiveBools;

};

class featureSpinBox : public QSpinBox
{
	Q_OBJECT;
public:
  featureSpinBox(unsigned i, unsigned j, QWidget *parent = 0)
  {
    ChannelIndex = i; AssociationIndex = j; 
    connect(this, SIGNAL(valueChanged(unsigned)), this, SLOT(featureSpinBoxChanged(unsigned)));
  };	
signals:
  void featureValueChanged(unsigned value, unsigned i, unsigned j, bool from_slider);

private slots:
  void featureSpinBoxChanged(unsigned value){ emit featureValueChanged(value, ChannelIndex, AssociationIndex, false); };

private:
  unsigned ChannelIndex, AssociationIndex;
};


class featureSlider : public QSlider
{
  Q_OBJECT;

public:
  featureSlider(unsigned i, unsigned j, QWidget *parent = 0)
  { 
    ChannelIndex = i; AssociationIndex = j; 
    connect(this, SIGNAL(valueChanged(unsigned)), this, SLOT(featureSliderChanged(unsigned)));
  };	
	
signals:
  void featureValueChanged(unsigned value, unsigned i, unsigned j, bool from_slider);

private slots:
  void featureSliderChanged(unsigned value){ emit featureValueChanged(value, ChannelIndex, AssociationIndex, true); };

private:
  unsigned ChannelIndex, AssociationIndex;
};

#endif
