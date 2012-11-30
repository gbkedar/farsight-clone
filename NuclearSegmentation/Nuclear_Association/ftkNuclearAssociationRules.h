/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Module:    $RCSfile: ftkNuclearAssociationRules.h,v $
  Language:  C++
  Date:      $Date: 2008/11/27 1:00:00 $
  Version:   $Revision: 1 $
 
=========================================================================*/
#ifndef __ftkNuclearAssociationRules_h
#define __ftkNuclearAssociationRules_h

#include <limits.h>

#include <ftkObject.h>
#include <ftkFeatures/ftkObjectAssociation.h>
#include "ftkImage.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h>
#include <vtkTable.h>

#ifdef _OPENMP
#include "omp.h"
#endif

namespace ftk
{ 

class AssociativeFeatureCalculator
{
public:
	AssociativeFeatureCalculator();
	~AssociativeFeatureCalculator();
	void SetInputFile(std::string filename);
	void SetInputs(ftk::Image::Pointer inp_labeled_image, int inp_channel_number, ftk::Image::Pointer seg_labeled_image, int seg_channel_number, ftk::AssociationRule *associationrule);
	void SetFeaturePrefix(std::string prefix);			//Set Prefix for feature names
	void Append(vtkSmartPointer<vtkTable> table);		//Compute features that are ON and append them to the existing table

private:
	template< typename InputImageType, typename LabelImageType >
	void Compute( vtkSmartPointer<vtkTable> table,
		typename itk::SmartPointer< InputImageType > InputITKImage,
		typename itk::SmartPointer< LabelImageType > LabelITKImage );
	ftk::AssociationRule *input_association;
	std::string inFilename;
	std::string fPrefix;
	ftk::Image::Pointer LabelImage;
	int LabelChannel;
	ftk::Image::Pointer InputImage;
	int InputChannel;
	bool inputs_set;
};


/** \class NuclearAssociation
 *  \brief To define a set of association rules between nuclei and other objects and to compute the corresponding associative measures. 
 *  
 * \author Yousef Al-Kofahi, Rensselear Polytechnic Institute (RPI), Troy, NY USA
 */

typedef itk::Image< unsigned short, 3 > LabImageType;
typedef itk::Image< unsigned short, 3 > TargImageType;
template < typename InputImageType = TargImageType, typename LabelImageType = LabImageType >
class NuclearAssociationRules : public ObjectAssociation
{
public:
	/* Contsructor */
	NuclearAssociationRules(std::string AssocFName, int numOfRules,
		typename itk::SmartPointer<LabelImageType> lImage,
		typename itk::SmartPointer<InputImageType> iImage);
	~NuclearAssociationRules();

	/* This method computes all the associative measurements for all the objects */
	void Compute();

	/* Get the number of objects */
	typename LabelImageType::PixelType GetNumOfObjects() {return numLabels;};
	std::vector< typename LabelImageType::PixelType > GetLabels(){return LabelsList;};
private:
	/* Private member variables */
	typedef itk::Image< double, 3 > FloatImageType;
	typedef itk::SignedDanielssonDistanceMapImageFilter< FloatImageType, FloatImageType > DTFilter;
	typedef		 itk::ImageRegionIteratorWithIndex< FloatImageType > IteratorTypeFloat;
	typedef typename itk::ImageRegionIteratorWithIndex< LabelImageType > IteratorTypeLabel;
	typedef typename itk::ImageRegionIteratorWithIndex< InputImageType > IteratorTypeInput;
	typedef typename itk::LabelStatisticsImageFilter< InputImageType, LabelImageType > LabelStatisticsImageFilterType;
	typedef typename LabelStatisticsImageFilterType::ValidLabelValuesContainerType
									ValidLabelValuesType;

	typename itk::SmartPointer<LabelImageType> labImage;
	typename itk::SmartPointer<InputImageType> inpImage;
	itk::SizeValueType x_Size;
	itk::SizeValueType y_Size;
	itk::SizeValueType z_Size;
	itk::SizeValueType imDim;
	std::vector< typename LabelImageType::PixelType > LabelsList;
	typename LabelImageType::PixelType numLabels;

	typename InputImageType::PixelType thresh;
	typename LabelStatisticsImageFilterType::Pointer LabelStatisticsImageFilter;

private:
	/* This method is used to compute a single associative measurement for one cell */
	float ComputeOneAssocMeasurement(int ruleID, typename LabelImageType::PixelType objID);

	int num_rois;
}; // end NuclearAssociation
}  // end namespace ftk

#include "AssociationAuxFns.txx"
#include "ftkNuclearAssociationRules.txx"

#endif	// end __ftkNuclearAssociationRules_h
