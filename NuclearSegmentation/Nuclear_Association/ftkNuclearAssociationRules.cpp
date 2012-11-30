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
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#include "ftkNuclearAssociationRules.h"

#include <iostream>
#include <fstream>
#include <iomanip>

//just for testing

namespace ftk 
{	

AssociativeFeatureCalculator::AssociativeFeatureCalculator()
{
	inFilename = "";
	fPrefix = "";
	inputs_set = false;	
}

AssociativeFeatureCalculator::~AssociativeFeatureCalculator()
{
}

void AssociativeFeatureCalculator::SetInputFile(std::string filename)
{
	inFilename = filename;
}

void AssociativeFeatureCalculator::SetInputs(ftk::Image::Pointer inp_labeled_image, int inp_channel_number, ftk::Image::Pointer seg_labeled_image, int seg_channel_number, ftk::AssociationRule *associationrule){
	input_association = associationrule;
	if( seg_channel_number == -1 || inp_channel_number == -1 )
		return;

	LabelImage   = seg_labeled_image;
	LabelChannel = seg_channel_number;
	InputImage   = inp_labeled_image;
	InputChannel = inp_channel_number;

	inputs_set = true;
}

void AssociativeFeatureCalculator::SetFeaturePrefix(std::string prefix)
{
	fPrefix = prefix;
}

//Update the features in this table whose names match (sets doFeat)
void AssociativeFeatureCalculator::Append(vtkSmartPointer<vtkTable> table)
{
	typedef itk::Image< unsigned short, 3 > UShortImageType;
	typedef itk::Image< unsigned, 3 > UnsignedImageType;
	//Compute features:
	if( inputs_set ){
		const ftk::Image::Info *InputImageInfo, *LabelImageInfo;
		InputImageInfo = InputImage->GetImageInfo();
		LabelImageInfo = LabelImage->GetImageInfo();
		if( (InputImageInfo->dataType == itk::ImageIOBase::UCHAR) ||
		    (InputImageInfo->dataType == itk::ImageIOBase::USHORT) ){
			if( LabelImageInfo->dataType == itk::ImageIOBase::USHORT ){
				UShortImageType::Pointer inp_im = InputImage->GetItkPtr
					<unsigned short>(0,InputChannel,ftk::Image::DEEP_COPY);
				UShortImageType::Pointer lbl_im = LabelImage->GetItkPtr
					<unsigned short>(0,LabelChannel,ftk::Image::DEEP_COPY);
				this->Compute( table, inp_im, lbl_im );
				return;
			}
			if( LabelImageInfo->dataType == itk::ImageIOBase::UINT ){
				UShortImageType::Pointer inp_im = InputImage->GetItkPtr
					<unsigned short>(0,InputChannel,ftk::Image::DEEP_COPY);
				UnsignedImageType::Pointer lbl_im = LabelImage->GetItkPtr
					<unsigned>(0,LabelChannel,ftk::Image::DEEP_COPY);
				this->Compute( table, inp_im, lbl_im );
				return;
			}
		}
		else{
			std::cout<<"Unknown input datatype. Perhaps "<<
			"AssociativeFeatureCalculator::Append needs to be extended\n";
			return;
		}
	}
	std::cout<<"Inputs not set for associations\n";
	return;
}


} //end namespace ftk
