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

namespace ftk{

template< typename InputImageType, typename LabelImageType >
	void AssociativeFeatureCalculator::Compute( vtkSmartPointer<vtkTable> table,
	typename itk::SmartPointer< InputImageType > InputITKImage,
	typename itk::SmartPointer< LabelImageType > LabelITKImage ){

	ftk::NuclearAssociationRules< InputImageType, LabelImageType > *assoc;
	assoc = new ftk::NuclearAssociationRules< InputImageType, LabelImageType >
						("", 0, LabelITKImage, InputITKImage);


	assoc->AddAssociation( input_association->GetRuleName(), "",
			input_association->GetOutDistance(), input_association->GetInDistance(),
			input_association->IsUseWholeObject(), input_association->IsUseBackgroundSubtraction(),
			input_association->IsUseMultiLevelThresholding(), input_association->GetNumberOfThresholds(),
			input_association->GetNumberIncludedInForeground(), input_association->GetAssocType(),
			input_association->get_path() );
	assoc->PrintSelf();
	assoc->Compute();

	//Init the table (headers):
	for (int i=0; i < assoc->GetNumofAssocRules(); ++i)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( (fPrefix+assoc->GetAssociationRules().at(i).GetRuleName()).c_str() );
		column->SetNumberOfValues( table->GetNumberOfRows() );
		table->AddColumn(column);
	}

	//Now update the table:
	typename std::vector<typename LabelImageType::PixelType > labels = assoc->GetLabels();
	float** vals = assoc->GetAssocFeaturesList();
	//#pragma omp parallel for num_threads(4)
	for (long i=0; i<labels.size(); ++i)
	{
		typename LabelImageType::PixelType id = labels.at(i);
		if(id == 0) continue;

		long row = -1;
		for(typename LabelImageType::PixelType r=0; r<table->GetNumberOfRows(); ++r)
		{
			if( table->GetValue(r,0) == id )
			{
				row = r;
				break;
			}
		}

		if(row == -1) continue;

		for (int f=0; f<assoc->GetNumofAssocRules(); ++f)
		{
			table->SetValueByName(row,(fPrefix+assoc->GetAssociationRules().at(f).GetRuleName()).c_str(), vtkVariant(vals[f][i]));
		}
	}
	delete assoc;
	return;
}

template < typename InputImageType, typename LabelImageType >
NuclearAssociationRules< InputImageType, LabelImageType >::NuclearAssociationRules(std::string AssocFName, int numOfRules,
		typename itk::SmartPointer<LabelImageType> lImage,
		typename itk::SmartPointer<InputImageType> iImage):ObjectAssociation(AssocFName, numOfRules){
	labImage = lImage;
	inpImage = iImage;
	x_Size = y_Size = z_Size = 0;
	imDim=3;	
	objectType = "Nucleus";
	num_rois = 8;
}

template < typename InputImageType, typename LabelImageType >
NuclearAssociationRules< InputImageType, LabelImageType >::~NuclearAssociationRules()
{
}

/* This is the main function for computing associative features */
template < typename InputImageType, typename LabelImageType >
void NuclearAssociationRules< InputImageType, LabelImageType >::Compute()
{
	std::cout<<"Starting Associative Features Computation\n";

	x_Size=labImage->GetLargestPossibleRegion().GetSize()[0];
	y_Size=labImage->GetLargestPossibleRegion().GetSize()[1];
	z_Size=labImage->GetLargestPossibleRegion().GetSize()[2];
	if(z_Size < 2)
		imDim = 2;

	//2. Use the labelStatistics filter to get the bounding box for each cell
	LabelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
	LabelStatisticsImageFilter->SetLabelInput( labImage );
	LabelStatisticsImageFilter->SetInput( inpImage );
	LabelStatisticsImageFilter->UseHistogramsOff();
	try{ LabelStatisticsImageFilter->Update(); }
	catch( itk::ExceptionObject & excp ){ std::cerr << excp << std::endl; }

	const typename LabelImageType::PixelType NumLabels =
		LabelStatisticsImageFilter->GetNumberOfObjects();

	//Get the list of labels
	for( typename LabelImageType::PixelType LabelIndex=1; LabelIndex<NumLabels; ++LabelIndex )
		if( !LabelStatisticsImageFilter->GetCount(LabelIndex) ) continue;
		else LabelsList.push_back( LabelIndex );

	numLabels = (typename LabelImageType::PixelType)LabelsList.size();

	//allocate memory for the features list
	assocMeasurementsList = new float*[GetNumofAssocRules()];
	
	//3. then, for each type of the associations get the requested reigion based on the type and the value of the inside and outside distances.
	for(int i=0; i<GetNumofAssocRules(); i++)
	{
		assocMeasurementsList[i] = new float[numLabels];
		//read the ith target image (the image from which we need to compute the ith assoc. rule

		if( assocRulesList[i].IsUseBackgroundSubtraction() ){
			if( assocRulesList[i].IsUseMultiLevelThresholding() )
				if( assocRulesList[i].GetNumberOfThresholds() >= 
				    assocRulesList[i].GetNumberIncludedInForeground() )
					thresh=returnthresh( inpImage, assocRulesList[i].GetNumberOfThresholds(),
								assocRulesList[i].GetNumberIncludedInForeground() );
				else
					thresh=returnthresh( inpImage, assocRulesList[i].GetNumberOfThresholds(),
									assocRulesList[i].GetNumberOfThresholds() );
			else
				thresh=returnthresh( inpImage, 1, 1 );
			//Write Binary Mask
			/*std::string out_filename;
			out_filename = assocRulesList[i].get_path()+assocRulesList[i].GetRuleName();
			if( assocRulesList[i].GetAssocType() == ASSOC_SURROUNDEDNESS )
				out_filename = out_filename + "binary_surroundedness.tif";
			if( assocRulesList[i].GetAssocType() == ASSOC_MIN )
				out_filename = out_filename + "binary_min.tif";
			if( assocRulesList[i].GetAssocType() == ASSOC_MAX )
				out_filename = out_filename + "binary_max.tif";
			if( assocRulesList[i].GetAssocType() == ASSOC_TOTAL )
				out_filename = out_filename + "binary_total.tif";
			if( assocRulesList[i].GetAssocType() == ASSOC_AVERAGE )
				out_filename = out_filename + "binary_average.tif";
			BinaryThresholdType::Pointer threshfilt = BinaryThresholdType::New();
			threshfilt->SetInput( inpImage );
			threshfilt->SetLowerThreshold(thresh);
			threshfilt->SetUpperThreshold(USHRT_MAX);
			threshfilt->SetInsideValue( USHRT_MAX );
			threshfilt->SetOutsideValue( 0 );
			threshfilt->Update();
			WriterType::Pointer writer1 = WriterType::New();
			writer1->SetInput( threshfilt->GetOutput() );
			writer1->SetFileName( out_filename.c_str() );
			writer1->Update();*/
		}
		else
			thresh = 0; //Assuming unsigned datatypes
		if( assocRulesList[i].GetAssocType() == ASSOC_SURROUNDEDNESS ){
			int inside_distance = assocRulesList[i].GetInDistance();
			if( assocRulesList[i].IsUseWholeObject() )
				inside_distance = INT_MAX; //Larger than the radius of the cell
			std::vector<float> ec_feat_vals = compute_ec_features( inpImage, labImage, num_rois,
								thresh, assocRulesList[i].GetOutDistance(),
								inside_distance );
			for(typename LabelImageType::PixelType j=0; j<numLabels; ++j)
				assocMeasurementsList[i][j] = ec_feat_vals[j];
		} else {
			typename LabelImageType::PixelType counterLabels = 0;
#ifdef _OPENMP
#ifndef _MSC_VER
			int n_thr = 0.95*omp_get_max_threads();
			std::cout<<"Number of threads "<<n_thr<<std::endl;
			itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
			#pragma omp parallel for num_threads(n_thr) 
#endif
#endif
			for( typename LabelImageType::PixelType j=0; j<numLabels; ++j ) {
				typename LabelImageType::PixelType lbl = LabelsList[j];
				if(lbl == 0) continue;
				assocMeasurementsList[i][j] = ComputeOneAssocMeasurement(i, lbl);
			}
#ifdef _OPENMP
#ifndef _MSC_VER
			itk::MultiThreader::SetGlobalDefaultNumberOfThreads(n_thr);
#endif
#endif
		}
		std::cout<<"\tdone"<<std::endl;
	}
}

/* use this function to compute one associative measurement */
template < typename InputImageType, typename LabelImageType >
float NuclearAssociationRules< InputImageType, LabelImageType >::ComputeOneAssocMeasurement(int ruleID, typename LabelImageType::PixelType objID)
{	
	//fisrt, get the bounding box around the object
	//The bounding box is defined by the area around the object such that the object+outerdistance are included
	typename LabelStatisticsImageFilterType::BoundingBoxType bbox =
				LabelStatisticsImageFilter->GetBoundingBox( objID );
	FloatImageType::Pointer inp_lab = FloatImageType::New();
	FloatImageType::PointType origint; origint[0] = 0; origint[1] = 0; origint[2] = 0;
	inp_lab->SetOrigin( origint );

	//Iterate through the label image while keeping the bounds
	typename LabelImageType::IndexType startT;
	startT[0] = bbox[0]; startT[1] = bbox[2]; startT[2] = bbox[4];
	typename LabelImageType::SizeType sizeT;
	sizeT[0] = bbox[1]-bbox[0];
	sizeT[1] = bbox[3]-bbox[2];
	sizeT[2] = (imDim==2) ? 1 : (bbox[5]-bbox[4]) ;
	typename LabelImageType::RegionType regionT;
	regionT.SetSize( sizeT ); regionT.SetIndex( startT );

	//We don't really care about the bounds in distance transform image
	FloatImageType::IndexType start; start[0] = 0; start[1] = 0; start[2] = 0;
	FloatImageType::SizeType  size;
	size[0] = bbox[1]-bbox[0]+2*assocRulesList[ruleID].GetOutDistance()+2;
	size[1] = bbox[3]-bbox[2]+2*assocRulesList[ruleID].GetOutDistance()+2;
	size[2] = (imDim==2) ? 1 : (bbox[5]-bbox[4]+2*assocRulesList[ruleID].GetOutDistance()+2);
	FloatImageType::RegionType region;
	region.SetSize( size ); region.SetIndex( start );
	inp_lab->SetRegions( region );
	inp_lab->Allocate();
	inp_lab->FillBuffer(0.0);
	inp_lab->Update();

	IteratorTypeFloat iteratorDistInp ( inp_lab, inp_lab->GetRequestedRegion() );
	IteratorTypeLabel iteratorlab ( labImage, regionT );

	for( iteratorlab.GoToBegin(); !iteratorlab.IsAtEnd(); ++iteratorlab ){
		if( iteratorlab.Get()==objID ){
			typename LabelImageType::IndexType CurIn = iteratorlab.GetIndex();
			CurIn[0] += (assocRulesList[ruleID].GetOutDistance()+1-bbox[0]);
			CurIn[1] += (assocRulesList[ruleID].GetOutDistance()+1-bbox[2]);
			if( imDim == 3 )
				CurIn[2] += (assocRulesList[ruleID].GetOutDistance()+1-bbox[4]);
			iteratorDistInp.SetIndex( CurIn );
			iteratorDistInp.Set( 255.0 );
		}
	}

	//Compute the distance transform in the sub-segmentation image region
	DTFilter::Pointer dt_obj = DTFilter::New() ;
	dt_obj->SetInput( inp_lab ) ;
	dt_obj->SquaredDistanceOff();
	dt_obj->InsideIsPositiveOff();
	try{
		dt_obj->Update() ;
	}
	catch( itk::ExceptionObject & err ){
		std::cout << "Error in Distance Transform: " << err << std::endl; 
		return 0;	
	}

	//It is cheaper to compute all three Associtative measures at once
	double sum = 0;
	itk::SizeValueType pixcount = 0;
	float average=0, min=0, max=0;
	min = itk::NumericTraits<typename InputImageType::PixelType>::max();
	max = itk::NumericTraits<typename InputImageType::PixelType>::min();

	IteratorTypeFloat iterDistance ( dt_obj->GetOutput(), dt_obj->GetOutput()->GetRequestedRegion() );
	IteratorTypeInput iterGray     ( inpImage, inpImage->GetRequestedRegion() );

	for( iterDistance.GoToBegin(); !iterDistance.IsAtEnd(); ++iterDistance ){
		double current_distance = iterDistance.Get();
		if( ( current_distance <= (double)assocRulesList[ruleID].GetOutDistance() ) ||
		    ( ( current_distance >= (double)(-1.0*assocRulesList[ruleID].GetInDistance()) ) ||
		      ( assocRulesList[ruleID].IsUseWholeObject() && ( current_distance < 0 )     )
		    )
		  )
		{
			typename InputImageType::IndexType cur_in = iterDistance.GetIndex();
			cur_in[0] += (bbox[0]-1-assocRulesList[ruleID].GetOutDistance());
			cur_in[1] += (bbox[2]-1-assocRulesList[ruleID].GetOutDistance());
			if( imDim == 3 )
				cur_in[2] += (bbox[4]-1-assocRulesList[ruleID].GetOutDistance());
			if( cur_in[0] < 0 || cur_in[1] < 0 || cur_in[2] < 0 ) continue;
			if( cur_in[0] >= x_Size || cur_in[1] >= y_Size || cur_in[2] >= z_Size ) continue;
			iterGray.SetIndex( cur_in );
			typename InputImageType::PixelType pixel_intensity;
			pixel_intensity = iterGray.Get();
			if( pixel_intensity < thresh ) continue;
			if( pixel_intensity < min ) min = pixel_intensity;
			if( pixel_intensity > max ) max = pixel_intensity;
			sum += pixel_intensity;
			++pixcount;
		}
	}

	if( assocRulesList[ruleID].GetAssocType() == ASSOC_MIN )
		return min;

	if( assocRulesList[ruleID].GetAssocType() == ASSOC_MAX )
		return max;
	
	if( assocRulesList[ruleID].GetAssocType() == ASSOC_TOTAL )
		return (float)sum;

	if( sum>0 ) sum/=pixcount;
	return ( (float) sum );
}

}//end namespace ftk
