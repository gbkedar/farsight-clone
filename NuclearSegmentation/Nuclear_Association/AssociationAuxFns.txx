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
#ifndef _ASC_FEAT_AUX_FN_txx_
#define _ASC_FEAT_AUX_FN_txx_

#define _USE_MATH_DEFINES
#include <math.h>

#include <algorithm>

#include "itkImage.h"
#include "itkIntTypes.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"
#include "itkExtractImageFilter.h"
#include "itkOtsuMultipleThresholdsCalculator.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_real.h"
#include "vnl/algo/vnl_real_eigensystem.h"
#include "vnl/vnl_double_3x3.h"

#include "ftkNuclearAssociationRules.h"

template < typename InputImageType , typename LabelImageType >
std::vector<float> compute_ec_features( typename itk::SmartPointer<InputImageType> input_image,
					typename itk::SmartPointer<LabelImageType> inp_labeled,
					int number_of_rois, typename InputImageType::PixelType thresh, int surr_dist,
					int inside_dist )
{
	std::vector< float > qfied_num;
	std::vector< typename LabelImageType::PixelType > labelsList;

	typedef itk::Image< float, 3 > FloatImageType;
	typedef itk::ImageRegionIteratorWithIndex< LabelImageType > IteratorType;
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > InpIteratorType;
	typedef itk::SignedDanielssonDistanceMapImageFilter< FloatImageType, FloatImageType > DTFilter;
	typedef itk::ImageRegionIteratorWithIndex< FloatImageType > IteratorTypeFloat;
	typedef itk::LabelStatisticsImageFilter< InputImageType, LabelImageType > LabelStatisticsImageFilterType;
	typedef typename LabelStatisticsImageFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
	typedef itk::ImageRegionConstIterator< FloatImageType > ConstIteratorTypeFloat;

	typename LabelImageType::SizeType sizee;
	sizee[0] = input_image->GetLargestPossibleRegion().GetSize()[0];
	sizee[1] = input_image->GetLargestPossibleRegion().GetSize()[1];
	sizee[2] = input_image->GetLargestPossibleRegion().GetSize()[2];
	bool image_is_3d = true;
	if( sizee[2] == 1 ) image_is_3d = false;

	typename LabelStatisticsImageFilterType::Pointer LabelStatisticsImageFilter =
						LabelStatisticsImageFilterType::New();
	LabelStatisticsImageFilter->SetLabelInput( inp_labeled );
	LabelStatisticsImageFilter->SetInput( input_image );
	LabelStatisticsImageFilter->UseHistogramsOff();
	try{ LabelStatisticsImageFilter->Update(); }
	catch( itk::ExceptionObject & excp ){ std::cerr << excp << std::endl; }

	//Store the label indices and centroids
	std::cout<<"Getting valid labels\n";
	for( typename ValidLabelValuesType::const_iterator
		vIt=LabelStatisticsImageFilter->GetValidLabelValues().begin();
		vIt != LabelStatisticsImageFilter->GetValidLabelValues().end(); ++vIt )
	{
		typename LabelImageType::PixelType LabelIndex=(*vIt);
		if( !LabelStatisticsImageFilter->HasLabel(*vIt) ) continue;
		labelsList.push_back( LabelIndex );
	}
	std::sort( labelsList.begin(), labelsList.end() );

	std::cout<<std::endl<<"The number of labels present: "<<labelsList.size()<<std::endl;

	itk::SizeValueType roi_list_size =
				((itk::SizeValueType)number_of_rois*labelsList.size());
	std::vector<double> quantified_numbers_cell( roi_list_size, 0.0 );
	std::cout<<"Bounding boxes computed"<<std::endl;

#ifdef _OPENMP
	int n_thr = 0.95*omp_get_max_threads();  //Use 95% of the cores by default n save a little for the OS
	std::cout<<"Number of threads "<<n_thr<<std::endl;
	itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);
	#pragma omp parallel for num_threads(n_thr)
#if _OPENMP >= 200805L
	for( typename LabelImageType::PixelType i=0; i<labelsList.size(); ++i )
#else
	for( itk::IndexValueType i=0; i<labelsList.size(); ++i )
#endif
	for( typename LabelImageType::PixelType i=0; i<labelsList.size(); ++i )
#endif
	{
		typename LabelImageType::PixelType LabelIndex = labelsList.at(i);
		typename LabelStatisticsImageFilterType::BoundingBoxType BoundBox =
				LabelStatisticsImageFilter->GetBoundingBox( labelsList.at(i) ); 
		typename LabelImageType::IndexType Start;
		Start[0] = BoundBox[0]; Start[1] = BoundBox[2]; Start[2] = BoundBox[4];
		typename LabelImageType::SizeType Size;
		Size[0] = BoundBox[1]-BoundBox[0]+1;
		Size[1] = BoundBox[3]-BoundBox[2]+1;
		Size[2] = image_is_3d ? (BoundBox[5]-BoundBox[4]+1) : 1;
		typename LabelImageType::RegionType CroppedRegion;
		CroppedRegion.SetSize ( Size  );
		CroppedRegion.SetIndex( Start );
		IteratorType PixBuf( inp_labeled, CroppedRegion );

		//Store the indices and the centroids
		double centroid_x, centroid_y, centroid_z;
		std::vector< typename LabelImageType::IndexType > LabelIndices;
		itk::SizeValueType pixel_count=0;
		for ( PixBuf.GoToBegin(); !PixBuf.IsAtEnd(); ++PixBuf )
		{
			if( PixBuf.Get() == LabelIndex )
			{
				typename LabelImageType::IndexType CurIn = PixBuf.GetIndex();
				centroid_x += CurIn[0];
				centroid_y += CurIn[1];
				centroid_z += CurIn[2];
				++pixel_count;
				LabelIndices.push_back( CurIn );				
			}
		}
		centroid_x = floor((centroid_x/(double)pixel_count)+0.5);
		centroid_y = floor((centroid_y/(double)pixel_count)+0.5);
		centroid_z = image_is_3d ? floor((centroid_z/(double)pixel_count)+0.5) : 0;

		//Create vnl array 3xN( label indicies )
		vnl_matrix<double> B( 3, pixel_count );

		FloatImageType::Pointer inp_lab_float = FloatImageType::New();
		FloatImageType::PointType origint;
		origint[0] = 0; origint[1] = 0; origint[2] = 0;
		inp_lab_float->SetOrigin( origint );
		FloatImageType::IndexType startt;
		startt[0] = 0;  // first index on X
		startt[1] = 0;  // first index on Y
		startt[2] = 0;  // first index on Z
		FloatImageType::SizeType  sizet;
		sizet[0] = BoundBox[1]-BoundBox[0]+2*surr_dist+2;  // size along X
		sizet[1] = BoundBox[3]-BoundBox[2]+2*surr_dist+2;  // size along Y
		sizet[2] = image_is_3d ? (BoundBox[5]-BoundBox[4]+2*surr_dist+2) : 1;  // size along Z
		FloatImageType::RegionType regiont;
		regiont.SetSize( sizet );
		regiont.SetIndex( startt );
		inp_lab_float->SetRegions( regiont );
		inp_lab_float->Allocate();
		inp_lab_float->FillBuffer(0.0);
		inp_lab_float->Update();

		IteratorTypeFloat iterator444 ( inp_lab_float, inp_lab_float->GetRequestedRegion() );
		//Populate matrix with deviations from the centroid for principal axes and
		//at the same time set up distance-transform computation
		itk::SizeValueType ind1=0;
		for( typename std::vector< typename LabelImageType::IndexType >::iterator itPixind =
				LabelIndices.begin();
				itPixind!=LabelIndices.end(); ++itPixind )
		{
			B(0,(ind1)) = (*itPixind)[0]-centroid_x;
			B(1,(ind1)) = (*itPixind)[1]-centroid_y;
			B(2,(ind1)) = (*itPixind)[2]-centroid_z;
			FloatImageType::IndexType cur_in;
			cur_in[0] = (*itPixind)[0]-BoundBox[0]+1+surr_dist;
			cur_in[1] = (*itPixind)[1]-BoundBox[2]+1+surr_dist;
			cur_in[2] = image_is_3d ? ((*itPixind)[2]-BoundBox[4]+1+surr_dist) : 0;
			iterator444.SetIndex( cur_in );
			iterator444.Set( 255.0 );
			++ind1;
		}

		//Compute distance transform for the current object
		DTFilter::Pointer dt_obj= DTFilter::New() ;
		dt_obj->SetInput( inp_lab_float );
		dt_obj->SquaredDistanceOff();
		dt_obj->InsideIsPositiveOff();
		try
		{
			dt_obj->Update() ;
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "Error in Distance Transform: " << err << std::endl;
		}
		FloatImageType::Pointer dist_im = dt_obj->GetOutput();

		//Use KLT to compute pricipal axes
		if(image_is_3d)
		{
			vnl_matrix<double> B_transp((int)pixel_count,3);
			B_transp = B.transpose();
			vnl_matrix<double>  COV(3,3);
			COV = B * B_transp;
			double norm = 1.0/(double)pixel_count;
			COV = COV * norm;
			//Eigen decomposition
			vnl_real_eigensystem Eyegun( COV );
			vnl_matrix<vcl_complex<double> > EVals = Eyegun.D;
			double Eval1 = vnl_real(EVals)(0,0);
			double Eval2 = vnl_real(EVals)(1,1);
			double Eval3 = vnl_real(EVals)(2,2);
			vnl_double_3x3 EVectMat = Eyegun.Vreal;
			double V1[3],V2[3],EP_norm[3];
			if( Eval1 >= Eval3 && Eval2 >= Eval3 )
			{
				if( Eval1 >= Eval2 )
				{
					V1[0] = EVectMat(0,0); V1[1] = EVectMat(1,0); V1[2] = EVectMat(2,0);
					V2[0] = EVectMat(0,1); V2[1] = EVectMat(1,1); V2[2] = EVectMat(2,1);
				}
				else
				{
					V2[0] = EVectMat(0,0); V2[1] = EVectMat(1,0); V2[2] = EVectMat(2,0);
					V1[0] = EVectMat(0,1); V1[1] = EVectMat(1,1); V1[2] = EVectMat(2,1);
				}
			}
			else if( Eval1 >= Eval2 && Eval3 >= Eval2 )
			{
				if( Eval1 >= Eval3 )
				{
					V1[0] = EVectMat(0,0); V1[1] = EVectMat(1,0); V1[2] = EVectMat(2,0);
					V2[0] = EVectMat(0,2); V2[1] = EVectMat(1,2); V2[2] = EVectMat(2,2);
				}
				else
				{
					V2[0] = EVectMat(0,0); V2[1] = EVectMat(1,0); V2[2] = EVectMat(2,0);
					V1[0] = EVectMat(0,2); V1[1] = EVectMat(1,2); V1[2] = EVectMat(2,2);
				}
			}
			else
			{
				if( Eval2 >= Eval3 )
				{
					V1[0] = EVectMat(0,1); V1[1] = EVectMat(1,1); V1[2] = EVectMat(2,1);
					V2[0] = EVectMat(0,2); V2[1] = EVectMat(1,2); V2[2] = EVectMat(2,2);
				}
				else
				{
					V2[0] = EVectMat(0,1); V2[1] = EVectMat(1,1); V2[2] = EVectMat(2,1);
					V1[0] = EVectMat(0,2); V1[1] = EVectMat(1,2); V1[2] = EVectMat(2,2);
				}
			}
			double n_sum = sqrt( V1[0]*V1[0]+V1[1]*V1[1]+V1[2]*V1[2] );
			V1[0] /= n_sum; V1[1] /= n_sum; V1[2] /= n_sum;
			n_sum = sqrt( V2[0]*V2[0]+V2[1]*V2[1]+V2[2]*V2[2] );
			V2[0] /= n_sum; V2[1] /= n_sum; V2[2] /= n_sum;
			//Get the normal to the plane formed by the biggest two EVs
			EP_norm[0] = V1[1]*V2[2]-V1[2]*V2[1];
			EP_norm[1] = V1[2]*V2[0]-V1[0]*V2[2];
			EP_norm[2] = V1[0]*V2[1]-V1[1]*V2[0];
			//Reassign V2 so that it is orthogonal to both EP_norm and V1
			V2[0] = V1[1]*EP_norm[2]-V1[2]*EP_norm[1];
			V2[1] = V1[2]*EP_norm[0]-V1[0]*EP_norm[2];
			V2[2] = V1[0]*EP_norm[1]-V1[1]*EP_norm[0];
			//Now we have the point normal form; EP_norm is the normal and
			//centroid_x, centroid_y, centroid_z is the point
			//The equation to the plane is
			//EP_norm[0](x-centroid_x)+EP_norm[1](y-centroid_y)+EP_norm[2](z-centroid_z)=0
			double dee = (centroid_x*EP_norm[0]+centroid_y*EP_norm[1]+centroid_z*EP_norm[2])*(-1.00);

			//Iterate through and assign values to each region
			ConstIteratorTypeFloat pix_buf2( dist_im, dist_im->GetRequestedRegion() );
			InpIteratorType iterator44( input_image, input_image->GetRequestedRegion() );

			for ( pix_buf2.GoToBegin(); !pix_buf2.IsAtEnd(); ++pix_buf2 )
			{
				//Use pixels that are only within the defined radius from the nucleus
				double current_distance = pix_buf2.Get();
				if( (current_distance <= (double)surr_dist) && (current_distance>=(-1*inside_dist)) )
				{
					typename LabelImageType::IndexType cur_in;
					double n_vec[3];
					cur_in[0] = pix_buf2.GetIndex()[0]+BoundBox[0]-1-surr_dist;
					cur_in[1] = pix_buf2.GetIndex()[1]+BoundBox[2]-1-surr_dist;
					cur_in[2] = pix_buf2.GetIndex()[2]+BoundBox[4]-1-surr_dist;
					if( cur_in[0] < 0 || cur_in[1] < 0 || cur_in[2] < 0 )
						continue;
					if( cur_in[0] >= sizee[0] || cur_in[1] >= sizee[1] || cur_in[2] >= sizee[2] )
						continue;
					iterator44.SetIndex( cur_in );
					typename InputImageType::PixelType pixel_intensity;
					pixel_intensity = iterator44.Get();
					if( pixel_intensity < thresh )
						continue;

					//The projection of the point on the plane formed by the fist two major axes
					double xxx, yyy, zzz;
					xxx = cur_in[0] - EP_norm[0]*
						((EP_norm[0]*cur_in[0]+EP_norm[1]*cur_in[1]+EP_norm[2]*cur_in[2]+dee)
						/(EP_norm[0]*EP_norm[0]+EP_norm[1]*EP_norm[1]+EP_norm[2]*EP_norm[2]));
					yyy = cur_in[1] - EP_norm[1]*
						((EP_norm[0]*cur_in[0]+EP_norm[1]*cur_in[1]+EP_norm[2]*cur_in[2]+dee)
						/(EP_norm[0]*EP_norm[0]+EP_norm[1]*EP_norm[1]+EP_norm[2]*EP_norm[2]));
					zzz = cur_in[2] - EP_norm[2]*
						((EP_norm[0]*cur_in[0]+EP_norm[1]*cur_in[1]+EP_norm[2]*cur_in[2]+dee)
						/(EP_norm[0]*EP_norm[0]+EP_norm[1]*EP_norm[1]+EP_norm[2]*EP_norm[2]));
					//The vector from the centroid to the projected point
					n_vec[0] = centroid_x-xxx;
					n_vec[1] = centroid_y-yyy;
					n_vec[2] = centroid_z-zzz;
					n_sum = sqrt( n_vec[0]*n_vec[0] + n_vec[1]*n_vec[1] + n_vec[2]*n_vec[2] );
					n_vec[0] /= n_sum; n_vec[1] /= n_sum; n_vec[2] /= n_sum;
					//n_vec is the normalized vect in the direction of the projected point
					//V1 is the largest eigenvector
					//Get the dot and cross product between the two
					double doooot, crooos,fin_est_angle;
					doooot = n_vec[0]*V1[0]+n_vec[1]*V1[1]+n_vec[2]*V1[2];
					crooos = n_vec[0]*V2[0]+n_vec[1]*V2[1]+n_vec[2]*V2[2];

					fin_est_angle = atan2( crooos, doooot );
					typename LabelImageType::PixelType bin_num;
					//Compute bin num
					if( fin_est_angle<0 )
						fin_est_angle += (2*M_PI);
					bin_num = floor(fin_est_angle*((double)number_of_rois/2)/(2*M_PI));

					//Check which side of the plane the point lies on
					double v_norm = (cur_in[0]-centroid_x)*(cur_in[0]-centroid_x)
									+(cur_in[1]-centroid_y)*(cur_in[1]-centroid_y)
									+(cur_in[2]-centroid_z)*(cur_in[2]-centroid_z);
					v_norm = sqrt( v_norm );
					double doot   = (cur_in[0]-centroid_x)*EP_norm[0]/v_norm 
							+ (cur_in[1]-centroid_y)*EP_norm[1]/v_norm 
							+ (cur_in[2]-centroid_z)*EP_norm[2]/v_norm;

					if( doot<0 )
						bin_num += ((double)number_of_rois/2);
					quantified_numbers_cell.at((i*(2*number_of_rois)+bin_num)) += pixel_intensity;
				}
			}
		}
		else
		{
			ConstIteratorTypeFloat pix_buf2( dist_im, dist_im->GetRequestedRegion() );
			InpIteratorType iterator44( input_image, input_image->GetRequestedRegion() );

			for ( pix_buf2.GoToBegin(); !pix_buf2.IsAtEnd(); ++pix_buf2 )
			{
				//Use pixels that are only within the defined radius from the nucleus
				double current_distance = pix_buf2.Get();
				if( (current_distance <= (double)surr_dist) && (current_distance>=(-1*inside_dist)) )
				{
					typename LabelImageType::IndexType cur_in;
					double n_vec[3];
					cur_in[0] = pix_buf2.GetIndex()[0]+BoundBox[0]-1-surr_dist;
					cur_in[1] = pix_buf2.GetIndex()[1]+BoundBox[2]-1-surr_dist;
					cur_in[2] = 0;
					if( cur_in[0] < 0 || cur_in[1] < 0 || cur_in[2] < 0 )
						continue;
					if( cur_in[0] >= sizee[0] || cur_in[1] >= sizee[1] || cur_in[2] >= sizee[2] )
						continue;
					iterator44.SetIndex( cur_in );
					typename InputImageType::PixelType pixel_intensity;
					pixel_intensity = iterator44.Get();
					if( pixel_intensity < thresh )
						continue;
					double angle = atan2((centroid_y-cur_in[1]),fabs(centroid_x-cur_in[0]));
					if( angle<0 )
						angle += (2*M_PI);
					typename LabelImageType::PixelType bin_num = floor(angle*number_of_rois/(2*M_PI));
					quantified_numbers_cell.at(i*number_of_rois+bin_num) += pixel_intensity;
				}
			}

		}
	}
#ifdef _OPENMP
#ifndef _MSC_VER
omp_set_nested(0);
itk::MultiThreader::SetGlobalMaximumNumberOfThreads(n_thr);
#endif
#endif
	std::cout<<"Starting k-means\n";
	if( labelsList.size() == 1 )
	{
		qfied_num.clear();
		return qfied_num;
	}
	std::vector<double> quantified_numbers_cell_cpy(roi_list_size);
	std::copy(quantified_numbers_cell.begin(), quantified_numbers_cell.end(), quantified_numbers_cell_cpy.begin() );
	//Run k-means
	//Most of the code is adapted from mul/mbl/mbl_k_means.cxx
	std::sort(quantified_numbers_cell.begin(), quantified_numbers_cell.end());
	unsigned k = 2;
	//Vectors and matrices for k-means
	std::vector< typename LabelImageType::PixelType > partition( roi_list_size, 0 );
	std::vector< double > sums   ( k, 0.0 );
	std::vector< double > centers( k, 0.0 );
	typename std::vector< typename LabelImageType::PixelType > nNearest( k, 0 );
	//Use the elements that are evenly spaced to get the intial centers
	for( unsigned i=0; i<k; ++i ){
		double index = ((double)(i+1)*roi_list_size)/(k+1);
		centers.at(i) = quantified_numbers_cell.at((itk::SizeValueType)index);
		bool duplicated;
		if(i){
			if( centers.at((i-1)) == centers.at(i) ){
				duplicated = true;
				itk::SizeValueType ind=i+1;
				while( centers.at((i-1))==quantified_numbers_cell.at(ind) )
					++ind;
				centers.at(i) = quantified_numbers_cell.at(ind);
				sums.at(i)    = quantified_numbers_cell.at(ind);
			}
		}
		if( !duplicated )
			sums.at(i) = quantified_numbers_cell.at((i+1)/(k+1));
		++nNearest[i];
	}

	bool changed = true;
	while(changed){
		changed = false;
		for(itk::SizeValueType i=0; i<roi_list_size; ++i){
			unsigned bestCentre = 0;
			double bestDist = fabs((centers.at(0)-quantified_numbers_cell.at(i)));
			for(unsigned j=1; j<k; ++j){
				double dist = fabs((centers.at(j)-quantified_numbers_cell.at(i)));
				if( dist < bestDist ){
					bestDist = dist; bestCentre = j;
				}
			}
			sums[bestCentre] += quantified_numbers_cell.at(i);
			++ nNearest[bestCentre];
			if( bestCentre != partition.at(i) ){
				changed = true;
				partition.at(i) = bestCentre;
			}
		}
		for( unsigned j=0; j<k; ++j) centers.at(j)  = sums.at(j)/nNearest.at(j);
		for( unsigned j=0; j<k; ++j) sums.at(j)     = 0;
		for( unsigned j=0; j<k; ++j) nNearest.at(j) = 0;
	}
	for( unsigned i=0; i<k; ++i )
		std::cout<<"Center "<<i<<" "<<centers.at(i)<<"\n";

	double Positive_thresh = ((centers.at(0)+centers.at(1))/2) < (255.0*thresh)?
				 ((centers.at(0)+centers.at(1))/2) : (255.0*thresh); //Arbitrary upper thresh

	std::cout<<"Positive_thresh "<<Positive_thresh<<"\n";

	std::cout<<"Done k-means\n";
	itk::SizeValueType ind = 0;
	for( typename LabelImageType::PixelType i=0; i<labelsList.size(); ++i ){
		int num_positive_rois = 0;
		for( unsigned j=0; j<number_of_rois; ++j ){
			itk::SizeValueType index_of_roi = ind*number_of_rois+j;
			if( quantified_numbers_cell_cpy.at(index_of_roi)>Positive_thresh )
				++num_positive_rois;
		}
		qfied_num.push_back(num_positive_rois);
		++ind;
	}
	std::cout<<"Done surroundedness\n"<<std::flush;
	return qfied_num;
}

template <class InputImageType>
typename InputImageType::PixelType returnthresh( typename itk::SmartPointer<InputImageType> input_image,
			     int num_bin_levs, int num_in_fg ){
	//Instantiate the different image and filter types that will be used
	typedef typename itk::ImageRegionConstIterator< InputImageType > ConstIteratorType;
	typedef itk::Statistics::Histogram< float > HistogramType;
	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;

	std::cout<<"Starting threshold computation\n";

	//Create a temporary histogram container:
	const int numBins = itk::NumericTraits<typename InputImageType::PixelType>::max();
	double *tempHist;
	tempHist = (double*) malloc( sizeof(double) * numBins );
	for(typename InputImageType::PixelType i=0; i<numBins; ++i)
		tempHist[i] = 0;

	typename InputImageType::PixelType maxval = itk::NumericTraits<typename InputImageType::PixelType>::ZeroValue();
	//Populate the histogram (assume pixel type is actually is some integer type):
	ConstIteratorType it( input_image, input_image->GetRequestedRegion() );
	for ( it.GoToBegin(); !it.IsAtEnd(); ++it ){
		typename InputImageType::PixelType pix = it.Get();
		++tempHist[pix];
		if( pix > maxval ) maxval = pix;
	}
	const int numBinsPresent = maxval+1;
	
	//Find max value in the histogram
	double floatIntegerMax = itk::NumericTraits<typename InputImageType::PixelType>::max();
	double max = 0.0;
	for(typename InputImageType::PixelType i=0; i<numBinsPresent; ++i)
		if( tempHist[i] > max )
			max = tempHist[i];

	double scaleFactor = 1;
	if(max >= floatIntegerMax)
		scaleFactor = floatIntegerMax / max;

	HistogramType::Pointer histogram = HistogramType::New() ;
	// initialize histogram
	HistogramType::SizeType size;
	HistogramType::MeasurementVectorType lowerBound;
	HistogramType::MeasurementVectorType upperBound;

	lowerBound.SetSize(1);
	upperBound.SetSize(1);
	size.SetSize(1);

	lowerBound.Fill(0.0);
	upperBound.Fill((double)maxval);
	size.Fill(numBinsPresent);

	histogram->SetMeasurementVectorSize(1);
	histogram->Initialize(size, lowerBound, upperBound ) ;

	int i=0;
	for (HistogramType::Iterator iter = histogram->Begin(); iter != histogram->End(); ++iter ){
		float norm_freq = (float)(tempHist[i] * scaleFactor);
		iter.SetFrequency(norm_freq);
		++i;
	}

	std::cout<<"Histogram computed\n";

	CalculatorType::Pointer calculator = CalculatorType::New();
	calculator->SetNumberOfThresholds( num_bin_levs );
	calculator->SetInputHistogram( histogram );
	calculator->Update();
	const CalculatorType::OutputType &thresholdVector = calculator->GetOutput(); 
	CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();

	float thresh;

	for(typename InputImageType::PixelType i=0; i < num_in_fg; ++itNum, ++i)
		thresh = (static_cast<float>(*itNum));

	std::cout<<"Threshold computed: "<<thresh<<std::endl;

	return (typename InputImageType::PixelType)(thresh+0.5);
}

#endif
