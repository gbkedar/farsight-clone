#ifndef ADAPTIVE_BINARIZATION_H
#define ADAPTIVE_BINARIZATION_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<limits>
#include<exception>
#include <algorithm>

#include "new_graph.h"
#include "itkMinErrorThresholdImageFilter.h"


#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
////////////////////

typedef unsigned short inputPixelType;
typedef unsigned char UCHARPixelType;
typedef itk::Image< inputPixelType, 3 >         inputImageType;
typedef itk::Image< UCHARPixelType, 3 >         binaryImageType;
typedef itk::ImageFileReader<inputImageType> ReaderType;
typedef itk::ImageFileWriter<binaryImageType> WriterType;

double computePoissonProb(int intensity, double alpha);
binaryImageType::Pointer Adaptive_Binarization(inputImageType::Pointer _inputImage);
//void subtractGradientImage(unsigned char* IM, int r, int c, int z, int sampl_ratio);

typedef itk::Image< unsigned short, 2 > US2ImageType;
typedef itk::Image< unsigned short, 3 > US3ImageType;
typedef itk::Image< unsigned char, 3 > UC3ImageType;
typedef itk::Image< double, 2 > CostImageType;
template<typename InputImageType> typename InputImageType::Pointer
  CreateDefaultCoordsNAllocateSpace( typename InputImageType::SizeType size );
template<typename InputImageType, typename OutputImageType> void
  RescaleCastNWriteImage( typename InputImageType::Pointer inputImage,
    std::string &outFileName );
template<typename InputImageType> void WriteITKImage
  ( typename InputImageType::Pointer inputImagePointer, std::string outputName );
void ComputeHistogram(
	US2ImageType::Pointer medFiltImages,
	std::vector< double > &histogram,
	US2ImageType::IndexType &start, US2ImageType::PixelType valsPerBin );
void UpdateHistogram(
	US2ImageType::Pointer medFiltImages, std::vector< double > &histogram,
	US2ImageType::IndexType &start, US2ImageType::IndexType &prevStart,
	US3ImageType::PixelType valsPerBin );
double ComputePoissonProbability(double intensity, double alpha);
void ComputeCut( US2ImageType::Pointer  medFiltImages,
		 CostImageType::Pointer flourCosts,
		 UC3ImageType::Pointer outputImage,
		 UC3ImageType::PixelType foregroundValue
		);
void ComputeCosts( US2ImageType::Pointer  medFiltImages,
		   CostImageType::Pointer flourCosts,
		   CostImageType::Pointer flourCostsBG,
		   US2ImageType::PixelType valsPerBin
		   );
UC3ImageType::Pointer AdaptiveBinarization2D( US3ImageType::Pointer inputImage3d );
US2ImageType::Pointer Get2DFrom3D( US3ImageType::Pointer readImage );

#endif
