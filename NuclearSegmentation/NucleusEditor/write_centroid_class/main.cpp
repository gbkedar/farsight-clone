#include <ftkUtils.h>

#include "itkIntTypes.h"
#include "itkMultiThreader.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include <ftkImage.h>

typedef unsigned char PixelType;

int main(int argc, char *argv[])
{
	typedef itk::Image< unsigned char, 3 > InputImageType;
	typedef itk::RGBPixel< unsigned char > RGBPixelType;
	typedef itk::Image< RGBPixelType, 3 > RGBImageType;
	typedef itk::ImageFileReader< InputImageType > ReaderType;
	typedef itk::ImageFileWriter< InputImageType > WriterType;
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > CentroidIteratorType;
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > RGBIteratorType;

	if( argc<4 )
	{
		std::cout<<"Usage:"<<argv[0]<<" table sample_image label_image class_column_name1 class_column_name2 ...\n";
		return EXIT_FAILURE;
	}
	std::string tableFilename = argv[1];
	std::string sampleFilename = argv[2];

	std::vector< std::string > ClassCols; 
	for( int i=3; i<argc; ++i )
	{
		std::string currentCol = argv[i];
		ClassCols.push_back( currentCol );
	}

	vtkSmartPointer<vtkTable> table = NULL;
	if( ftk::FileExists(tableFilename) )
	{
		table = ftk::LoadTable(tableFilename);
	}
	else
	{
		std::cout<<"Could not load table\n";
		return EXIT_FAILURE;
	}

	//Try to load the input image:
	ftk::Image::Pointer myImg = NULL;
	if( ftk::GetExtension(inputFilename) == "xml" )
	{
		myImg = ftk::LoadXMLImage(inputFilename);
	}
	else
	{
		myImg = ftk::Image::New();
		if( !myImg->LoadFile(inputFilename) )
		{
			std::cout<<"Could not load input image\n";
			return EXIT_FAILURE;
		}
	}

	//Compose an itk-rgb image with the channels and the colors
	const ftk::Image::Info *info;
	info = channelImg->GetImageInfo();
	RGBImageType::SizeType  size;
	size[0] = sampleImage->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
	size[1] = sampleImage->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
	size[2] = sampleImage->GetOutput()->GetLargestPossibleRegion().GetSize()[2];

	//Generating an itk-rgb image
	RGBImageType::Pointer image = RGBImageType::New();
	RGBImageType::IndexType start;
	start[0] = start[1] = start[2] = 0;
	RGBImageType::PointType origin;
	origin[0] = origin[1] = origin[2] = 0;
	image->SetOrigin( origin );
	RGBImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	image->SetRegions( region );
	image->Allocate();
	image->FillBuffer(0);
	image->Update();

	for( int i=0; i < (*info).numChannels; ++i )
	{
		InputImageType::Pointer currentCh = myImg->GetItkPtr< BinaryImageType::PixelType >( 0, i, ftk::Image::DEFAULT );
		i
	}

	InputImageType::SizeType  size;
	{
		ReaderType::Pointer sampleImage = ReaderType::New();
		sampleImage->SetFileName( sampleFilename.c_str() );
		try
		{
			sampleImage->Update();
		}
		catch( itk::ExceptionObject & excp )
		{
			std::cerr << excp << std::endl;
			return EXIT_FAILURE;
		}
	}
//	#pragma omp parallel for shared( table, size, ClassCols )
	for( unsigned i=0; i<ClassCols.size(); ++i )
	{
		std::cout<<"Generating class map for "<<ClassCols.at(i)<<std::endl;
		InputImageType::Pointer image = InputImageType::New();
		InputImageType::IndexType start;
		start[0] = start[1] = start[2] = 0;
		InputImageType::PointType origin;
		origin[0] = origin[1] = origin[2] = 0;
		image->SetOrigin( origin );
		InputImageType::RegionType region;
		region.SetSize( size );
		region.SetIndex( start );
		image->SetRegions( region );
		image->Allocate();
		image->FillBuffer(0);
		image->Update();
		vtkAbstractArray * arr = table->GetColumnByName( ClassCols.at(i).c_str() );
		vtkAbstractArray * cen_x = table->GetColumnByName( "centroid_x" );
		vtkAbstractArray * cen_y = table->GetColumnByName( "centroid_y" );
		vtkAbstractArray * cen_z = table->GetColumnByName( "centroid_z" );
		for( vtkIdType j = 0; j < table->GetNumberOfRows(); ++j )
		{
			if( arr->GetVariantValue(j).ToInt() )
			{
				InputImageType::IndexType centroid;
				centroid[0] = cen_x->GetVariantValue(j).ToLong();
				centroid[1] = cen_y->GetVariantValue(j).ToLong();
				centroid[2] = cen_z->GetVariantValue(j).ToLong();
				CentroidIteratorType iter( image, image->GetRequestedRegion() );
				InputImageType::IndexType startC, endC;
				startC[0] = (centroid[0]-7)>0 ? (centroid[0]-7):0;
				startC[1] = (centroid[1]-7)>0 ? (centroid[1]-7):0;
				startC[2] = (centroid[2]-7)>0 ? (centroid[2]-7):0;
				endC[0] = (centroid[0]+7)<size[0] ? (centroid[0]+7):(size[0]-1);
				endC[1] = (centroid[1]+7)<size[1] ? (centroid[1]+7):(size[1]-1);
				endC[2] = (centroid[2]+7)<size[2] ? (centroid[2]+7):(size[2]-1);
				if( !endC[2] ) endC[2]=1;
				for( itk::IndexValueType col = startC[0]; col<endC[0]; ++col )
				  for( itk::IndexValueType row = startC[1]; row<endC[1]; ++row )
				    for( itk::IndexValueType stack = startC[2]; stack<endC[2]; ++stack )
				    {
				    	InputImageType::IndexType Index;
					Index[0] = col; Index[1] = row; Index[2] = stack;
					iter.SetIndex( Index );
					iter.Set( 255 );
				    }
			}
		}
		std::cout<<"Writing file\n"<<std::flush;
		std::string myFilename = sampleFilename;
		unsigned extension = ftk::GetExtension(myFilename).size()+1;
		std::string::iterator it;
		it = myFilename.end() - extension;
		myFilename.erase(it, it+extension);
		std::string classImageFileName = myFilename + "_" + ClassCols.at(i) + ".nrrd";
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( classImageFileName.c_str() );
		writer->SetInput( image );
		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & excp )
		{
			std::cerr << excp << std::endl;
		}
	}

	return EXIT_SUCCESS;
}
