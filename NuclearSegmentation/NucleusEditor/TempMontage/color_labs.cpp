#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabelToRGBImageFilter.h"

typedef itk::Image<unsigned int, 3> InputLabelImageType;
typedef itk::ImageFileReader<InputLabelImageType> ReaderType;
typedef itk::Image<itk::RGBPixel<unsigned char>, 3> RGBImageType;
typedef itk::LabelToRGBImageFilter< InputLabelImageType, RGBImageType > LabelToRGBFilterType;
typedef itk::ImageFileWriter<RGBImageType> WriterType;

int main(int argc, char *argv[])
{
  if( argc<3 )
  {
    std::cout << "USAGE:\n";
    std::cout << " " << argv[0] << " InputLabelImage OutputColorLabelImage\n";
  }

  std::string inputFilename = argv[1];
  std::string outputFilename = argv[2];

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputFilename.c_str() );

  LabelToRGBFilterType::Pointer labelToRGBFilter = LabelToRGBFilterType::New();
  labelToRGBFilter->SetInput( reader->GetOutput() );

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFilename.c_str() );
  writer->SetInput( labelToRGBFilter->GetOutput() );

  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Error in montage compute stats:"<< excp << std::endl;
    exit (EXIT_FAILURE);
  }

  return  EXIT_SUCCESS;

}
