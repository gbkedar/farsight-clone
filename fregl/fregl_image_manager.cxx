/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

#include <iostream>
#include <fregl/fregl_image_manager.h>
#include <fregl/fregl_joint_register.h>
#include <fregl/fregl_space_transformer.h>
#include <fregl/fregl_util.h>
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkTranslationTransform.h"
#include "itkResampleImageFilter.h"


typedef itk::ImageRegionIterator< ImageType > RegionIterator;
typedef itk::ImageRegionConstIterator< ImageType > RegionConstIterator;

static std::string ToString(double val);

fregl_image_manager::fregl_image_manager(std::string const & xml_filename, std::string const & image_path, std::string const & anchor_image,
        bool use_NN) {
    //    std::cout << "Anchor Image " << anchor_image << std::endl;
    use_NN_interpolator = use_NN;
    global_anchor = anchor_image;
    global_image_path = image_path;
    global_joint_register = new fregl_joint_register(xml_filename);
    //    std::cout << "Creating Space Transformer" << std::endl;
    global_space_transformer = new fregl_space_transformer(global_joint_register);
    //    std::cout << "Created Space Transformer" << std::endl;
    std::cout << "Setting Anchor" << std::endl;
    global_space_transformer->set_anchor(global_anchor, false, false);
    //    std::cout << "Anchor Set" << std::endl;
    image_names = global_space_transformer->image_names();
    global_origin = global_space_transformer->origin();
    global_size = global_space_transformer->montage_size();
    roi_origin = global_origin;
    roi_size = global_size;
    std::cout << "GLOBAL Origin = " << global_origin[0] << "," << global_origin[1] << "," << global_origin[2] << std::endl;
    std::cout << "GLOBAL Size = " << global_size[0] << " x " << global_size[1] << " x " << global_size[2] << std::endl;
    global_space_transformer->set_roi(roi_origin, roi_size);
    global_channel = 0;
    global_use_channel = false;

}

//: Set the Region of Interest
//
//  Index and size are updated accordingly.

void fregl_image_manager::set_regionofinterest(PointType origin, SizeType size) {
    roi_origin = origin;
    roi_size = size;
    if (roi_origin[0] < global_origin[0]) roi_origin[0] = global_origin[0];
    if (roi_origin[1] < global_origin[1]) roi_origin[1] = global_origin[1];
    if (roi_origin[2] < global_origin[2]) roi_origin[2] = global_origin[2];
    if (roi_origin[0] > global_origin[0] + global_size[0]) roi_origin[0] = global_origin[0] + global_size[0];
    if (roi_origin[1] > global_origin[1] + global_size[1]) roi_origin[1] = global_origin[1] + global_size[1];
    if (roi_origin[2] > global_origin[2] + global_size[2]) roi_origin[2] = global_origin[2] + global_size[2];
    // Make sure size is in range or fix
    if (roi_origin[0] + roi_size[0] > global_origin[0] + global_size[0])
        roi_size[0] = (global_origin[0] + global_size[0]) - roi_origin[0];
    if (roi_origin[1] + roi_size[1] > global_origin[1] + global_size[1])
        roi_size[1] = (global_origin[1] + global_size[1]) - roi_origin[1];
    if (roi_origin[2] + roi_size[2] > global_origin[2] + global_size[2])
        roi_size[2] = (global_origin[2] + global_size[2]) - roi_origin[2];
    global_space_transformer->set_roi(roi_origin, roi_size);
}

//: Set the Region of Interest
//
//  Index and size are updated accordingly.

void fregl_image_manager::set_regionofinterest(IndexType origin, SizeType size) {
    roi_origin[0] = origin[0];
    roi_origin[1] = origin[1];
    roi_origin[2] = origin[2];
    roi_size = size;
    set_regionofinterest(roi_origin, roi_size);
}
//:  Build the Montage based on the current Region
//   "s_origin" is the new roi origin point
//   "s_size" is the new roi size

void fregl_image_manager::Update() {

    std::string image_name = global_image_path + std::string("/") + image_names[0];
    ImageType::Pointer image, xformed_image;
    image = fregl_util_read_image(image_name, global_use_channel, global_channel, false);
    std::cout << "Composing the final image ..." << std::endl;
    montage_image = global_space_transformer->transform_image_roi(image, 0, 0, use_NN_interpolator);
    for (unsigned int i = 1; i < image_names.size(); i++) {
        if (global_space_transformer->image_in_roi(i)) {
            image_name = global_image_path + std::string("/") + image_names[i];
            image = fregl_util_read_image(image_name, global_use_channel, global_channel, false);
            xformed_image = global_space_transformer->transform_image_roi(image, i, 0, use_NN_interpolator);
        } else
            continue;
        if (!xformed_image)
            continue;

        //fuse the image
        RegionConstIterator inputIt(xformed_image, xformed_image->GetRequestedRegion());
        RegionIterator outputIt(montage_image, montage_image->GetRequestedRegion());

        for (inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd();
                ++inputIt, ++outputIt) {
            outputIt.Set(vnl_math_max(outputIt.Get(), inputIt.Get()));
        }
    }
}

//: Return an ITK image pointer to the current Montage
//

fregl_image_manager::ImageType::Pointer fregl_image_manager::GetOutput() {
    typedef itk::ImageDuplicator< ImageType > DuplicatorType;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(montage_image);
    duplicator->Update();
    ImageType::Pointer image = duplicator->GetOutput();
    	
	typedef itk::TranslationTransform<double, ImageType::ImageDimension> TranslationTransformType;
	TranslationTransformType::Pointer transform = TranslationTransformType::New();
	TranslationTransformType::OutputVectorType translation;
	std::cout << "global_origin: " << global_origin[0] << " " << global_origin[1] << " " << global_origin[2] << std::endl;
	translation[0] = global_origin[0];
	translation[1] = global_origin[1];
	translation[2] = global_origin[2];
	transform->Translate(translation);

	typedef itk::ResampleImageFilter<ImageType, ImageType, double> ResampleImageFilterType;
	ResampleImageFilterType::Pointer resampleFilter = ResampleImageFilterType::New();
	resampleFilter->SetTransform(transform.GetPointer());
	resampleFilter->SetInput(image);
	
	ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
	resampleFilter->SetSize(size);
	
	try
	{
		resampleFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cout << "Error in ResampleFilter" << std::endl;
		std::cout << err << std::endl;
	}

	
	return resampleFilter->GetOutput();
}

//: Set the image channel to montage
//

void fregl_image_manager::set_channel(int channel) {
    global_channel = channel;
    global_use_channel = true;
}

//: Return the global space origin
//

fregl_image_manager::PointType fregl_image_manager::get_global_origin() {
    return global_origin;
}

//: Return return the global space size
//

fregl_image_manager::SizeType fregl_image_manager::get_global_size() {
    return global_size;
}

//: Return the current Region origin
//

fregl_image_manager::PointType fregl_image_manager::get_region_origin() {
    return roi_origin;
}

//: Return the current Region Size
//

fregl_image_manager::SizeType fregl_image_manager::get_region_size() {
    return roi_size;
}

//: Return the anchor image name

std::string const
fregl_image_manager::get_anchor_name() {
    return global_anchor;
}

//: Return a pointer the space transformer

fregl_space_transformer::Pointer fregl_image_manager::get_space_transformer() {
    return global_space_transformer;
}