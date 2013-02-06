template<typename pixelType> MontageView::NucleusEditorLabelType::Pointer
	MontageView::RelabelImage(ftk::Image::Pointer InputIm)
{
  typedef typename itk::Image<pixelType, 3> InputLabelType;
  typedef typename itk::CastImageFilter< InputLabelType, MontageView::NucleusEditorLabelType > CastFilterType;

  typename InputLabelType::Pointer InputImage = InputIm->GetItkPtr<pixelType>(0,0);
  pixelType* buf = (pixelType*)InputImage->GetBufferPointer();

  const ftk::Image::Info * imInfo = InputIm->GetImageInfo();
  itk::SizeValueType BufferSize = imInfo->numColumns*imInfo->numRows*imInfo->numZSlices;

#ifdef _OPENMP
  #pragma omp parallel for
#if _OPENMP < 200805L
  for( itk::IndexValueType i=0; i<BufferSize; ++i )
#else
  for(  itk::SizeValueType i=0; i<BufferSize; ++i )
#endif
#else
  for(  itk::SizeValueType i=0; i<BufferSize; ++i )
#endif
  {
    if( buf[i] )
    {
      std::map< itk::SizeValueType , itk::SizeValueType >::iterator it;
      it = LabelToRelabelMap.find( buf[i] );
      if( it!=LabelToRelabelMap.end() ) buf[i] = it->second;
      else
      {
#pragma omp critical
	buf[i] = this->InsertNewLabelToRelabelMap( buf[i] );
      }
    }
  }

  typename CastFilterType::Pointer caster = CastFilterType::New();
  caster->SetInput( InputImage );
  try
  {
    caster->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Error is casting label for display:"
		<< excep << std::endl;
    return NULL;
  }

  this->NumberOfLabelsFound = LabelToRelabelMap.size();

  return caster->GetOutput();
}
