#ifndef MONTAGE_REGION_SELECTION
#define MONTAGE_REGION_SELECTION

#include "ftkImage/ftkImage.h"
#include "vtkTable.h"

/*Utility class for storing the objects related to the current ROI selected
  and all the edits related to it from the Nucleus Editor
*/
class MontageRegionSelection
{
  public:
  MontageRegionSelection();
  ~MontageRegionSelection();

  private:
  ftk::Image::Pointer ChannelImage;
  ftk::Image::Pointer LabelImage;
  vtkSmartPointer<vtkTable> Table;

  public:
  void SetChannelImage(ftk::Image::Pointer Image) { ChannelImage = Image; }
  void SetLabelImage(ftk::Image::Pointer Image)   { LabelImage = Image; }
  void SetTable(vtkSmartPointer<vtkTable> table)  { Table = table; }

  ftk::Image::Pointer GetChannelImage(void) { return ChannelImage; }
  ftk::Image::Pointer GetLabelImage(void)   { return LabelImage; }
  vtkSmartPointer<vtkTable> GetTable(void)  { return Table; }

  //Edits returned from the nucleus editor
};

#endif //MONTAGE_REGION_SELECTION
