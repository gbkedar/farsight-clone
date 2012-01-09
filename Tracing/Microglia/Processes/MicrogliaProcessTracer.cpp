#include "MicrogliaProcessTracer.h"
  
bool CompareNeighbors(std::pair< double, Node *> n1,
                      std::pair< double, Node *> n2)
{
  if(n1.first < n2.first)
    {
    return true;
    }
  return false;
}

Node::Node()
  {
    this->parent = NULL;
    this->IsOpen = true;
    this->type = 3;
  }

MicrogliaProcessTracer::MicrogliaProcessTracer()
{
  this->NodeCounter = 1;
  this->Padding = 1;
  this->ProcessRadius = 0.3;
  this->SeparateFilePerCell = false;

  //no default distance threshold
  //this means all nodes will be added to a tree
  this->MaxDistance = std::numeric_limits<double>::max();
}

MicrogliaProcessTracer::~MicrogliaProcessTracer()
{
  std::vector< Node * >::iterator nodeItr;
  for(nodeItr = this->Open.begin(); nodeItr != this->Open.end(); ++nodeItr)
    {
    delete (*nodeItr);
    }
  for(nodeItr = this->Closed.begin(); nodeItr != this->Closed.end(); ++nodeItr)
    {
    delete (*nodeItr);
    }
} 

void MicrogliaProcessTracer::LoadInputImage(std::string fname)
{
  std::cout << "Reading input file "<< fname << std::endl;
  ReaderType::GlobalWarningDisplayOff();
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fname);
  ImageType3D::Pointer image = reader->GetOutput();
  image->Update();

  this->LoadInputImage(image);
}

void MicrogliaProcessTracer::LoadInputImage(ImageType3D::Pointer &image)
{
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetOutputMinimum(0.0);
  rescaler->SetOutputMaximum(1.0);
  rescaler->SetInput(image);
  this->InputImage = rescaler->GetOutput();
  rescaler->Update();

  //pad z slices
  itk::Size<3> isz = this->InputImage->GetBufferedRegion().GetSize();
  itk::Size<3> osz = isz;
  osz[2] += 2*this->Padding;
  itk::Index<3> indx, ondx;
  
  this->PaddedInputImage = ImageType3D::New();
  this->PaddedInputImage->SetRegions(osz);
  this->PaddedInputImage->Allocate();
  this->PaddedInputImage->SetSpacing(this->InputImage->GetSpacing());
  
  for(ondx[2] = 0; ondx[2] < (int)osz[2]; ++ondx[2]) 
  {
    indx[2] = (ondx[2] < (int)this->Padding) ? 0 : ondx[2] - this->Padding;
    indx[2] = (ondx[2] >= (int)(osz[2]-this->Padding)) ? isz[2]-1 : indx[2];
    for(ondx[1] = 0; ondx[1] < (int)osz[1]; ++ondx[1]) 
    {
      indx[1] = ondx[1];
      for(ondx[0] = 0; ondx[0] < (int)osz[0]; ++ondx[0]) 
      {
        indx[0] = ondx[0];
        this->PaddedInputImage->SetPixel(ondx, this->InputImage->GetPixel(indx));
      }
    }
  }

  std::cout << "Input file size (after zero padding) is " << this->PaddedInputImage->GetBufferedRegion().GetSize() << std::endl;
}

///////////////////////////////////////////////////////////////////////
std::vector< Node * > MicrogliaProcessTracer::ReadListOfPoints(std::string fname)
{
  std::vector< Node * > listOfPoints;
  std::string temp, num;
  std::ifstream infile;
  infile.open(fname.c_str());
  size_t x1, x2;

  while(!infile.eof()) 
  {
    std::getline(infile,temp);
    if (temp.length() < 1)
      continue;
    
    x1 = temp.find_first_of("0123456789.");
    x2 = temp.find_first_not_of("0123456789.",x1);
    if ((x2 - x1) > 10)
      continue;
    
    num = temp.substr(x1,x2-x1);
    int x = atoi(num.c_str());

    x1 = temp.find_first_of("0123456789.",x2+1);
    x2 = temp.find_first_not_of("0123456789.",x1);
    if ((x2 - x1) > 10)
      continue;
    
    num = temp.substr(x1,x2-x1);
    int y = atoi(num.c_str());

    x1 = temp.find_first_of("0123456789.",x2+1);
    x2 = temp.find_first_not_of("0123456789.",x1);
    if (x2 > temp.length())
      x2 = temp.length();
    
    if ((x2 - x1) > 10)
      continue;
    
    num = temp.substr(x1,x2-x1);
    int z = atoi(num.c_str());

    itk::Index<3> idx;
    idx[0] = x;
    idx[1] = y;
    idx[2] = z;
    Node *n = new Node();
    n->ID = this->NodeCounter++;
    n->index = idx;
    listOfPoints.push_back(n);
  }
  infile.close();

  return listOfPoints;
}
  
////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::LoadSeedPoints(std::string fname)
{
  this->Closed = this->ReadListOfPoints(fname);

  //nodes are created open by default, so close all these ones.
  std::vector< Node * >::iterator nodeItr;
  for(nodeItr = this->Closed.begin(); nodeItr != this->Closed.end(); ++nodeItr)
    {
    (*nodeItr)->IsOpen = false;
    (*nodeItr)->type = 1;
    }

  this->Roots = this->Closed;
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::LoadNodes(std::string fname)
{
  this->Open = this->ReadListOfPoints(fname);
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::LoadSomaImage(std::string somaFileName)
{
  typedef itk::ImageFileReader<CharImageType3D> SomaReaderType;
  SomaReaderType::Pointer somaReader = SomaReaderType::New();
  somaReader->SetFileName(somaFileName);
  this->SomaImage = somaReader->GetOutput();
  somaReader->Update();
  this->SomaImage->SetSpacing(this->InputImage->GetSpacing());
}

///////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::RunTracing()
{
  if(!this->InputImage)
    {
    std::cerr << "You must call LoadImage() before RunTracing()" << std::endl;
    return;
    }
  if(this->Closed.empty())
    {
    std::cerr << "You must call LoadSeedPoints() before RunTracing()" << std::endl;
    return;
    }

  this->CalculateCriticalPoints();

  itk::ImageRegionConstIterator<ImageType3D> Nit(this->NDXImage, this->NDXImage->GetBufferedRegion());
  for (Nit.GoToBegin(); !Nit.IsAtEnd(); ++Nit) 
    {
    if (Nit.Get() > 0) 
      {
      Node *n = new Node();
      n->ID = this->NodeCounter++;
      n->index = Nit.GetIndex();
      this->Open.push_back(n);
      this->IndexToNodeMap[n->index] = n;
      }
    }
 
  //add the centroids to the Index -> Node map too
  std::vector< Node * >::iterator nodeItr;
  for(nodeItr = this->Closed.begin(); nodeItr != this->Closed.end(); ++nodeItr)
    {
    Node *n = (*nodeItr);
    this->IndexToNodeMap[n->index] = n;
    }
  
  std::cout << "open list size: " << this->Open.size() << std::endl;
  std::cout << "closed list size: " << this->Closed.size() << std::endl;

  //use RATS to threshold the input image.
  //this is used in GetDistanceBetweenPoints 
  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< ImageType3D, ImageType3D > GradientType;
  GradientType::Pointer gradient = GradientType::New();
  gradient->SetInput( this->InputImage );
  gradient->SetSigma( 0.25 );
  gradient->Update();
  typedef itk::RobustAutomaticThresholdImageFilter< ImageType3D, ImageType3D > ThresholdFilterType;
  ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
  thresholdFilter->SetInput( this->InputImage );
  thresholdFilter->SetGradientImage( gradient->GetOutput() );
  thresholdFilter->SetPow( 1.0 );
  this->ThresholdedImage = thresholdFilter->GetOutput();
  thresholdFilter->Update();
 
  std::cout << "computing adjacency matrix for seed points" << std::endl;
  this->ComputeAdjacencies(this->Closed);
  
  std::cout << "computing adjacency matrix for critical points" << std::endl;
  this->ComputeAdjacencies(this->Open);
  
  this->BuildTrees();

  std::cout << "cleaning up nodes that fall within a soma" << std::endl;
  this->PruneSomaNodes();
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::CalculateCriticalPoints(void)
{
  std::cout << std::endl<< "Searching for critical points" << std::endl;
  this->NDXImage = ImageType3D::New();
  this->NDXImage->SetRegions(this->PaddedInputImage->GetBufferedRegion());
  this->NDXImage->Allocate();
  this->NDXImage->FillBuffer(0.0f);
  
  float power = -0.25;

  for (unsigned int i = 0; i < 6; ++i)
  {
    power = power + 0.25;
    float sigma = vcl_pow(2, power);
    std::cout << "Analysis at sigma " << sigma << std::endl;
    CalculateCriticalPointsAtScale( sigma );
  }
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::CalculateCriticalPointsAtScale( float sigma ) 
{
  typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType3D , ImageType3D> GFilterType;
  GFilterType::Pointer gauss = GFilterType::New();
  gauss->SetInput( this->PaddedInputImage );
  gauss->SetSigma( sigma );
  gauss->SetNormalizeAcrossScale(false);
  gauss->GetOutput()->Update();
  
  itk::ImageRegionIterator<ImageType3D> ittemp(gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());
  float tnorm = vcl_pow(sigma, 1.6f);
  for(ittemp.GoToBegin(); !ittemp.IsAtEnd(); ++ittemp) 
  {
    float q = ittemp.Get()*tnorm;
    ittemp.Set(-1.0f*q);
  }

  // set the diagonal terms in neighborhood iterator
  itk::Offset<3>
    xp =  {{2 ,  0 ,   0}},
    xn =  {{-2,  0,    0}},
    yp =  {{0,   2,   0}},
    yn =  {{0,  -2,    0}},
    zp =  {{0,   0,    2}},
    zn =  {{0,   0,   -2}};

  itk::Size<3> rad = {{1,1,1}};
  itk::NeighborhoodIterator<ImageType3D> nit(rad , gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());
  itk::ImageRegionIterator<ImageType3D> it(gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());

  unsigned int
    xy1 =  17, //{ 1 ,   1 ,  0 },
    xy2 =  9,  //{ -1,  -1 ,  0 },
    xy3 =  15, //{ -1,   1 ,  0 },
    xy4 =  11, //{ 1 ,  -1 ,  0 },

    yz1 =  25, //{ 0 ,   1 ,  1 },
    yz2 =  1,  //{ 0 ,  -1 , -1 },
    yz3 =  19, //{ 0 ,  -1 ,  1 },
    yz4 =  7,  //{ 0 ,   1 , -1 },

    xz1 =  23, //{ 1 ,   0 ,  1 },
    xz2 =  3,  //{-1 ,   0 , -1 },
    xz3 =  21, //{-1 ,   0 ,  1 },
    xz4 =  5;  //{ 1 ,   0 , -1 };

  typedef itk::FixedArray< double, 3 > EigenValuesArrayType;
  typedef itk::Matrix< double, 3, 3 > EigenVectorMatrixType;
  typedef itk::SymmetricSecondRankTensor<double,3> TensorType;

  itk::Size<3> sz = this->PaddedInputImage->GetBufferedRegion().GetSize();
  sz[0] = sz[0] - 3;
  sz[1] = sz[1] - 3; 
  sz[2] = sz[2] - 3;

  it.GoToBegin();
  nit.GoToBegin();
  itk::Vector<float,3> sp = this->PaddedInputImage->GetSpacing();

  long ctCnt = 0;
  ImageType3D::PointType origin = this->InputImage->GetOrigin();
  while(!nit.IsAtEnd()) 
  {
    itk::Index<3> ndx = it.GetIndex();
    if ( (ndx[0] < 2) || (ndx[1] < 2) || (ndx[2] < 2) ||
      (ndx[0] > (int)sz[0]) || (ndx[1] > (int)sz[1]) ||
      (ndx[2] > (int)sz[2]) )
    {
      ++it;
      ++nit;
      continue;
    }

    float a1 = 0.0;
    for (unsigned int i=0; i < 13; ++i)
    {
      a1 += vnl_math_max(nit.GetPixel(i), nit.GetPixel(26 - i));
    }
    
    float val = nit.GetPixel(13) ;

    const float thresh1 = 0.03;   // 3% of maximum theshold from Lowe 2004
    const float thresh2 = 0.001;  // -0.1 percent of range

    if ( ((val - a1/13.0f) > thresh2 ) && ( val > thresh1 ))  
      {
      TensorType hessian;
      hessian[0] = gauss->GetOutput()->GetPixel( ndx + xp ) +
                   gauss->GetOutput()->GetPixel( ndx + xn ) -
                   2*nit.GetPixel( 13 );
      hessian[3] = gauss->GetOutput()->GetPixel( ndx + yp ) +
                   gauss->GetOutput()->GetPixel( ndx + yn ) -
                   2*nit.GetPixel( 13 );
      hessian[5] = gauss->GetOutput()->GetPixel( ndx + zp ) +
                   gauss->GetOutput()->GetPixel( ndx + zn ) -
                   2*nit.GetPixel( 13 );
      hessian[1] = nit.GetPixel(xy1) + nit.GetPixel(xy2) -
                   nit.GetPixel(xy3) - nit.GetPixel(xy4);
      hessian[2] = nit.GetPixel(xz1) + nit.GetPixel(xz2) -
                   nit.GetPixel(xz3) - nit.GetPixel(xz4);
      hessian[4] = nit.GetPixel(yz1) + nit.GetPixel(yz2) -
                   nit.GetPixel(yz3) - nit.GetPixel(yz4);
      EigenValuesArrayType ev;
      EigenVectorMatrixType em;
      hessian.ComputeEigenAnalysis (ev, em);
      unsigned int w = GetEigenvalueL1(ev);

      float value = 0.0;
      switch(w)
        {
        case 0:
          value = -1 * (ev[1] + ev[2]);
          break;
        case 1:
          value = -1 * (ev[0] + ev[2]);
          break;
        case 2:
          value = -1 * (ev[0] + ev[1]);
          break;
        default:
          std::cout << "impossible switch value" << std::endl;
          break;
        }
      value -= ev[w];

      if (RegisterIndex(value, ndx, sz)) 
        {
        this->NDXImage->SetPixel(ndx,value);
        ctCnt++;
        }
      }
    ++it;
    ++nit;
  }
  std::cout <<"Number of CTs at this stage: " << ctCnt <<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
unsigned int MicrogliaProcessTracer::
  GetEigenvalueL1(const itk::FixedArray<float, 3> &ev)
{
  unsigned int w = 0;
  float ev1 = abs(ev[0]);
  float ev2 = abs(ev[1]);
  float ev3 = abs(ev[2]);

  if ( (ev2 < ev1) && (ev2 < ev3) )
    {
    w = 1;
    }
  
  if ( (ev3 < ev1) && (ev3 < ev2) )
    {
    w = 2;
    }

  return w;
}

////////////////////////////////////////////////////////////////////////////////
bool MicrogliaProcessTracer::RegisterIndex(const float value,
                                           itk::Index<3> &ndx,
                                           itk::Size<3>& sz) 
{
  itk::Index<3> n;
  bool higherPresent = false;

  unsigned int radX = (unsigned int)
    ceil( (this->ProcessRadius / this->InputImage->GetSpacing()[0] ) );
  unsigned int radY = (unsigned int)
    ceil( (this->ProcessRadius / this->InputImage->GetSpacing()[1] ) );
  unsigned int radZ = (unsigned int)
    ceil( (this->ProcessRadius / this->InputImage->GetSpacing()[2] ) );

  for (n[0] = ndx[0]-radX; n[0] <= (int)(ndx[0]+radX); ++n[0]) 
  {
    for (n[1] = ndx[1]-radY; n[1] <= (int)(ndx[1]+radY); ++n[1]) 
    {
      for (n[2] = ndx[2]-radZ; n[2] <= (int)(ndx[2]+radZ); ++n[2]) 
      {
        if ( (n[0] < 2) || (n[1] < 2) || (n[2] < 2) || (n[0] > (int)sz[0]) ||
          (n[1] > (int)sz[1]) || (n[2] > (int)sz[2]) )
        {
          continue;
        }

        float curval = this->NDXImage->GetPixel(n);
        if (value > curval) 
        {
          this->NDXImage->SetPixel(n,0.0f);
        }
        else if (value < curval) 
        {
          higherPresent = true;       
        }
      }
    }
  }
  if (higherPresent == true) 
  {
    return false;
  }
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////
float MicrogliaProcessTracer::GetRadius(itk::Index<3> idx) 
{
  float r = 2.0f;
  itk::Vector<float,3> m1, m2, m, pos;
  itk::Index<3> ndx;

  pos[0] = idx[0];
  pos[1] = idx[1];
  pos[2] = idx[2];
      
  for (int iter = 0; iter < 20; ++iter) 
  {
    for( int i = 0; i<3; i++) 
    {
      m1[i] = pos[i] - vnl_math_max(2.0f*r, 5.0f);
      m2[i] = pos[i] + vnl_math_max(2.0f*r, 5.0f);
    }

    std::vector<float> c;
    c.reserve(4*4*int(r*r));
    float i1 = 0.0f, i2 = 0.0f, i1s = 0.0f, i2s = 0.0f;
    for (m[2] = m1[2]; m[2] <= m2[2]; m[2]++) 
    {
      for (m[1] = m1[1]; m[1] <= m2[1]; m[1]++) 
      {
        for (m[0] = m1[0]; m[0] <= m2[0]; m[0]++) 
        {
          ndx.CopyWithRound(m);
          itk::Vector<float,3> mm = pos - m;
          float d = mm.GetNorm();
          if (this->PaddedInputImage->GetBufferedRegion().IsInside(ndx)) 
          {
            float val = this->PaddedInputImage->GetPixel(ndx);
            if (d < r) 
            {
              i1 += val;
              ++i1s;
            }
            else 
            {
              i2 += val;
              ++i2s;
            }

            if (vnl_math_abs(d - r) < 0.7f) 
              c.push_back(val);           
          }
        }
      }
    }
    i1 /= i1s;
    i2 /= i2s;
    float dr = 0.0f;
    for (std::vector<float>::iterator it = c.begin(); it < c.end(); ++it) 
      dr += vnl_math_abs((*it) - i1) - vnl_math_abs((*it) - i2);
    
    dr *= 1.0f / float(c.size()); //rate
    dr = vnl_math_max(dr , -1.0f);
    dr = vnl_math_min(dr , 1.0f);
    r -= dr;
    r = vnl_math_max(r , 1.0f) ; 
  }
  return r;
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::BuildTrees()
{
  std::pair< Node *, Node * > parentAndChild;
  Node *parent = NULL;
  Node *child = NULL;

  while( !( this->Open.empty() ) )
    {
    parentAndChild = this->FindClosestOpenNode();
    if(parentAndChild.first == NULL)
      {
      //We've added all the nodes within our distance threshold.
      std::cout << this->Open.size() << " nodes were too far away to get added to a cell model" << std::endl;
      break;
      }
    
    parent = parentAndChild.first;
    child = parentAndChild.second;

    //setup the parent/child relationship between these two nodes
    parent->children.push_back(child);
    child->parent = parent;
    
    //remove the child node from the Open list
    bool itsThere = false;
    if(std::find(this->Open.begin(), this->Open.end(), child) != this->Open.end())
      {
      itsThere = true;
      }
    if(!itsThere)
      {
      std::cout << "(Debug) PROBLEM: Before deletion, child is NOT in Open." << std::endl;
      }
    std::vector< Node* >::iterator newEnd =
      std::remove(this->Open.begin(), this->Open.end(), child);
    this->Open.erase(newEnd, this->Open.end());
    if(std::find(this->Open.begin(), this->Open.end(), child) != this->Open.end())
      {
      itsThere = true;
      }
    else
      {
      itsThere = false;
      }
    if(itsThere)
      {
      std::cout << "(Debug) PROBLEM: After deletion, child is still in Open." << std::endl;
      }
  
    //and add it to the Closed list
    child->IsOpen = false;
    this->Closed.push_back(child);
    }
}

////////////////////////////////////////////////////////////////////////////////
// This function returns (NULL, NULL) when the closest open node is further away
// than our distance threshold
std::pair< Node *, Node * > MicrogliaProcessTracer::FindClosestOpenNode()
{
  std::pair< Node *, Node * > parentAndChild;
  parentAndChild.first = NULL;
  parentAndChild.second = NULL;

  std::pair< double, Node * > distanceAndNeighbor;
  std::list< std::pair< double, Node * > > *neighborList;
  double minDistance = std::numeric_limits<double>::max();
  std::vector< Node * >::iterator nodeItr, nbrItr;
 
  for(nodeItr = this->Closed.begin(); nodeItr != this->Closed.end(); ++nodeItr)
    {
    neighborList = &(this->AdjacencyMap[(*nodeItr)]);
    bool openNeighborFound = false;
    while(!openNeighborFound && !neighborList->empty())
      {
      distanceAndNeighbor = neighborList->front();
      if(!distanceAndNeighbor.second->IsOpen)
        {
        neighborList->pop_front();
        }
      else
        {
        openNeighborFound = true;
        if(distanceAndNeighbor.first < minDistance)
          {
          minDistance = distanceAndNeighbor.first;
          parentAndChild.first = (*nodeItr);
          parentAndChild.second = distanceAndNeighbor.second;
          }
        }
      }
    }

 //the child node is closed in the caller function
 return parentAndChild;
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::WriteToSWC( std::string fname )
{
  std::cout << "writing output to disk" << std::endl;

  if(this->SeparateFilePerCell)
    {
    this->WriteMultipleSWCFiles(fname);
    }
  else
    {
    this->WriteSingleSWCFile(fname);
    }
}
  
////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::WriteSingleSWCFile( std::string fname )
{
  std::ofstream output( fname.c_str() );
  
  std::vector< Node * >::iterator nodeItr;
  
  ImageType3D::PointType origin = this->InputImage->GetOrigin();
  for(nodeItr = this->Closed.begin(); nodeItr != this->Closed.end(); ++nodeItr)
    {
    itk::Index<3> idx = (*nodeItr)->index;
    float x = idx[0] + origin[0];
    float y = idx[1] + origin[1];
    float z = idx[2] + origin[2];
    long parentID;

    //compute vessel radius at this point
    float r = this->GetRadius( idx );
  
    if ( (*nodeItr)->parent == NULL )
      {
      parentID = -1;
      }
    else
      {
      parentID = (*nodeItr)->parent->ID;
      }
    
    output << (*nodeItr)->ID << " " << (*nodeItr)->type << " " << x << " " << y << " " << z << " " << r << " " << parentID << std::endl; 
    }
  output.close();
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::WriteMultipleSWCFiles( std::string fname )
{
  std::vector< Node * >::iterator rootItr, childItr;
  int n = 1;
  //print out each cell to a separate file
  for(rootItr = this->Roots.begin(); rootItr != this->Roots.end(); ++rootItr)
    {
    //open a new file for this cellular model
    size_t pos = fname.find(".swc");
    std::string baseName = fname.substr(0,pos);
    std::stringstream ss;
    ss << baseName << n << ".swc";
    std::string filename =  ss.str();
    std::ofstream outFile( filename.c_str() );
    std::cout << "writing out " << filename << std::endl;

    Node *root = (*rootItr);
    //this function recursively writes out all the descendents of root to
    //the same file.
    this->WriteNodeToSWCFile(root, &outFile); 
    outFile.close();
    ++n;
    }
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::WriteNodeToSWCFile( Node *n, std::ofstream *outFile)
{
  if(n->IsOpen || n->ID < 1 || n->ID > this->NodeCounter)
    {
    return;
    }
  ImageType3D::PointType origin = this->InputImage->GetOrigin();
  itk::Index<3> idx = n->index;
  float x = idx[0] + origin[0];
  float y = idx[1] + origin[1];
  float z = idx[2] + origin[2];
  long parentID;

  //compute vessel radius at this point
  float r = this->GetRadius( idx );

  if ( n->parent == NULL )
    {
    parentID = -1;
    }
  else
    {
    parentID = n->parent->ID;
    }
  
  *outFile << n->ID << " " << n->type << " " << x << " " << y << " " << z <<
             " " << r << " " << parentID << std::endl; 
    
  std::vector< Node * >::iterator childItr;
  for(childItr = n->children.begin(); childItr != n->children.end();
      ++childItr)
    {
    Node *child = (*childItr);
    if(child->parent->ID != n->ID)
      {
      continue;
      }
    this->WriteNodeToSWCFile(child, outFile);
    }
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::ComputeAdjacencies( std::vector< Node * > nodes )
{
  std::vector< Node * >::iterator nodeItr;
  Node *nbrNode;
  
  ImageType3D::SizeType radius;
  radius[0] = (unsigned int)
    ceil( (this->MaxDistance / this->InputImage->GetSpacing()[0] ) ) * 2;
  radius[1] = (unsigned int)
    ceil( (this->MaxDistance / this->InputImage->GetSpacing()[1] ) ) * 2;
  radius[2] = (unsigned int)
    ceil( (this->MaxDistance / this->InputImage->GetSpacing()[2] ) ) * 2;
 
  itk::ConstNeighborhoodIterator<ImageType3D>
    nbrItr(radius, this->NDXImage, this->NDXImage->GetBufferedRegion());

  for(nodeItr = nodes.begin(); nodeItr != nodes.end(); ++nodeItr)
    {
    itk::Index<3> start = (*nodeItr)->index;
    itk::Point<double,3> startPoint;
    this->InputImage->TransformIndexToPhysicalPoint( start, startPoint );
    std::list< std::pair< double, Node *> > neighbors;

    nbrItr.SetLocation(start);
    itk::Size<3> nbrhdSize = nbrItr.GetSize();
    unsigned int totalNbrhdSize = nbrhdSize[0] * nbrhdSize[1] * nbrhdSize[2]; 
    
    for (unsigned int i=0; i < totalNbrhdSize; ++i)
      {
      if(nbrItr.GetPixel(i) == 0.0)
        {
        continue;
        }

      itk::Index<3> end = nbrItr.GetIndex(i);
      if(start == end)
        {
        continue;
        }

      itk::Point<double,3> endPoint;
      this->InputImage->TransformIndexToPhysicalPoint( end, endPoint );

      double euclideanDistance = startPoint.EuclideanDistanceTo(endPoint);
      double scaledDistance = this->GetDistanceBetweenPoints(start, end);
     
      if( scaledDistance > this->MaxDistance )
        {
        continue;
        }

      std::map<itk::Index<3>, Node *, CompareIndices>::iterator mapItr =
        this->IndexToNodeMap.find(end);
      if(mapItr == this->IndexToNodeMap.end())
        {
        continue;
        }

      nbrNode = (*mapItr).second; 
      if( (*nodeItr)->ID == nbrNode->ID )
        {
        std::cout << "this should never happen" << std::endl;
        continue;
        }

      std::pair< double, Node * > nbrPair;
      nbrPair.first = euclideanDistance;
      nbrPair.second = nbrNode;
      neighbors.push_front(nbrPair);
      }

    neighbors.sort(CompareNeighbors);
    this->AdjacencyMap[(*nodeItr)] = neighbors;
    }
}

////////////////////////////////////////////////////////////////////////////////
double MicrogliaProcessTracer::GetDistanceBetweenPoints(itk::Index<3> start,
                                                        itk::Index<3> end)
{
  ImageType3D::SizeType imageSize =
    this->ThresholdedImage->GetLargestPossibleRegion().GetSize();

  //compute Euclidean distance between the two points
  itk::Point<double,3> startPoint;
  itk::Point<double,3> endPoint;
  this->InputImage->TransformIndexToPhysicalPoint( start, startPoint );
  this->InputImage->TransformIndexToPhysicalPoint( end, endPoint );
  double euclideanDistance = startPoint.EuclideanDistanceTo(endPoint);
 
  //analyze the pixels along the line between the start & end points.
  //the distance returned by this function is typically less than the actual
  //Euclidean distance.
  //How much less is based on how many "foreground" pixels lie along this line.
  std::vector< itk::Index<3> > indices = this->Line.BuildLine(start, end);
  float scaledDistance = 0;
  for(unsigned int i = 0; i < indices.size(); i++)
    {
    float distanceWeight = 0.5;
    if(indices[i] == start)
      {
      continue;
      }
    if(indices[i] == end)
      {
      break;
      }

    //search for foreground pixels in each 3x3x3 neighborhood
    bool foregroundPixelFound = false;
    for(int x = indices[i][0] - 1; x < indices[i][0] + 2; x++)
      {
      if(x < 0 || x > (int)imageSize[0] - 1) 
        {
        continue;
        }
      for(int y = indices[i][1] - 1; y < indices[i][1] + 2; y++)
        {
        if(y < 0 || y > (int)imageSize[1] - 1) 
          {
          continue;
          }
        for(int z = indices[i][2] - 1; z < indices[i][2] + 2; z++)
          {
          if(z < 0 || z > (int)imageSize[2] - 1)
            {
            continue;
            }
          //here's where we check for foreground pixels
          itk::Index<3> nbrIdx = {x, y, z};
          if(this->ThresholdedImage->GetPixel(nbrIdx) != 0)
            {
            foregroundPixelFound = true;
            break;
            }
          }
        if(foregroundPixelFound)
          {
          break;
          }
        }
      if(foregroundPixelFound)
        {
        break;
        }
      }
    if(!foregroundPixelFound)
      {
      distanceWeight = 1.0;
      }

    //take image spacing into account
    itk::Offset<3> diff = indices[i] - indices[i-1];
    if(diff[2] != 0)
      {
      scaledDistance += ( this->InputImage->GetSpacing()[2] * distanceWeight ); 
      }
    else if(diff[1] != 0)
      {
      scaledDistance += ( this->InputImage->GetSpacing()[1] * distanceWeight ); 
      }
    else
      {
      scaledDistance += ( this->InputImage->GetSpacing()[0] * distanceWeight ); 
      }
    }

  return scaledDistance;
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::PruneSomaNodes()
{

  //first thing: need to treat the soma image as a shape label object map
  typedef itk::ShapeLabelObject< unsigned long, 3 > LabelObjectType;
  typedef itk::LabelMap< LabelObjectType > LabelMapType;

  typedef itk::BinaryImageToShapeLabelMapFilter<CharImageType3D, LabelMapType>
    BinaryToLabelType;
  BinaryToLabelType::Pointer binaryToLabelConverter = BinaryToLabelType::New();
  binaryToLabelConverter->SetInput(this->SomaImage);
  LabelMapType::Pointer labelMap = binaryToLabelConverter->GetOutput();
  binaryToLabelConverter->Update();

  std::vector<Node *> somaBranches;

  // Loop over each region
  for(unsigned int i = 0; i < labelMap->GetNumberOfLabelObjects(); i++)
    {
    // Get the ith region
    LabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(i);
    ImageType3D::RegionType region = labelObject->GetBoundingBox();

    // Get its centroid and the correponding node
    const LabelObjectType::CentroidType centroid = labelObject->GetCentroid();
    ImageType3D::IndexType centroidIdx;
    this->InputImage->TransformPhysicalPointToIndex( centroid, centroidIdx );
    int prevSize = this->IndexToNodeMap.size();
    Node *centroidNode = this->IndexToNodeMap[centroidIdx];
    if(this->IndexToNodeMap.size() != prevSize)
      {
      std::cout << "ERROR with soma centroid lookup..." << std::endl;
      }

    //iterate over this region in NDX image
    //
    itk::ImageRegionConstIteratorWithIndex<ImageType3D>
      itr(this->NDXImage, region);
    while(!itr.IsAtEnd())
      {
      if(itr.Get() == 0.0)
        {
        ++itr;
        continue;
        }
      itk::Index<3> idx = itr.GetIndex();
      if(idx == centroidIdx)
        {
        ++itr;
        continue;
        }
      Node *currentNode = this->IndexToNodeMap[idx];
      std::vector< Node * >::iterator childItr;

      //we only really care if this node is in the Closed list,
      //as this means it's been added to a tree
      childItr =
        std::find(this->Closed.begin(), this->Closed.end(), currentNode);
      if(childItr == this->Closed.end())
        {
        ++itr;
        continue;
        }

      bool keepThisNode = false;
      for(childItr = currentNode->children.begin();
          childItr != currentNode->children.end();
          childItr++)
        {
        //we only keep soma nodes that are the parent of a non-soma node
        if( !region.IsInside( (*childItr)->index ) )
          {
          keepThisNode = true;
          break;
          }
        }

      //keepers are reassigned to be children of the centroid
      if(keepThisNode)
        {
        childItr =
        std::find(centroidNode->children.begin(), centroidNode->children.end(),
                  currentNode);
        if(childItr == centroidNode->children.end())
          {
          centroidNode->children.push_back(currentNode);
          currentNode->parent = centroidNode;
          }
        currentNode->type = 1;
        somaBranches.push_back(currentNode);
        }

      //the rest are carefully deleted
      else
        {
        //assign all the children to their grandparent
        Node *grandparent = currentNode->parent;
        for(childItr = currentNode->children.begin();
            childItr != currentNode->children.end();
            childItr++)
          {
          //don't bother if it's been adopted by someone else
          if( (*childItr)->parent->ID != currentNode->ID )
            {
            continue;
            }
          std::vector< Node *>::iterator rmItr =
            std::find(grandparent->children.begin(), grandparent->children.end(),
                      (*childItr));
          if(rmItr == grandparent->children.end())
            {
            (*childItr)->parent = grandparent;
            grandparent->children.push_back( (*childItr) );
            }
          }

        //remove the node from its parent's list of children
        std::vector< Node* >::iterator rmItr =
          std::remove(grandparent->children.begin(), grandparent->children.end(), currentNode);
        grandparent->children.erase(rmItr, grandparent->children.end());

        //remove the node from Open & Closed
        rmItr =
          std::remove(this->Open.begin(), this->Open.end(), currentNode);
        this->Open.erase(rmItr, this->Open.end());
        rmItr =
          std::remove(this->Closed.begin(), this->Closed.end(), currentNode);
        this->Closed.erase(rmItr, this->Closed.end());
        
        //now we can actually call delete
        delete currentNode;
        //also set the pointer equal to NULL
        currentNode = NULL;
        }
      ++itr;
      }
    }

  this->PruneSomaBranches(somaBranches);
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::PruneSomaBranches(std::vector< Node *> branches)
{
  std::vector< Node * >::iterator branchItr;
  Node *branch;
  for(branchItr = branches.begin(); branchItr != branches.end(); ++branchItr)
    {
    branch = (*branchItr);
    bool anyBranchPts = this->AnyBranchPoints( branch );
    if(!anyBranchPts)
      {
      unsigned int depth = this->GetPathDepth( branch );
      if( depth < 10)
        {
        this->DeleteBranch(branch, true);
        }
      }
    }
}

////////////////////////////////////////////////////////////////////////////////
bool MicrogliaProcessTracer::AnyBranchPoints( Node *n )
{
  if(n->children.size() > 1)
    {
    return true;
    }

  std::vector< Node * >::iterator childItr;
  for(childItr = n->children.begin(); childItr != n->children.end(); ++childItr)
    {
    if(this->AnyBranchPoints( (*childItr) ))
      {
      return true;
      }
    }
 
  return false;
}

////////////////////////////////////////////////////////////////////////////////
unsigned int MicrogliaProcessTracer::GetPathDepth( Node *n )
{
  if(n->children.size() == 0)
    {
    return 1;
    }

  unsigned int maxChildDepth = 1;
  std::vector< Node * >::iterator childItr;
  for(childItr = n->children.begin(); childItr != n->children.end(); ++childItr)
    {
    unsigned int childDepth = this->GetPathDepth( (*childItr) );
    if(childDepth > maxChildDepth)
      {
      maxChildDepth = childDepth;
      }
    }
  return maxChildDepth + 1;
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::DeleteBranch( Node *n, bool parentSurvives )
{
  std::vector< Node *>::iterator rmItr;

  //if your parent is surviving, remove yourself from their list of children
  Node *parent = n->parent;
  if(parentSurvives)
    {
    rmItr =
      std::remove(parent->children.begin(), parent->children.end(), n);
    parent->children.erase(rmItr, parent->children.end());
    }

  //delete your kids
  std::vector< Node * >::iterator childItr;
  for(childItr = n->children.begin(); childItr != n->children.end(); ++childItr)
    {
    this->DeleteBranch( (*childItr), false );
    }

  //remove yourself from the Open and/or Closed lists
  rmItr =
    std::remove(this->Open.begin(), this->Open.end(), n);
  this->Open.erase(rmItr, this->Open.end());
  rmItr =
    std::remove(this->Closed.begin(), this->Closed.end(), n);
  this->Closed.erase(rmItr, this->Closed.end());

  //bye bye
  delete n;
  n = NULL;
}
