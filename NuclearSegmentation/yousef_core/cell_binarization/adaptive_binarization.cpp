/// latest code



#include "adaptive_binarization.h"


double computePoissonProb(int intensity, double alpha)
{
	/*this is the equation
	P = (alpha^intensity)*exp(-alpha)/factorial(intensity);
	however, since alpha and the intensity could be large, computing P in that
	way will result in infinity values from (alpha^intensity) and
	factorial(intensity) as a result of matlab's limitations of the data types*/

	//here is the solution
	double A, P;
	A = exp(-alpha);
	P = 1;
	for (int i=1; i<= intensity; ++i)
	{
		P = P * (alpha/i);
	}

	P = P*A;

	if(P < std::numeric_limits<long double>::epsilon())
		P = std::numeric_limits<long double>::epsilon();

	return P;
}


binaryImageType::Pointer Adaptive_Binarization(inputImageType::Pointer _inputImage)
{

	// READ THE IMAGE
	//if( argc < 2 )
	//{
	//	std::cerr << "Usage: " << std::endl;
	//	std::cerr << argv[0] << " inputImageFile" << std::endl;
	//	return EXIT_FAILURE;
	//}

	// ALSO WORKD FOR UNSIGNED SHORT



	//ReaderType::Pointer reader = ReaderType::New();
	////reader->SetFileName("C:\\DATA\\TEST_REGISTRATION_16_BIT\\montage_kt11306_w410DAPIdsu_new_BS_norm_crop_crop.mhd");
	//reader->SetFileName("C:\\DATA\\TEST_REGISTRATION_16_BIT\\montage_kt11306_w410DAPIdsu_croptest2.mhd");
	//
	//reader->Update();
	//inputImageType::Pointer montage = reader->GetOutput();
	//itk::SizeValueType size_montage[3];
	//size_montage[0] = montage->GetLargestPossibleRegion().GetSize()[0];
	//size_montage[1] = montage->GetLargestPossibleRegion().GetSize()[1];
	//size_montage[2] = montage->GetLargestPossibleRegion().GetSize()[2];

	//typedef  itk::ImageFileWriter< inputImageType  > WriterType3;
	//WriterType3::Pointer writer3 = WriterType3::New();
	//writer3->SetFileName("out.mhd");
	//writer3->SetInput(montage);
	//writer3->Update();

	//binaryImageType::Pointer binMontage = binaryImageType::New();
	//itk::Size<3> im_size;
	//im_size[0] = size_montage[0]; 
	//im_size[1] = size_montage[1];    
	//im_size[2] = size_montage[2];  
	//binaryImageType::IndexType start;
	//start[0] =   0;  // first index on X
	//start[1] =   0;  // first index on Y    
	//start[2] =   0;  // first index on Z  
	//binaryImageType::PointType origin;
	//origin[0] = 0; 
	//origin[1] = 0;    
	//origin[2] = 0;    
	//binMontage->SetOrigin( origin );
	//binaryImageType::RegionType region_montage;
	//region_montage.SetSize( im_size );
	//region_montage.SetIndex( start );
	//binMontage->SetRegions( region_montage );
	//binMontage->Allocate();
	//binMontage->FillBuffer(0);
	//binMontage->Update();
	//binaryImageType::PixelType * binArray = binMontage->GetBufferPointer();
	//int slice_size = im_size[1] * im_size[0];
	//int row_size = im_size[0];

	//int num_rows = (int)ceil((double)size_montage[1]/450);
	//int num_cols = (int)ceil((double)size_montage[0]/450);

	//for(int row=0; row<num_rows; ++row)
	//{
	//	for(int col=0; col<num_cols; ++col)
	//	{

	//		std::cout<<"Extracting Tile " << col*450 << "_" << row*450 << "\n";			
	//		inputImageType::IndexType start_tile;
	//		start_tile[0] = col*450;
	//		start_tile[1] = row*450;
	//		start_tile[2] = 0;

	//		inputImageType::SizeType size_tile;
	//		size_tile[0] = (((col*450)+500) < size_montage[0]) ? 500 : (size_montage[0] - (col*450));
	//		size_tile[1] = (((row*450)+500) < size_montage[1]) ? 500 : (size_montage[1] - (row*450));
	//		size_tile[2] = size_montage[2];

	//		inputImageType::RegionType region_tile;
	//		region_tile.SetSize(size_tile);
	//		region_tile.SetIndex(start_tile);

	//		ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
	//		ROIfilter->SetRegionOfInterest(region_tile);
	//		ROIfilter->SetInput(montage);
	//		ROIfilter->Update();
	//		inputImageType::Pointer _inputImage = ROIfilter->GetOutput();

			//ReaderType::Pointer reader = ReaderType::New();
			//reader->SetFileName(argv[1]);
			//reader->Update();

			//inputImageType::Pointer _inputImage = reader->GetOutput();

			inputImageType::PixelType 		*_inputImageArray; 
			_inputImageArray = _inputImage->GetBufferPointer();


			int _numColumns = _inputImage->GetLargestPossibleRegion().GetSize()[0];
			int _numRows = _inputImage->GetLargestPossibleRegion().GetSize()[1];
			int _numStacks = _inputImage->GetLargestPossibleRegion().GetSize()[2];
			long long _totNumPixels = (long long)_numRows*(long long)_numColumns*(long long)_numStacks;





			// PARAMETERS
			// 	int _numberBins_mixPoisson = 65536;
			int _numberBins_mixPoisson = 256;
			int _sigmaNeighCost = 2000; // !! This should be changed
			int _wNeigh = 2;

			int _windowSize = 60; // the size will be _windowSize*2+1

			int win = 3;

			int _connComponentsConnectivity = 26;
			int _minObjectSize = 10;

			int banderaHistoGlobal = 0;
			int banderaShowGlobalHistogram = 0;
			int bandera3Poisson = 2; // 2 3 distribution, 2 2 distributions

			inputPixelType _maxValueInputPixelType = std::numeric_limits<inputPixelType>::max();
			unsigned char _maxValueBinaryPixelType = std::numeric_limits<unsigned char>::max();

			double _spacing_Z = 2;


			//std::cout << std::endl << _wNeigh*exp(-pow(10000,2)/(2*pow(_sigmaNeighCost,2)));

			// 	int yyt;
			// 	std::cin >> yyt;

//
			// CREATE THE BINARY IMAGE

			binaryImageType::Pointer _binaryImage = binaryImageType::New();
			binaryImageType::IndexType start;
			start[0] = 0;
			start[1] = 0;
			start[2] = 0;
			binaryImageType::SizeType size;
			size[0] = _numColumns;
			size[1] = _numRows;
			size[2] = _numStacks;
			binaryImageType::RegionType region;
			region.SetIndex( start );
			region.SetSize( size );

			double spacing[3];
			spacing[0] = 1; //spacing along x
			spacing[1] = 1; //spacing along y
			spacing[2] = _spacing_Z; //spacing along z

			_binaryImage->SetRegions( region );
			_binaryImage->SetSpacing( spacing );

			try{
				_binaryImage->Allocate();
			}
			catch( itk::ExceptionObject & err ){
				std::cerr << "ExceptionObject caught!" << std::endl;
				std::cerr << err << std::endl;
			}
			binaryImageType::PixelType ceros = 0;
			_binaryImage->FillBuffer( ceros );
			_binaryImage->Update();

			binaryImageType::PixelType 		*_binaryImageArray; 
			_binaryImageArray = _binaryImage->GetBufferPointer();




			// CREATE A DOUBLE IMAGE FOR THE COST FUNCTIONS
			typedef itk::Image< double, 3 >   costImageType;

			costImageType::Pointer _costImageF = costImageType::New();
			costImageType::IndexType startF;
			startF[0] = 0;
			startF[1] = 0;
			startF[2] = 0;
			costImageType::SizeType sizeF;
			sizeF[0] = _numColumns;
			sizeF[1] = _numRows;
			sizeF[2] = _numStacks;
			costImageType::RegionType regionF;
			regionF.SetIndex( startF );
			regionF.SetSize( sizeF );

			double spacingF[3];
			spacingF[0] = 1; //spacing along x
			spacingF[1] = 1; //spacing along y
			spacingF[2] = _spacing_Z; //spacing along z

			_costImageF->SetRegions( regionF );
			_costImageF->SetSpacing( spacingF );

			try{
				_costImageF->Allocate();
			}
			catch( itk::ExceptionObject & err ){
				std::cerr << "ExceptionObject caught!" << std::endl;
				std::cerr << err << std::endl;
			}
			costImageType::PixelType cerosF = 0;
			_costImageF->FillBuffer( cerosF );
			_costImageF->Update();

			costImageType::PixelType 		*_costImageFArray; 
			_costImageFArray = _costImageF->GetBufferPointer();



			costImageType::Pointer _costImageB = costImageType::New();
			costImageType::IndexType startB;
			startB[0] = 0;
			startB[1] = 0;
			startB[2] = 0;
			costImageType::SizeType sizeB;
			sizeB[0] = _numColumns;
			sizeB[1] = _numRows;
			sizeB[2] = _numStacks;
			costImageType::RegionType regionB;
			regionB.SetIndex( startB );
			regionB.SetSize( sizeB );

			double spacingB[3];
			spacingB[0] = 1; //spacing along x
			spacingB[1] = 1; //spacing along y
			spacingB[2] = _spacing_Z; //spacing along z

			_costImageB->SetRegions( regionB );
			_costImageB->SetSpacing( spacingB );

			try{
				_costImageB->Allocate();
			}
			catch( itk::ExceptionObject & err ){
				std::cerr << "ExceptionObject caught!" << std::endl;
				std::cerr << err << std::endl;
			}
			costImageType::PixelType cerosB = 0;
			_costImageB->FillBuffer( cerosB );
			_costImageB->Update();

			costImageType::PixelType 		*_costImageBArray; 
			_costImageBArray = _costImageB->GetBufferPointer();


			// 	costImageType::Pointer _costImageAB = costImageType::New();
			// 	 costImageType::IndexType startAB;
			// 	startAB[0] = 0;
			// 	startAB[1] = 0;
			// 	startAB[2] = 0;
			// 	 costImageType::SizeType sizeAB;
			// 	sizeAB[0] = _numColumns;
			// 	sizeAB[1] = _numRows;
			// 	sizeAB[2] = _numStacks;
			// 	 costImageType::RegionType regionAB;
			// 	regionAB.SetIndex( startAB );
			// 	regionAB.SetSize( sizeAB );
			// 	
			// 	double spacingAB[3];
			// 	spacingAB[0] = 1; //spacing along x
			// 	spacingAB[1] = 1; //spacing along y
			// 	spacingAB[2] = _spacing_Z; //spacing along z
			// 
			// 	_costImageAB->SetRegions( regionAB );
			// 	_costImageAB->SetSpacing( spacingAB );
			// 	
			// 	try{
			// 		_costImageAB->Allocate();
			// 	}
			// 	catch( itk::ExceptionObject & err ){
			// 		std::cerr << "ExceptionObject caught!" << std::endl;
			// 		std::cerr << err << std::endl;
			// 	}
			// 	   costImageType::PixelType cerosAB = 0;
			// 	_costImageB->FillBuffer( cerosAB );
			// 	_costImageB->Update();
			// 	
			// 	 costImageType::PixelType 		*_costImageABArray; 
			// 	_costImageABArray = _costImageAB->GetBufferPointer();








			// COMMON VARIABLES

			int num_nodes;
			int num_edges;





			// 	long long curr_node_z_inside;
			// 	long long curr_node_yz_inside;
			// 	long long curr_node_xyz_inside;
			// 	
			// 	
			// 	
			// 	long long rght_node;
			// 	long long down_node;
			// 	long long diag_node;
			// 	
			// 	long long diag_node_1;
			// 	long long rght_node_1;
			// 	long long down_node_1;
			// 	long long diag_node_2;
			// 	
			// 	
			// 	double Df;
			// 	double Db;
			// 	double Dr;
			// 	double Dd; 
			// 	double Dg; 
			// 	double F_H[_maxValueInputPixelType+1];
			// 	double B_H[_maxValueInputPixelType+1];
			// 	
			// 	double _alpha_1;
			// 	double _P_I;
			// 	double _alpha_2;
			// 	double _P_I2;
			// 	double _alpha_3;
			// 	int _numPoissonDist;



			typedef Graph_B<int,int,int> GraphType;
			// No much diference in the result when using double, the memory increases significatively. 
			// 	typedef Graph<double,double,double> GraphType;

			// 	Set the number of edges and the number of nodes
			num_nodes = _numRows*_numColumns*_numStacks;
			if( _numStacks == 1 )
			{
				num_edges = (_numRows-1)*(_numColumns-1)*3;
			}
			else
			{
				num_edges = (_numRows-1)*(_numColumns-1)*(_numStacks-1)*3;
				// In case of using 26 neighborhood
				// 		num_edges = (_numRows-1)*(_numColumns-1)*(_numStacks-1)*7;
			}


			// 	std::cout << std::endl << "	Allocating:	" << num_nodes << " " << num_nodes*sizeof(int) / (1024.0 * 1024.0) << " MB of memory, for nodes.";
			// 	std::cout << std::endl << "	And, allocating:	" << num_edges << " " << num_edges*sizeof(int) / (1024.0 * 1024.0) << " MB of memory, for edges.";

			GraphType *g = new GraphType(/*estimated # of nodes*/ num_nodes, /*estimated # of edges*/ num_edges); 






			double P0_Global, U0_Global, P1_Global, U1_Global, P2_Global, U2_Global, U_Global, J_Global, min_J_Global;
			min_J_Global = 1000000.0;


			double _alpha_1_Global;
			double _P_I_Global;
			double _alpha_2_Global;
			double _P_I2_Global;
			double _alpha_3_Global;
			int _numPoissonDist_Global;

			//double histoGramGlobal[_numberBins_mixPoisson];
			std::vector<double> histoGramGlobal;
			histoGramGlobal.resize(_numberBins_mixPoisson);

			for( unsigned int mm=0; mm<_numberBins_mixPoisson; ++mm )
			{
				histoGramGlobal[mm] = 0.0;
			}



			// GLOBAL HISTOGRAM
			if( banderaHistoGlobal == 1 )
			{
				std::cout << std::endl << "Lets Compute Histogram";
				std::cout << std::endl << "Lets Compute Histogram";

				for( int k=0; k<_numStacks; ++k )
				{
					// 			std::cout << std::endl << "Stack: " << k << " of " << _numStacks << " Checked";
					for( int i=0; i<_numRows; ++i )
					{
						// 				printf('\r 5\n');
						// 				printf("\r|%f",i);
						// 				std::cout << std::endl << "	Row: " << i << " of " << _numRows;

						//#pragma omp parallel for
						for(int j=0; j<_numColumns; ++j )
						{	
							long long curr_node_z_outside = k*_numColumns*_numRows;
							long long curr_node_yz_outside = (i*_numColumns)+curr_node_z_outside;
							long long curr_node_xyz_outside = j + curr_node_yz_outside;

							double factorBins = ((double)_maxValueInputPixelType+1)/((double)_numberBins_mixPoisson); // +1 because is unsigned char (255), inputPixelType (65535), for other data types have to be different (unsigned in general +1)
							int binIndex = floor(((double)_inputImageArray[curr_node_xyz_outside]/factorBins));	// !!!!!! DIVIDE THE IMAGE INTO TILES
							if( binIndex >= _numberBins_mixPoisson )
							{
								binIndex = _numberBins_mixPoisson-1;
							}
							histoGramGlobal[binIndex]++;
						}
						// 				printf("\r");
						// 				std::cout << "\r";
					}
				}

				// 		std::cout << std::endl << "Histogram computed";
				std::cout << std::endl << "Histogram computed";

				for( unsigned int nn=0; nn<_numberBins_mixPoisson; ++nn )
				{
					histoGramGlobal[nn] /= _numStacks*_numColumns*_numRows;
					double factorBins = ((double)_maxValueInputPixelType+1)/((double)_numberBins_mixPoisson); // +1 because is unsigned char (255), inputPixelType (65535), for other data types have to be different (unsigned in general +1)
					if( banderaShowGlobalHistogram == 1 )
					{
						std::cout << "\n" << nn << "\t" << (double)histoGramGlobal[nn] << " " << factorBins;
					}

				}


				for(int pp=0; pp<(_numberBins_mixPoisson-1); ++pp)//to set the first threshold
				{
					//compute the current parameters of the first component
					P0_Global = U0_Global = 0.0;		
					for(int l=0; l<=pp; l++)
					{
						P0_Global+=histoGramGlobal[l];
						U0_Global+=(l+1)*histoGramGlobal[l];
					}
					U0_Global /= P0_Global;

					for(int qq=pp+1; qq<_numberBins_mixPoisson; qq++)//to set the second threshold
					{
						//compute the current parameters of the second component
						P1_Global = U1_Global = 0.0;		
						for(int l=pp+1; l<=qq; l++)
						{
							P1_Global+=histoGramGlobal[l];
							U1_Global+=(l+1)*histoGramGlobal[l];
						}
						U1_Global /= P1_Global;

						//compute the current parameters of the third component
						P2_Global = U2_Global = 0.0;		
						for(int l=qq+1; l<=255; l++)
						{
							P2_Global+=histoGramGlobal[l];
							U2_Global+=(l+1)*histoGramGlobal[l];
						}
						U2_Global /= P2_Global;

						//compute the overall mean
						U_Global = P0_Global*U0_Global + P1_Global*U1_Global + P2_Global*U2_Global;

						//Compute the current value of the error criterion function
						J_Global =  U_Global - (P0_Global*(log(P0_Global)+U0_Global*log(U0_Global))+ P1_Global*(log(P1_Global)+U1_Global*log(U1_Global)) + P2_Global*(log(P2_Global)+U2_Global*log(U2_Global)));
						//Add the penalty term
						// 						J +=PenTerm3;

						if(J_Global<min_J_Global)
						{

							min_J_Global = J_Global;
							_alpha_1_Global = U0_Global;
							_P_I_Global = P0_Global;
							_alpha_2_Global = U1_Global;
							_P_I2_Global = P1_Global;
							_alpha_3_Global = U2_Global; //Just a negative number to let the program knows that two levels will be used		
							_numPoissonDist_Global = 3;
						}
					}
				}

				std::cout << std::endl << "Global hist params: " << _alpha_1_Global << ", " << _alpha_2_Global << ", " << _alpha_3_Global << ", " << _P_I_Global << ", " << _P_I2_Global;
			}






			// FOR EACH PIXEL



			int modBig = (2*win+1);
			int modSmall = (win+1);
			int k_end_big = win + modBig*((int)(_numStacks-win-1)/modBig);
			int i_end_big = win + modBig*((int)(_numRows-win-1)/modBig);
			int j_end_big = win + modBig*((int)(_numColumns-win-1)/modBig);


			// std::cout << std::endl << _numStacks << " " << _numRows << " " << _numColumns;
			// std::cout << std::endl << modBig << " " << modSmall << " " << k_end_big << " " << i_end_big << " " << j_end_big;
			// int rt;
			// std::cin >> rt;
			//#pragma omp parallel for collapse
			for( int k=win; k<=k_end_big; k=k+2*win+1 )
			{
				// 		std::cout << std::endl << "Stack: " << k << " of " << _numStacks;
				// 		long long curr_node_z_outside = k*_numColumns*_numRows;
#pragma omp parallel for
				for(int i=win; i<=i_end_big; i=i+2*win+1)
				{
					// Temporary it should go up
					// 			long long curr_node_z_outside = k*_numColumns*_numRows;

					// 			std::cout << std::endl << "Stack: " << k << " of " << _numStacks;
					// 			std::cout << std::endl << "	Row: " << i << " of " << _numRows;

					// 			long long curr_node_yz_outside = (i*_numColumns)+curr_node_z_outside;
#pragma omp parallel for
					for(int j=win; j<=j_end_big; j=j+2*win+1)
					{	
						// 				std::cout << std::endl << "Stack: " << k << " of " << _numStacks;
						// 				std::cout << std::endl << "	Row: " << i << " of " << _numRows;
						// 				std::cout << std::endl << "	Col: " << j << " of " << _numColumns;


						long long curr_node_z_outside = k*_numColumns*_numRows;
						long long curr_node_yz_outside = (i*_numColumns)+curr_node_z_outside;
						long long curr_node_xyz_outside = j + curr_node_yz_outside;




						// 				int counter_small = 0;
						// 				for( int k_small = k-win; k_small<=k+win; ++k_small )
						// 				{
						// 					long long curr_node_z_outside_small = k_small*_numColumns*_numRows;
						// 					for( int i_small = i-win; i_small<=i+win; ++i_small )
						// 					{
						// 						long long curr_node_yz_outside_small = (i_small*_numColumns)+curr_node_z_outside_small;
						// 						for( int j_small = j-win; j_small<=j+win; ++j_small )
						// 						{
						// 							long long curr_node_xyz_outside_small = j_small + curr_node_yz_outside_small;
						// 				
						// 							inputPixelType value_test = _inputImageArray[curr_node_xyz_outside_small];
						// 							
						// 							
						// 							if( value_test < 15 )
						// 							{
						// 								_costImageBArray[curr_node_xyz_outside_small] = 0;
						// 								_costImageFArray[curr_node_xyz_outside_small] = 1000;
						// 								counter_small++;
						// // 								continue;
						// 							}
						// 							else if( value_test > 190 )
						// 							{
						// 								_costImageBArray[curr_node_xyz_outside_small] = 1000;
						// 								_costImageFArray[curr_node_xyz_outside_small] = 0;
						// 								counter_small++;
						// // 								continue;
						// 							}
						// // 							if( (k_small == 9) && (i_small == 45) && (j_small == 40) )
						// 							if( (int)value_test == 26 )
						// 							{
						// 								std::cout << std::endl << "value_1: " << (int)value_test;
						// 								std::cout << std::endl << _costImageBArray[curr_node_xyz_outside_small] << ", " << _costImageFArray[curr_node_xyz_outside_small];
						// 							}
						// 				
						// 						}
						// 					}
						// 				}
						// 				
						// 				if( counter_small == 27 )
						// 					continue;



						long long curr_node_z_inside;
						long long curr_node_yz_inside;
						long long curr_node_xyz_inside;





						long long diag_node_1;
						long long rght_node_1;
						long long down_node_1;
						long long diag_node_2;



						double P0, U0, P1, U1, P2, U2, U, J, min_J;
						min_J = 1000000.0;


						double _alpha_1;
						double _P_I;
						double _alpha_2;
						double _P_I2;
						double _alpha_3;
						int _numPoissonDist;


						double factorBins = ((double)_maxValueInputPixelType+1)/((double)_numberBins_mixPoisson); // +1 because is unsigned char (255), inputPixelType (65535), for other data types have to be different (unsigned in general +1)
						unsigned int binIndex;


						// Create a normalized histogram
						std::vector<double> histoGram;
						histoGram.resize(_numberBins_mixPoisson);
						//double histoGram[_numberBins_mixPoisson];

						if( banderaHistoGlobal == 0 )
						{

							for( unsigned int mm=0; mm<_numberBins_mixPoisson; ++mm )
							{
								histoGram[mm] = 0.0;
							}

							// Create histogram





							int k_ini = k - _windowSize;
							if( k_ini<0 )
								k_ini = 0;
							int k_end = k + _windowSize;
							if( k_end>_numStacks )
								k_end = _numStacks;

							int i_ini = i - _windowSize;
							if( i_ini<0 )
								i_ini = 0;
							int i_end = i + _windowSize;
							if( i_end>_numRows )
								i_end = _numRows;

							int j_ini = j - _windowSize;
							if( j_ini<0 )
								j_ini = 0;
							int j_end = j + _windowSize;
							if( j_end>_numColumns )
								j_end = _numColumns;


							long long _totNumPixels_inside = (k_end-k_ini)*(i_end-i_ini)*(j_end-j_ini);

							for( int k_inside=k_ini; k_inside<k_end; ++k_inside )
							{
								curr_node_z_inside = k_inside*_numColumns*_numRows;
								for(int i_inside=i_ini; i_inside<i_end; ++i_inside)
								{
									curr_node_yz_inside = (i_inside*_numColumns)+curr_node_z_inside;
									for(int j_inside=j_ini; j_inside<j_end; ++j_inside)
									{	
										// 							std::cout << std::endl << k_inside << " " << i_inside << " " << j_inside;
										curr_node_xyz_inside = j_inside + curr_node_yz_inside;

										// 							std::cout << std::endl << (double)_inputImageArray[curr_node_xyz_inside] << " " << floor(((double)_inputImageArray[curr_node_xyz_inside]/factorBins));
										binIndex = floor(((double)_inputImageArray[curr_node_xyz_inside]/factorBins));	// !!!!!! DIVIDE THE IMAGE INTO TILES
										if( binIndex >= _numberBins_mixPoisson )
										{
											binIndex = _numberBins_mixPoisson-1;
										}
										histoGram[binIndex]++;
										// 								_totNumPixels_inside++;
									}
								}
							}

							// 				std::cout << std::endl << "Histogram Computed";
							// 				std::cout << std::endl << "Histogram Computed";


							// Normalize histogram
							// 				std::cout << std::endl;
							for( unsigned int nn=0; nn<_numberBins_mixPoisson; ++nn )
							{
								// 					std::cout << std::endl << histoGram[nn] << " " << _totNumPixels_inside << " " << factorBins << " " << (double)_maxValueInputPixelType;
								histoGram[nn] /= _totNumPixels_inside;
								// 					162 187 312
								// 						if( (k == 117) && (i == 262) && (j == 192) )
								// 						{
								// 							std::cout << "\n" << nn << "\t" << (double)histoGram[nn] << " " << factorBins;
								// 						}

							}

							// 				std::cout << std::endl << "Histogram Normalize";
							// 				std::cout << std::endl << "Histogram Normalize";
							// 				int yt;
							// 				std::cin>>yt;

							if( bandera3Poisson == 2 )
							{


								//The three-level min error thresholding algorithm

								std::vector<double> histoGramSum;
								histoGramSum.resize(_numberBins_mixPoisson);
								//double histoGramSum[_numberBins_mixPoisson];
								histoGramSum[0] = histoGram[0];

								std::vector<double> histoGramSumWeight;
								histoGramSumWeight.resize(_numberBins_mixPoisson);
								//double histoGramSumWeight[_numberBins_mixPoisson];
								histoGramSumWeight[0] = (0+1)*histoGram[0];

								std::vector<double> histoGramSumWeight2;
								histoGramSumWeight2.resize(_numberBins_mixPoisson);
								//double histoGramSumWeight2[_numberBins_mixPoisson];
								histoGramSumWeight2[_numberBins_mixPoisson-1] = (_numberBins_mixPoisson)*histoGram[_numberBins_mixPoisson-1];

								for(int pp=1; pp<_numberBins_mixPoisson; ++pp){
									histoGramSum[pp] = histoGramSum[pp-1] + histoGram[pp];
									histoGramSumWeight[pp]= histoGramSumWeight[pp-1] + (pp+1)*histoGram[pp];
									histoGramSumWeight2[_numberBins_mixPoisson-1-pp]=histoGramSumWeight2[_numberBins_mixPoisson-pp] + (_numberBins_mixPoisson-pp)*histoGram[_numberBins_mixPoisson-1-pp];
								}


								for(int pp=0; pp<(_numberBins_mixPoisson-1); ++pp)//to set the first threshold
								{
									P0 = histoGramSum[pp];
									P1 = 1-P0;

									U0 = histoGramSumWeight[pp]/P0;
									U1 = histoGramSumWeight2[pp+1]/P1;

									// // 						std::cout<<std::endl << P0 << "\t" << P1 << "\t" << U0 << "\t" << U1 << "\t" << pp;
									// // 						std::cout<<std::endl << histoGramSum[pp] << "\t" << P1 << "\t" << histoGramSumWeight[pp]/P0 << "\t" << histoGramSumWeight2[pp+1]/P1<<"\n";
									// 	// 					std::cout << std::endl << "P: " << pp;
									// 						//compute the current parameters of the first component
									// 						P0 = U0 = 0.0;		
									// 						for(int l=0; l<=pp; ++l)
									// 						{
									// 	// 						std::cout << std::endl << "L: " << l;
									// 							P0+=histoGram[l];
									// 							U0+=(l+1)*histoGram[l];
									// 						}
									// 	// 					if( P0 == 0 )
									// 	// 						std::cout << std::endl << "PP:		" << pp;
									// 						
									// 						
									// 						U0 /= P0;
									// 	
									// 						U1 = 0.0;		
									// 						P1 = 1-P0;
									// 						for(int qq=pp+1; qq<_numberBins_mixPoisson; ++qq)//to set the second threshold
									// 						{
									// 							//compute the current parameters of the second component
									// 							
									// 				// 			for(int l=j; l<=_numberBins_mixPoisson; ++l)
									// 				// 			{
									// // 								P1+=histoGram[qq];
									// 								U1+=(qq+1)*histoGram[qq];
									// 				// 			}
									// 							
									// 						}
									// 						U1 /= P1;

									//compute the overall mean
									U = P0*U0 + P1*U1;

									//Compute the current value of the error criterion function
									J =  U - (P0*(log(P0)+U0*log(U0))+ P1*(log(P1)+U1*log(U1)));
									//Add the penalty term
									// 		J +=PenTerm2;
									if(J<min_J)
									{
										min_J = J;
										_alpha_1 = U0;
										_P_I = P0;
										_alpha_2 = U1;
										_P_I2 = P1;
										_alpha_3 = U2; //Just a negative number to let the program knows that two levels will be used		
										_numPoissonDist = 2;

										// 						std::cout<<std::endl << P0 << "\t" << P1 << "\t" << U0 << "\t" << U1 << "\t" << pp;
										// 						std::cout<<std::endl << histoGramSum[pp] << "\t" << P1 << "\t" << histoGramSumWeight[pp]/P0 << "\t" << histoGramSumWeight2[pp+1]/P1<<"\n";

										// 						std::cout << " " << pp;


										// 						std::cout << std::endl << "PP:		" << P0 << " " << _P_I;
										// 						int yy;
										// 						std::cin>>yy;
									}
								}

							}
							else
							{

								std::vector<double> histoGramSum;
								histoGramSum.resize(_numberBins_mixPoisson);
								//double histoGramSum[_numberBins_mixPoisson];
								histoGramSum[0] = histoGram[0];

								std::vector<double> histoGramSum2;
								histoGramSum2.resize(_numberBins_mixPoisson);
								//double histoGramSum2[_numberBins_mixPoisson];
								histoGramSum2[_numberBins_mixPoisson-1] = histoGram[_numberBins_mixPoisson-1];

								std::vector<double> histoGramSumWeight;
								histoGramSumWeight.resize(_numberBins_mixPoisson);
								//double histoGramSumWeight[_numberBins_mixPoisson];
								histoGramSumWeight[0] = (0+1)*histoGram[0];

								std::vector<double> histoGramSumWeight2;
								histoGramSumWeight2.resize(_numberBins_mixPoisson);
								//double histoGramSumWeight2[_numberBins_mixPoisson];
								histoGramSumWeight2[_numberBins_mixPoisson-1] = (_numberBins_mixPoisson)*histoGram[_numberBins_mixPoisson-1];

								for(int pp=1; pp<_numberBins_mixPoisson; ++pp){
									histoGramSum[pp] = histoGramSum[pp-1] + histoGram[pp];
									histoGramSum2[_numberBins_mixPoisson-1-pp] = histoGramSum2[_numberBins_mixPoisson-pp] + histoGram[_numberBins_mixPoisson-1-pp];
									histoGramSumWeight[pp]= histoGramSumWeight[pp-1] + (pp+1)*histoGram[pp];
									histoGramSumWeight2[_numberBins_mixPoisson-1-pp]=histoGramSumWeight2[_numberBins_mixPoisson-pp] + (_numberBins_mixPoisson-pp)*histoGram[_numberBins_mixPoisson-1-pp];
								}
								double P0_test,P1_test,P2_test;
								double U0_test,U1_test,U2_test;
								double U0_ini, U1_ini, U2_ini;
								// CAREFUL WITH THIS TEST
								for(int pp=0; pp<(_numberBins_mixPoisson-1)/2; ++pp)//to set the first threshold
								{
									P0_test = histoGramSum[pp];
									U0_ini = histoGramSumWeight[pp];
									U0_test = U0_ini/P0_test;

									// 						for(int qq=pp+1; qq<_numberBins_mixPoisson; qq++)//to set the second threshold
									// 						{


									// 							P2_test = histoGramSum2[qq];
									// 							U2_ini = histoGramSumWeight2[qq];
									// 							U2_test = U2_ini/P2_test;
									// 							P1_test = 1 - P1_test - P2_test;
									// 							U1_ini = histoGramSumWeight[_numberBins_mixPoisson-1] - U2_test - U0_test;
									// 							U1_test = U1_ini/P1;


									// 						//compute the current parameters of the first component
									// 						P0 = U0 = 0.0;		
									// 						for(int l=0; l<=pp; l++)
									// 						{
									// 							P0+=histoGram[l];
									// 							U0+=(l+1)*histoGram[l];
									// 						}
									// 						U0 /= P0;

									for(int qq=pp+1; qq<_numberBins_mixPoisson-1; qq++)//to set the second threshold
									{
										P2_test = histoGramSum2[qq+1];
										U2_ini = histoGramSumWeight2[qq+1];
										U2_test = U2_ini/P2_test;
										P1_test = 1 - P0_test - P2_test;
										U1_ini = histoGramSumWeight[_numberBins_mixPoisson-1] - histoGramSumWeight[pp] - histoGramSumWeight2[qq+1];
										U1_test = U1_ini/P1_test;

										// 							//compute the current parameters of the second component
										// 							P1 = U1 = 0.0;		
										// 							for(int l=pp+1; l<=qq; l++)
										// 							{
										// 								P1+=histoGram[l];
										// 								U1+=(l+1)*histoGram[l];
										// 							}
										// 							U1 /= P1;

										// 							//compute the current parameters of the third component
										// 							P2 = U2 = 0.0;		
										// 							for(int l=qq+1; l<_numberBins_mixPoisson; l++)
										// 							{
										// 								P2+=histoGram[l];
										// 								U2+=(l+1)*histoGram[l];
										// 							}
										// 							U2 /= P2;

										//compute the overall mean
										U = P0_test*U0_test + P1_test*U1_test + P2_test*U2_test;

										//Compute the current value of the error criterion function
										J =  U - (P0_test*(log(P0_test)+U0_test*log(U0_test))+ P1_test*(log(P1_test)+U1_test*log(U1_test)) + P2_test*(log(P2_test)+U2_test*log(U2_test)));
										//Add the penalty term
										// 						J +=PenTerm3;

										if(J<min_J)
										{
											/*							min_J = J;
											Alpha1[0] = U0;
											P_I1[0] = P0;
											Alpha2[0] = U1;
											P_I2[0] = P1;
											Alpha3[0] = U2;			*/	

											min_J = J;
											_alpha_1 = U0_test;
											_P_I = P0_test;
											_alpha_2 = U1_test;
											_P_I2 = P1_test;
											_alpha_3 = U2_test; //Just a negative number to let the program knows that two levels will be used		
											_numPoissonDist = 3;

											// 								std::cout << std::endl << P0 << "\t" << U0 << "\t" << P1 << "\t" << U1 << "\t" << P2 << "\t" << U2;
											// 								std::cout << std::endl << P0_test << "\t" << U0_test << "\t" << P1_test << "\t" << U1_test << "\t" << P2_test << "\t" << U2_test << "\n";
										}
									}
								}
							}

						} // Only done one time

						else
						{
							// 					std::cout << std::endl << "NEVER";

							min_J = min_J_Global;
							_alpha_1 = _alpha_1_Global;
							_P_I = _P_I_Global;
							_alpha_2 = _alpha_2_Global;
							_P_I2 = _P_I2_Global;
							_alpha_3 = _alpha_3_Global; //Just a negative number to let the program knows that two levels will be used		
							_numPoissonDist = _numPoissonDist_Global;
						}
						// 				std::cout << std::endl << "SAlio";
						// 				std::cout << std::endl << "SAlio";
						// 				std::cout << std::endl << "SAlio";



						// 				
						int k_ini_small = k - win;
						int k_end_small = k + win;
						if( k == k_end_big )
							k_end_small = _numStacks-1;
						int i_ini_small = i - win;
						int i_end_small = i + win;
						if( i == i_end_big )
							i_end_small = _numRows-1;
						int j_ini_small = j - win;
						int j_end_small = j + win;
						if( j == j_end_big )
							j_end_small = _numColumns-1;

						// 				std::cout << std::endl << "Mix Estimated" << "i: " << i << ", j: " << j  << ", k: " << k << " " << k_ini_small << " " << k_end_small << " " << i_ini_small << " " << i_end_small << " " << j_ini_small << " " << j_end_small;
						// 				std::cout << std::endl << "Mix Estimated";


						if( bandera3Poisson == 2)
						{
							//#pragma omp parallel for
							double F_H_Alpha2 = _P_I2*computePoissonProb((int)_alpha_2,_alpha_2);
							double B_H_Alpha1 = _P_I*computePoissonProb((int)_alpha_1,_alpha_1);
							for( int k_small = k_ini_small; k_small<=k_end_small; ++k_small )
							{
								long long curr_node_z_outside_small = k_small*_numColumns*_numRows;
								for( int i_small = i_ini_small; i_small<=i_end_small; ++i_small )
								{
									long long curr_node_yz_outside_small = (i_small*_numColumns)+curr_node_z_outside_small;
									for( int j_small = j_ini_small; j_small<=j_end_small; ++j_small )
									{
										long long curr_node_xyz_outside_small = j_small + curr_node_yz_outside_small;

										inputPixelType value = _inputImageArray[curr_node_xyz_outside_small];
										double F_H;
										double B_H;

										// Not Sure if is P_! or P_I2
										if(floor(((double)value/factorBins))>=_alpha_2)
											F_H = F_H_Alpha2;
										else
											// 								F_H = (1-_P_I)*computePoissonProb(value,_alpha_2);
											F_H = _P_I2*computePoissonProb(floor(((double)value/factorBins)),_alpha_2);
										if(floor(((double)value/factorBins))<=_alpha_1)
											B_H = B_H_Alpha1;
										else
											// 								B_H = _P_I*computePoissonProb(value,_alpha_1);
											B_H = _P_I*computePoissonProb(floor(((double)value/factorBins)),_alpha_1);


										// 				std::cout << std::endl << F_H << " " << B_H;


										_costImageBArray[curr_node_xyz_outside_small] = -log(B_H);
										// 							if(_costImageBArray[curr_node_xyz_outside_small]>1000.0)
										// 								_costImageBArray[curr_node_xyz_outside_small] = 1000;
										// 							if(_costImageBArray[curr_node_xyz_outside_small]<0)
										// 								_costImageBArray[curr_node_xyz_outside_small] = 0;

										_costImageFArray[curr_node_xyz_outside_small] = -log(F_H);
										// 							if(_costImageFArray[curr_node_xyz_outside_small]>1000.0)
										// 								_costImageFArray[curr_node_xyz_outside_small] = 1000;
										// 							if(_costImageFArray[curr_node_xyz_outside_small]<0)
										// 								_costImageFArray[curr_node_xyz_outside_small] = 0;

										// Just to test that each pixels is only access one time
										// 							_costImageABArray[curr_node_xyz_outside_small] = _costImageABArray[curr_node_xyz_outside_small] + 1;
										if( _alpha_2 < 20 )
										{
											_costImageBArray[curr_node_xyz_outside_small] = 0;
											_costImageFArray[curr_node_xyz_outside_small] = 1000;
										}




										// // 							if( (k_small == 9) && (i_small == 45) && (j_small == 40) )
										// 							if( (int)value == 500 )
										// 							if( _costImageBArray[curr_node_xyz_outside_small] > _costImageFArray[curr_node_xyz_outside_small])
										// 							/*if( (k_small == 1) && (i_small == 372) && (j_small == 746) )
										// 							{
										// 							    std::cout << std::endl << k << " " << i << " " << j;
										// 								std::cout << std::endl << "value: " << (int)value << ", " << _alpha_1 << ", " << _alpha_2 << ", " << _alpha_3 << ", " << _P_I << ", " << _P_I2 << ", " << F_H << ", " << B_H;
										// 								std::cout << std::endl << _costImageFArray[curr_node_xyz_outside_small] << ", " << _costImageBArray[curr_node_xyz_outside_small];
										// 							int yu;
										// 							std::cin>>yu;
										// 							}*/

									}
								}
							}
						}
						else
						{

							//#pragma omp paralle for collapse (3)
							for( int k_small = k_ini_small; k_small<=k_end_small; ++k_small )
							{
								long long curr_node_z_outside_small = k_small*_numColumns*_numRows;
								for( int i_small = i_ini_small; i_small<=i_end_small; ++i_small )
								{
									long long curr_node_yz_outside_small = (i_small*_numColumns)+curr_node_z_outside_small;
									for( int j_small = j_ini_small; j_small<=j_end_small; ++j_small )
									{


										long long curr_node_xyz_outside_small = j_small + curr_node_yz_outside_small;

										inputPixelType value = _inputImageArray[curr_node_xyz_outside_small];
										double F_H;
										double B_H;

										// Not Sure if is P_! or P_I2
										if(floor(((double)value/factorBins))>=_alpha_3)
											F_H = (1-_P_I-_P_I2)*computePoissonProb((int)_alpha_3,_alpha_3);
										else
											// 								F_H = (1-_P_I)*computePoissonProb(value,_alpha_2);
											F_H = (1-_P_I-_P_I2)*computePoissonProb(floor(((double)value/factorBins)),_alpha_3);
										if(floor(((double)value/factorBins))<=_alpha_2)
											B_H = _P_I*computePoissonProb(int(_alpha_2),_alpha_2)+_P_I2*computePoissonProb(int(_alpha_1),_alpha_1);
										else
											// 								B_H = _P_I*computePoissonProb(value,_alpha_1);
											B_H = _P_I2*computePoissonProb(floor((double)value/factorBins),_alpha_2)+_P_I*computePoissonProb(floor((double)value/factorBins),_alpha_1);


										// 				std::cout << std::endl << F_H << " " << B_H;

										_costImageBArray[curr_node_xyz_outside_small] = -log(B_H);
										// 							if(_costImageBArray[curr_node_xyz_outside_small]>1000.0)
										// 								_costImageBArray[curr_node_xyz_outside_small] = 1000;
										// 							if(_costImageBArray[curr_node_xyz_outside_small]<0)
										// 								_costImageBArray[curr_node_xyz_outside_small] = 0;
										// 							
										_costImageFArray[curr_node_xyz_outside_small] = -log(F_H);
										// 							if(_costImageFArray[curr_node_xyz_outside_small]>1000.0)
										// 								_costImageFArray[curr_node_xyz_outside_small] = 1000;
										// 							if(_costImageFArray[curr_node_xyz_outside_small]<0)
										// 								_costImageFArray[curr_node_xyz_outside_small] = 0;

										// Just to test that each pixels is only access one time
										// 							_costImageABArray[curr_node_xyz_outside_small] = _costImageABArray[curr_node_xyz_outside_small] + 1;



										// // 							if( (k_small == 9) && (i_small == 45) && (j_small == 40) )
										// 							if( (int)value == 500 )
										// 							if( _costImageBArray[curr_node_xyz_outside_small] > _costImageFArray[curr_node_xyz_outside_small])
										// 							if( ((k_small == 162) && (i_small == 187) && (j_small == 128)) )// || ((k_small == 113) && (i_small == 228) && (j_small == 214)) || ((k_small == 113) && (i_small == 227) && (j_small == 227)) || ((k_small == 113) && (i_small == 238) && (j_small == 222)) )
										// 							{
										// 							    std::cout << std::endl << k << " " << i << " " << j;
										// 								std::cout << std::endl << "value: " << (int)value << ", " << _alpha_1 << ", " << _alpha_2 << ", " << _alpha_3 << ", " << _P_I << ", " << _P_I2 << ", " << B_H << ", " << F_H;
										// 								std::cout << std::endl << _costImageBArray[curr_node_xyz_outside_small] << ", " << _costImageFArray[curr_node_xyz_outside_small] << " " << min_J;
										// 							}
									}
								}
							}
						}
					}
				}
			}

			// // 	int nic;
			// // 	std::cin >>nic;
			// 
			// 	std::cout << std::endl << " Now lets check if all the pixels were visited";
			// 	std::cout << std::endl << " Now lets check if all the pixels were visited";
			// 	
			// 	
			// 	for( int k=0; k<_numStacks; ++k )
			// 	{
			// // 		std::cout << std::endl << "Stack: " << k << " of " << _numStacks << " Checked";
			// 		for( int i=0; i<_numRows; ++i )
			// 		{
			// 			//#pragma omp parallel for
			// 			for(int j=0; j<_numColumns; ++j )
			// 			{	
			// // 				std::cout << std::endl << "Stack: " << k << " of " << _numStacks;
			// // 				std::cout << std::endl << "	Row: " << i << " of " << _numRows;
			// // 				std::cout << std::endl << "	Col: " << j << " of " << _numColumns;
			// 				
			// 				
			// 				long long curr_node_z_outside = k*_numColumns*_numRows;
			// 				long long curr_node_yz_outside = (i*_numColumns)+curr_node_z_outside;
			// 				long long curr_node_xyz_outside = j + curr_node_yz_outside;
			// 				
			// 				if( _costImageABArray[curr_node_xyz_outside] != 1 )
			// 				  std::cout << std::endl << "BIG ERROR: " << _costImageABArray[curr_node_xyz_outside] << " " << k << " " << i << " " << j;
			// 			}
			// 		}
			// 	}


			// 	for( unsigned int i=0; i<_totNumPixels; ++i )
			// 	{
			// 		binIndex = floor(((double)_inputImageArray[i]/factorBins));	// !!!!!! DIVIDE THE IMAGE INTO TILES
			// 		if( binIndex >= _numberBins_mixPoisson )
			// 		{
			// 			binIndex = _numberBins_mixPoisson-1;
			// 		}
			// 		histoGram[binIndex]++;
			// 	}




			// 	// !!!!!! THIS PART IS STILL NOT INTEGRATED WITH THE DIVISION OF THE IMAGE TO CALCULATE THE COST B AND H IAMGES
			// 	
			// 	//Before entering the loop, compute the poisson probs 
			// 	for(int i=0; i<=_maxValueInputPixelType; ++i)
			// 	{
			// 		if(i>=_alpha_1)
			// 			F_H[i] = (1-_P_I)*computePoissonProb((int)_alpha_1,_alpha_1);
			// 		else
			// 			F_H[i] = (1-_P_I)*computePoissonProb(i,_alpha_1);
			// 		if(i<=_alpha_2)
			// 			B_H[i] = _P_I*computePoissonProb(int(_alpha_2),_alpha_2);
			// 		else
			// 			B_H[i] = _P_I*computePoissonProb(i,_alpha_2);
			// 	}







			// AFTER THIS WE SHOULD HAVE TWO IMAGES EACH ONE WITH THE COST OF ONE OR THE OTHER SEGMENTATION


			// 	int num_nodes;
			// 	int num_edges;
			// 	
			// 	long long curr_node_z;
			// 	long long curr_node_yz;
			// 	
			// 	long long curr_node_xy;
			// 	long long curr_node_xyz;
			// 	
			// 	long long rght_node;
			// 	long long down_node;
			// 	long long diag_node;
			// 	
			// 	long long diag_node_1;
			// 	long long rght_node_1;
			// 	long long down_node_1;
			// 	long long diag_node_2;
			// 	
			// 	
			// 	double Df;
			// 	double Db;
			// 	double Dr;
			// 	double Dd; 
			// 	double Dg; 
			// 	double F_H[_maxValueInputPixelType+1];
			// 	double B_H[_maxValueInputPixelType+1];







			std::cout << " Now lets set s-t costs";
			std::cout << std::endl;



			long long curr_node_z;
			long long curr_node_yz;

			long long curr_node_xy;
			long long curr_node_xyz;

			double Df;
			double Db;


			for( int k=0; k<_numStacks; ++k )
			{
				curr_node_z = k*_numColumns*_numRows;
				for(int i=0; i<_numRows; ++i)
				{
					curr_node_yz = (i*_numColumns)+curr_node_z;
					for(int j=0; j<_numColumns; ++j)
					{	
						curr_node_xyz = j + curr_node_yz;


						// !!!!!!!!!!!!!!!!!! THIS HAS TO BE CHANGED
						// 					Df = -log(F_H[intenPixel]);
						Df = _costImageFArray[curr_node_xyz];
						if(Df>1000.0)
							Df = 1000;
						// 					!!!!!!!!!!!!!!!!!! THIS HAS TO BE CHANGED
						// 					Db = -log(B_H[intenPixel]);
						Db = _costImageBArray[curr_node_xyz];
						if(Db>1000.0)
							Db=1000;
						// 				} 

						// 				if( (Df!=1000) && (Db!=1000) )
						// 				{
						// 					int n;

						// 					std::cout << std::endl << _alpha_1 << " " << _alpha_2 << " " << _P_I;
						// 					std::cout << std::endl << "DF: " << Df;
						// 					std::cout << std::endl << "DB: " << Db;

						// 					std::cin>>n;
						// 				}

						g -> add_node();
						g -> add_tweights( curr_node_xyz,   /* capacities */ Df,Db);  

						// 				if(( (k == 113) && (i == 228) && (j == 215)) || ((k == 113) && (i == 238) && (j == 222)) )
						// 					std::cout << std::endl << Df << " " << Db;

					}
				}
			}


			std::cout << " Now lets set neigh costs";
			std::cout << std::endl;

			long long rght_node;
			long long down_node;
			long long diag_node;

			long long diag_node_1;
			long long rght_node_1;
			long long down_node_1;
			long long diag_node_2;

			double Dr;
			double Dd; 
			double Dg; 

			double sigmaPow2 = 2*pow((double)_sigmaNeighCost,2);

			for( int k=0; k<_numStacks-1; ++k )
			{
				curr_node_z = k*_numColumns*_numRows;
				for(int i=0; i<_numRows-1; ++i)
				{
					curr_node_yz = (i*_numColumns)+curr_node_z;
					for(int j=0; j<_numColumns-1; ++j)
					{	
						// 				curr_node_xyz = j + curr_node_yz+1; // DUDA
						curr_node_xyz = j + curr_node_yz;

						rght_node = curr_node_xyz+1;
						down_node = curr_node_xyz+_numColumns;
						diag_node = curr_node_xyz+_numColumns*_numRows;



						Dr = _wNeigh*exp(-pow((double)_inputImageArray[curr_node_xyz]-(double)_inputImageArray[rght_node],2)/sigmaPow2);

						// 				std::cout << std::endl << curr_node_xyz << " " << rght_node << " " << down_node << " " << diag_node << " " << Dr;
						// 				std::cout << std::endl << curr_node_xyz << " " << rght_node << " " << down_node << " " << diag_node << " " << Dr;

						g->add_edge( curr_node_xyz, rght_node,    /* capacities */  Dr, Dr );	

						// 				std::cout << std::endl << curr_node_xyz << " " << rght_node << " " << down_node << " " << diag_node;
						// 				std::cout << std::endl << curr_node_xyz << " " << rght_node << " " << down_node << " " << diag_node;

						Dd = _wNeigh*exp(-pow((double)_inputImageArray[curr_node_xyz]-(double)_inputImageArray[down_node],2)/sigmaPow2);
						g->add_edge( curr_node_xyz, down_node,    /* capacities */  Dd, Dd );

						Dg = _wNeigh*exp(-pow((double)_inputImageArray[curr_node_xyz]-(double)_inputImageArray[diag_node],2)/sigmaPow2);
						g->add_edge( curr_node_xyz, diag_node,    /* capacities */  Dg, Dg );  

						// 				std::cout<< std::endl << Dr;


						// // 				In case of using 26 neighborhood system
						// 				diag_node_1 = down_node+1;
						// 				rght_node_1 = rght_node+_numColumns*_numRows;
						// 				down_node_1 = down_node+_numColumns*_numRows;
						// 				diag_node_2 = diag_node_1+_numColumns*_numRows;
						// 				
						// 				Dr = _wNeigh*exp(-pow((double)_inputImageArray[curr_node_xyz]-(double)_inputImageArray[diag_node_1],2)/sigmaPow2);
						// 				g->add_edge( curr_node_xyz, diag_node_1,    /* capacities */  Dr, Dr );	
						// 				
						// 				Dr = _wNeigh*exp(-pow((double)_inputImageArray[curr_node_xyz]-(double)_inputImageArray[rght_node_1],2)/sigmaPow2);
						// 				g->add_edge( curr_node_xyz, rght_node_1,    /* capacities */  Dr, Dr );	
						// 				
						// 				Dr = _wNeigh*exp(-pow((double)_inputImageArray[curr_node_xyz]-(double)_inputImageArray[down_node_1],2)/sigmaPow2);
						// 				g->add_edge( curr_node_xyz, down_node_1,    /* capacities */  Dr, Dr );	
						// // 				
						// 				Dr = _wNeigh*exp(-pow((double)_inputImageArray[curr_node_xyz]-(double)_inputImageArray[diag_node_2],2)/sigmaPow2);
						// 				g->add_edge( curr_node_xyz, diag_node_2,    /* capacities */  Dr, Dr );	
					}
				}
			} 

			// 	std::cout << std::endl << " Now lets compute maxflow";
			std::cout << " Now lets compute maxflow";
			std::cout << std::endl;


			//Compute the maximum flow:
			g->maxflow();



			//std::cin >> nic;

			// 	#pragma omp parallel for collapse(3)
			for( int k=0; k<_numStacks; ++k )
			{
				for(int i=0; i<_numRows; ++i)
				{
					for(int j=0; j<_numColumns; ++j)
					{	
						int curr_node_z2 = k*_numColumns*_numRows;
						int curr_node_yz2 = (i*_numColumns)+curr_node_z2;
						int curr_node_xyz2 = j + curr_node_yz2;

						if(g->what_segment(curr_node_xyz2) == GraphType::SOURCE)
						{
							_binaryImageArray[curr_node_xyz2] = 0;
						}
						else
						{
							_binaryImageArray[curr_node_xyz2] = _maxValueBinaryPixelType;
						}
					}
				}
			}
			std::cout << std::endl << "out";

			/*





			int SS,RR,CC;
			for(int i=0; i<num_nodes; ++i)
			{
			SS = ((long)i)%(_numColumns*_numRows);
			RR = (i-SS*_numColumns*_numRows)%_numColumns;
			CC = (i-SS*_numColumns*_numRows-);
			if(g->what_segment(i) == GraphType::SOURCE)
			{
			_binaryImageArray[RR*_numColumns + CC] = 0;
			}
			else
			{
			_binaryImageArray[RR*_numColumns + CC] = _maxValueBinaryPixelType;
			}
			}*/
			// 	std::cin >> nic;
			delete g;


			//typedef  itk::ImageFileWriter< binaryImageType  > WriterType;
			//WriterType::Pointer writer = WriterType::New();
			//writer->SetFileName("aout.tif");
			//writer->SetInput(_binaryImage);
			//writer->Update();

			//std::stringstream testNum;
			//testNum << row << "_" << col;
			//// 	testSubsBackNorm.tif
			//std::string filename2 = testNum.str() + ".tif";

			//typedef  itk::ImageFileWriter< inputImageType  > WriterType2;
			//WriterType2::Pointer writer2 = WriterType2::New();
			//writer2->SetFileName(filename2.c_str());
			//// 	writer->SetFileName("/data/nicolas/dapi_11_10/normNoBackground/kt11315_w410DAPIdsu_norm_back_crop_out.tif");
			//writer2->SetInput(_inputImage);
			//writer2->Update();



			//itk::Size<3> _binary_size = _binaryImage->GetLargestPossibleRegion().GetSize();
			//int _binary_slice_size_1 = _binary_size[1] * _binary_size[0];
			//int _binary_row_size_1 = _binary_size[0];
			//int x_offset = col*450;
			//int y_offset = row*450;

			//for(int z=0; z<_binary_size[2]; ++z)
			//{
			//	for(int y=0; y<_binary_size[1]; ++y)
			//	{					
			//		for(int x=0; x<_binary_size[0]; ++x)
			//		{
			//			unsigned char value = binArray[(slice_size*z) + (row_size*(y_offset + y)) + (x+x_offset)];
			//			unsigned char value_1 = _binaryImageArray[(_binary_slice_size_1*z) + (_binary_row_size_1*y) + (x)];
			//			if( (value == _maxValueBinaryPixelType) || (value_1 == _maxValueBinaryPixelType))
			//				binArray[(slice_size*z) + (row_size*(y_offset + y)) + (x+x_offset)] = _maxValueBinaryPixelType;			
			//		}
			//	}
			//}			
	//	}
	//}

	//WriterType::Pointer writer = WriterType::New();
	//writer->SetFileName("C:\\DATA\\TEST_REGISTRATION_16_BIT\\montage_kt11306_w410DAPIdsu_BS_binary.mhd");
	//writer->SetInput(binMontage);
	//writer->Update();


	// WRITE THE BINARY IMAGE






	return _binaryImage;
}

//2D Code ******************************************************
#define WinSz 256	//Histogram computed on this window
#define CWin  32	//This is half the inner window and must be a divisor of WinSz
#define NumBins 512	//Downsampled to these number of bins
#define HistPCMin 20.0  //Min percentage of histogram to search while updating poisson pars

typedef itk::Image< unsigned short, 2 > US2ImageType;
typedef itk::Image< unsigned short, 3 > US3ImageType;
typedef itk::Image< unsigned char, 3 > UC3ImageType;
typedef itk::Image< double, 2 > CostImageType;

template<typename InputImageType> typename InputImageType::Pointer
  CreateDefaultCoordsNAllocateSpace( typename InputImageType::SizeType size )
{
  typename InputImageType::Pointer inputImagePointer = InputImageType::New();
  typename InputImageType::PointType origin;
  typename InputImageType::IndexType start;
  const int imDims = InputImageType::ImageDimension;
  for( itk::IndexValueType i=0;
       i<imDims; ++i )
  {
    origin[i] = 0; start[i] = 0;
  }
  typename InputImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  inputImagePointer->SetOrigin( origin );
  inputImagePointer->SetRegions( region );
  inputImagePointer->Allocate();
  inputImagePointer->FillBuffer(0);
  try
  {
    inputImagePointer->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception caught in allocating space for image!" << excep << std::endl;
    exit (EXIT_FAILURE);
  }
  return inputImagePointer;
}

template<typename InputImageType, typename OutputImageType> void
RescaleCastNWriteImage
( typename InputImageType::Pointer inputImage,
    std::string &outFileName )
{
  typedef itk::RescaleIntensityImageFilter< InputImageType, OutputImageType > RescaleFilterType;
  if( InputImageType::ImageDimension!=OutputImageType::ImageDimension )
  {
    std::cout<<"This function needs equal input and output dimensions";
    return;
  }
  typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(inputImage);
  rescaleFilter->SetOutputMinimum( itk::NumericTraits< typename OutputImageType::PixelType >::min() );
  rescaleFilter->SetOutputMaximum( itk::NumericTraits< typename OutputImageType::PixelType >::max() );
  WriteITKImage< OutputImageType >( rescaleFilter->GetOutput(), outFileName );
  return;
}

template<typename InputImageType> void WriteITKImage
  ( typename InputImageType::Pointer inputImagePointer,
    std::string outputName )
{
  typedef typename itk::ImageFileWriter< InputImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputName.c_str() );
  writer->SetInput( inputImagePointer );
  try
  {
    writer->Update();
  }
  catch(itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
    exit( EXIT_FAILURE );
  }
  return;
}

void ComputeHistogram(
	US2ImageType::Pointer medFiltImages,
	std::vector< double > &histogram,
	US2ImageType::IndexType &start, US2ImageType::PixelType valsPerBin )
{
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType;
  US2ImageType::SizeType size; size[0] = WinSz; size[1] = WinSz;
  if( start[0]+WinSz >
      medFiltImages->GetLargestPossibleRegion().GetSize()[0] )
    size[0] = medFiltImages->GetLargestPossibleRegion().GetSize()[0]-start[0];
  if( start[1]+WinSz >
      medFiltImages->GetLargestPossibleRegion().GetSize()[1] )
    size[1] = medFiltImages->GetLargestPossibleRegion().GetSize()[1]-start[1];
  US2ImageType::RegionType region;
  region.SetSize( size ); region.SetIndex( start );
  ConstIterType constIter ( medFiltImages, region );
  for( constIter.GoToBegin(); !constIter.IsAtEnd(); ++constIter )
  {
    itk::SizeValueType currentBin = (itk::SizeValueType) 
		std::floor((double)constIter.Get()/(double)valsPerBin);
    if( histogram.size() <= currentBin ) currentBin = histogram.size()-1;
    ++histogram[currentBin];
  }
  double normalizeFactor = ((double)WinSz)*((double)WinSz);
  for( itk::SizeValueType j=0; j<histogram.size(); ++j )
  {
    histogram.at(j) /= normalizeFactor;
  }
  return;
}

void UpdateHistogram(
	US2ImageType::Pointer medFiltImages, std::vector< double > &histogram,
	US2ImageType::IndexType &start, US2ImageType::IndexType &prevStart,
	US3ImageType::PixelType valsPerBin )
{
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType;
  double pixelContrib = 1.00/(((double)WinSz)*((double)WinSz));
  //Remove the previous col from the histogram
  US2ImageType::IndexType startRem, startAdd;
  startRem[0] = prevStart[0]; startRem[1] = prevStart[1];
  US2ImageType::SizeType size;
  size[0] = start[0]-prevStart[0]; size[1] = start[1]-prevStart[1];
  US2ImageType::RegionType RemRegion, AddRegion;
  RemRegion.SetSize( size ); RemRegion.SetIndex( startRem );
  ConstIterType constIterRem ( medFiltImages, RemRegion );
  for( constIterRem.GoToBegin(); !constIterRem.IsAtEnd(); ++constIterRem )
  {
    itk::SizeValueType currentBin = (itk::SizeValueType) 
		std::floor((double)constIterRem.Get()/(double)valsPerBin);
    if( histogram.size() <= currentBin ) currentBin = histogram.size()-1;
    histogram[currentBin] -= pixelContrib;
  }
  startAdd[0] = start[0]; startAdd[1] = start[1]+WinSz-size[1];
  AddRegion.SetSize( size ); AddRegion.SetIndex( startAdd );
  ConstIterType constIterAdd ( medFiltImages, AddRegion );
  for( constIterAdd.GoToBegin(); !constIterAdd.IsAtEnd(); ++constIterAdd )
  {
    itk::SizeValueType currentBin = (itk::SizeValueType) 
		std::floor((double)constIterAdd.Get()/(double)valsPerBin);
    if( histogram.size() <= currentBin ) currentBin = histogram.size()-1;
    histogram[currentBin] += pixelContrib;
  }
  for( itk::SizeValueType i=0; i<histogram.size(); ++i ) 
    if( histogram.at(i) < pixelContrib ) histogram.at(i) = 0;
  return;
}

void computePoissonParams( std::vector< double > &histogram,
			   std::vector< double > &parameters, bool firstPass
			   )
{
  itk::SizeValueType max = histogram.size()-1;
  //The three-level min error thresholding algorithm
  double min_J = std::numeric_limits<double>::max();
  double P0, U0, P1, U1, U, J;
  // Try this: we need to define a penalty term that depends on the number of parameters
  //The penalty term is given as 0.5*k*ln(n)
  //where k is the number of parameters of the model and n is the number of samples
  double pc5=0, histPc=0;
  for( itk::SizeValueType i=0; i<histogram.size(); ++i )
  {
    histPc += histogram.at(i);
    if( histPc > 0.05 )
    {
      pc5 = i;
      break;
    }
  }
  if( pc5<1 ) pc5 = 1;

  itk::SizeValueType min_i, max_i;
  if( firstPass )
  {
    min_i = pc5;
    max_i = max;
  }
  else
  {
    double histPcChange = 2.0/WinSz*100;
    if( histPcChange<HistPCMin ) histPcChange=HistPCMin;
    histPcChange = std::ceil(histPcChange/100.0*((double)histogram.size()));
    min_i = (parameters.at(0)-histPcChange)<pc5 ? pc5 : (parameters.at(0)-histPcChange);
    max_i = (parameters.at(0)+histPcChange)>max ?
		max : (parameters.at(0)+histPcChange);
    if( parameters.at(2)>0.90 )//The window was in the background, estimate may be skewed
    {
      min_i = pc5;
      max_i = max;
    }
  }

  //Two level
  //The penalty term is given as sqrt(k)*ln(n)
  //In this case, k=4 and n=histogram.size()
  double PenTerm2 = 4.0*log(((double)(max+1)));
  for( itk::SizeValueType i=min_i; i<max_i; ++i )//to set the first threshold
  {
    //compute the current parameters of the first component
    P0 = U0 = 0.0;
    for( itk::SizeValueType l=0; l<=i; ++l )
    {
      P0 += histogram.at(l);
      U0 += (l+1)*histogram.at(l);
    }
    U0 /= P0;

    for( itk::SizeValueType j=i+1; j<max; ++j )//to set the second threshold
    {
      //compute the current parameters of the second component
      U1 = 0.0;
      P1 = 1-P0;
      for( itk::SizeValueType l=j; l<=max; ++l )
        U1 += (l+1)*histogram.at(l);
      U1 /= P1;

      //compute the overall mean
      U = P0*U0 + P1*U1;

      //Compute the current value of the error criterion function
      J = U - (P0*(log(P0)+U0*log(U0))+ P1*(log(P1)+U1*log(U1)));
      //Add the penalty term
      J += PenTerm2;

      if( J<min_J )
      {
        min_J = J;
        parameters.at(0) = U0; //Lowest mean
        parameters.at(1) = U1; //Intermediate mean
        parameters.at(2) = P0; //Prior for the lowest
        parameters.at(3) = P1; //Prior for the highest
      }
    }
  }
  return;
}

//Intialize pdf vector with max_intensity+1, 1
double ComputePoissonProbability(double intensity, double alpha)
{
	/*this is the equation
	P = (alpha^intensity)*exp(-alpha)/factorial(intensity);*/
	double A, P;
	A = std::exp(-alpha);
	P = 1;
	for (US2ImageType::PixelType i=1; i<= intensity; ++i)
		P = P * (alpha/i);
	P = P*A;
	if(P < std::numeric_limits<double>::epsilon())
		P = std::numeric_limits<double>::epsilon();
	return P;
}

void ComputeCosts( US2ImageType::Pointer  medFiltImages,
		   CostImageType::Pointer flourCosts,
		   CostImageType::Pointer flourCostsBG,
		   US2ImageType::PixelType valsPerBin
		   )
{
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType;
  typedef itk::ImageRegionIterator< CostImageType > CostIterType;
  typedef itk::ImageRegionIterator< US3ImageType > IterTypeUS3d;
  itk::IndexValueType numCol =  medFiltImages->
				GetLargestPossibleRegion().GetSize()[1];
  itk::IndexValueType numRow =  medFiltImages->
				GetLargestPossibleRegion().GetSize()[0];
  itk::IndexValueType WinSz2 = (itk::IndexValueType)floor(((double)WinSz)/2+0.5)
			      -(itk::IndexValueType)floor(((double)CWin)/2+0.5);
  //WinSz2*2+CWin defines the whole window
  itk::IndexValueType WinSzHalf = (itk::IndexValueType)floor(((double)WinSz)/2+0.5);

  for( itk::IndexValueType i=0; i<numRow; i+=CWin )
  { 
    std::vector< double > parameters( 4, 0 );
    //std::vector< double > histogram( NumBins, 0 );
    US2ImageType::IndexType prevStart; prevStart[0] = 0; prevStart[1] = 0;
    for( itk::IndexValueType j=0; j<numCol; j+=CWin )
    {
    std::vector< double > histogram( NumBins, 0 );
      //Compute histogram at point i,j with window size define WinSz
      US2ImageType::IndexType curPoint; curPoint[0] = i; curPoint[1] = j;
      US2ImageType::IndexType start; start[0] = i-WinSz2; start[1] = j-WinSz2;
      if( (start[0]+WinSz)>=numRow ) start[0] =  numRow-WinSz-1;
      if( (start[1]+WinSz)>=numCol ) start[1] =  numCol-WinSz-1;
      if( start[0]<0 ) start[0] = 0; if( start[1]<0 ) start[1] = 0;
      if( j==0 )
      {
        ComputeHistogram( medFiltImages, histogram, start, valsPerBin );
	computePoissonParams( histogram, parameters, true );
      }
      else if( (start[1]!=prevStart[1]) || (start[0]!=prevStart[0]) )
      {
        ComputeHistogram( medFiltImages, histogram, start, valsPerBin );
	//UpdateHistogram( medFiltImages, histogram, start, prevStart, valsPerBin );
	computePoissonParams( histogram, parameters, false );
      }
      prevStart[0] = start[0]; prevStart[1] = start[1];

	//Declare iterators for the two cost images
	US2ImageType::SizeType size; size[0] = CWin; size[1] = CWin;
	if( (i+CWin)>=numRow ) size[0] = numRow-i-1;
	if( (j+CWin)>=numCol ) size[1] = numCol-j-1;
	US2ImageType::RegionType region;
	region.SetSize( size ); region.SetIndex( curPoint );
	ConstIterType constIter ( medFiltImages, region );
	CostIterType costIterFlour	( flourCosts,	region );
	CostIterType costIterFlourBG	( flourCostsBG,	region );
	constIter.GoToBegin(); costIterFlour.GoToBegin(); costIterFlourBG.GoToBegin();
	for( ; !constIter.IsAtEnd(); ++constIter, ++costIterFlour, ++costIterFlourBG )
	{
	  itk::SizeValueType currentPixel = (US2ImageType::PixelType)std::floor
						( ((double)constIter.Get())/((double)valsPerBin) );

	  //Compute node costs for each type
	  double F, FBG;
	  if( currentPixel < 1 ) currentPixel = 1;
	  if( (parameters.at(1)-parameters.at(0))<1 )//|| currentPixel < 1 )
	  {
	    FBG = 10000.0;
	    F   = 0;
	  }
	  else
	  {
	    if( currentPixel >= parameters.at(1) )
	      F  =  parameters.at(3)*ComputePoissonProbability(parameters.at(1),
	    							parameters.at(1));
	    else
	      F  =  parameters.at(3)*ComputePoissonProbability(currentPixel,
	    							parameters.at(1));
	    if( currentPixel <= parameters.at(0) )
	      FBG = parameters.at(2)*ComputePoissonProbability(parameters.at(0),
	    							parameters.at(0));
	    else
	      FBG = parameters.at(2)*ComputePoissonProbability(currentPixel,
	    							parameters.at(0));
	      F    = -log( F    ); if( F    > 10000.0 ) F    = 10000.0;
	      FBG  = -log( FBG  ); if( FBG  > 10000.0 ) FBG  = 10000.0;
	  }
	  /*CostImageType::IndexType curInd = costIterFlour.GetIndex();
	  if( curInd[1] == 0  && curInd[0]>1120 && curInd[0]<1128 )
	  {
	    std::cout<<"FGC="<<F<<" BGC="<<FBG<<" MuB="<<parameters.at(0)
	    	<<" MuF="<<parameters.at(1)<<" PB="<<parameters.at(2)
		<<" PF="<<parameters.at(3)<<std::endl;
	  }*/
	  costIterFlour.Set( F ); costIterFlourBG.Set( FBG );
        }
    }
  }
  return;
}

void ComputeCut( US2ImageType::Pointer  medFiltImages,
		 CostImageType::Pointer flourCosts,
		 CostImageType::Pointer flourCostsBG,
		 UC3ImageType::Pointer outputImage,
		 UC3ImageType::PixelType foregroundValue
		)
{
  double sigma = 25.0; //What! A hard coded constant check Boykov's paper!! Also check 20 in weights
  double neighWt = 80.0;
  typedef itk::ImageRegionIteratorWithIndex< CostImageType > CostIterType;
  typedef itk::ImageRegionIteratorWithIndex< US2ImageType > US2IterType;
  typedef itk::ImageRegionIteratorWithIndex< UC3ImageType > UC3IterType;
  //Compute the number of nodes and edges
  itk::SizeValueType numRow   = medFiltImages->GetLargestPossibleRegion().GetSize()[0];
  itk::SizeValueType numCol   = medFiltImages->GetLargestPossibleRegion().GetSize()[1];
  itk::SizeValueType numNodes = numCol*numRow;
  itk::SizeValueType numEdges = 3*numCol*numRow //Down, right and diagonal
  				+ 1 - 2*numCol	//No Down At Bottom
				- 2*numRow;	//No Right At Edge
  //Declare some common iterators
  US2IterType medianIter( medFiltImages,
  			  medFiltImages->GetLargestPossibleRegion() );
  CostIterType AFCostIter( flourCosts, 
  			   flourCosts->GetLargestPossibleRegion() );
  CostIterType AFBGCostIter( flourCostsBG, 
  			     flourCostsBG->GetLargestPossibleRegion() );
  UC3ImageType::IndexType start;start[0] = 0;	  start[1] = 0;	    start[2] = 0;
  UC3ImageType::SizeType size;   size[0] = numRow; size[1] = numCol; size[2] = 1;
  UC3ImageType::RegionType region; region.SetSize( size ); region.SetIndex( start );
  UC3IterType outputIter( outputImage, region );
  typedef Graph_B < double, double, double > GraphType;
  GraphType *graph = new GraphType( numNodes, numEdges );
  //Iterate and add terminal weights
  for( itk::SizeValueType i=0; i<numRow; ++i )
  {
    for( itk::SizeValueType j=0; j<numCol; ++j )
    {
      CostImageType::IndexType index; index[0] = i; index[1] = j;
      AFCostIter.SetIndex( index ); AFBGCostIter.SetIndex( index );
      itk::SizeValueType indexCurrentNode = i*numCol+j; 
      graph->add_node();
      graph->add_tweights( indexCurrentNode, //rand()%100, rand()%100 );//Testing
			AFCostIter.Get(), AFBGCostIter.Get() );
    }
  }
  for( itk::SizeValueType i=0; i<numRow-1; ++i )
  {
    for( itk::SizeValueType j=0; j<numCol-1; ++j )
    {
      US2ImageType::IndexType index; index[0] = i; index[1] = j; medianIter.SetIndex( index );
      double currentVal = medianIter.Get();
      //Intensity discontinuity terms as edges as done in Yousef's paper
      itk::SizeValueType indexCurrentNode  = i*numCol+j;
      itk::SizeValueType indexRightNode    = indexCurrentNode+1;
      itk::SizeValueType indexBelowNode    = indexCurrentNode+numCol;
      itk::SizeValueType indexDiagonalNode = indexBelowNode+1;
      //Right
      index[0] = i; index[1] = j+1; medianIter.SetIndex( index );
      double rightCost = neighWt*exp(-pow(currentVal-medianIter.Get(),2)/(2*pow(sigma,2)));
      graph->add_edge( indexCurrentNode, indexRightNode, rightCost, rightCost );
      //Below
      index[0] = i+1; index[1] = j; medianIter.SetIndex( index );
      double downCost = neighWt*exp(-pow(currentVal-medianIter.Get(),2)/(2*pow(sigma,2)));
      graph->add_edge( indexCurrentNode, indexBelowNode, downCost, downCost );
      //Diagonal
      index[0] = i+1; index[1] = j+1; medianIter.SetIndex( index );
      double diagonalCost = neighWt*exp(-pow(currentVal-medianIter.Get(),2)/(2*pow(sigma,2)));
      graph->add_edge( indexCurrentNode, indexDiagonalNode, diagonalCost, diagonalCost );
    }
  }
  //Max flow:
  graph->maxflow();

  //Iterate and write out the pixels
  for( itk::SizeValueType i=0; i<numRow-1; ++i )
  {
    for( itk::SizeValueType j=0; j<numCol-1; ++j )
    {
      UC3ImageType::IndexType index; index[0] = i; index[1] = j; index[2] = 0;
      outputIter.SetIndex( index );
      itk::SizeValueType indexCurrentNode = i*numCol+j;
      if( graph->what_segment( indexCurrentNode ) != GraphType::SOURCE )
	outputIter.Set( foregroundValue );
    }
  }
  delete graph;
}

US2ImageType::Pointer Get2DFrom3D( US3ImageType::Pointer readImage )
{
  typedef itk::ExtractImageFilter< US3ImageType, US2ImageType > DataExtractType;
  DataExtractType::Pointer deFilter = DataExtractType::New();
  US3ImageType::RegionType dRegion  = readImage->GetLargestPossibleRegion();
  dRegion.SetSize (2,0);
  dRegion.SetIndex(2,0);
  deFilter->SetExtractionRegion(dRegion);
  deFilter->SetDirectionCollapseToIdentity();
  deFilter->SetInput( readImage );
  try
  {
    deFilter->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception caught in slice extraction filter!" << excep << std::endl;
    exit (EXIT_FAILURE);
  }
  US2ImageType::Pointer d2Im = deFilter->GetOutput();
  return d2Im;
}

UC3ImageType::Pointer AdaptiveBinarization2D( US3ImageType::Pointer inputImage3d )
{
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType;

  US2ImageType::Pointer inputImage = Get2DFrom3D( inputImage3d );
  UC3ImageType::SizeType size3d;
  US2ImageType::SizeType size;
  size3d[0] = size[0] = inputImage->GetLargestPossibleRegion().GetSize()[0];
  size3d[1] = size[1] = inputImage->GetLargestPossibleRegion().GetSize()[1];
  size3d[2] = 1;

  double valsPerBinDbl = -1;

  //Get max pixel val
  ConstIterType medianIter( inputImage, inputImage->GetLargestPossibleRegion() );
  for( medianIter.GoToBegin(); !medianIter.IsAtEnd(); ++medianIter )
    if( medianIter.Get() > valsPerBinDbl )
      valsPerBinDbl = medianIter.Get();
  valsPerBinDbl /= (double)NumBins;
  valsPerBinDbl = std::floor( valsPerBinDbl+0.5 );

  CostImageType::Pointer flourCosts   =
			CreateDefaultCoordsNAllocateSpace<CostImageType>( size );
  CostImageType::Pointer flourCostsBG =
			CreateDefaultCoordsNAllocateSpace<CostImageType>( size );
  UC3ImageType::Pointer labelImage =
			CreateDefaultCoordsNAllocateSpace<UC3ImageType>( size3d );

  ComputeCosts( inputImage, flourCosts, flourCostsBG,
  					(US2ImageType::PixelType)valsPerBinDbl );
  ComputeCut( inputImage, flourCosts, flourCostsBG, labelImage,
			std::numeric_limits<UC3ImageType::PixelType>::max() );
/*
  std::string nameFg = "bfg.tif";
  std::string nameBg = "bbg.tif";
  std::string nameBn = "bin.tif";
  RescaleCastNWriteImage< CostImageType, US2ImageType>( flourCosts, nameFg );
  RescaleCastNWriteImage< CostImageType, US2ImageType>( flourCostsBG, nameBg );
  RescaleCastNWriteImage< UC3ImageType, UC3ImageType>( labelImage, nameBn );*/

  return labelImage;
}
