#include <cellCenteredGrid.h>
#include <paramReg.h>
#include <mutualInformation.h>

namespace aol {
 
/**
 * Based on https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
 */
class CommandLineArgumentParser {
public:
  CommandLineArgumentParser ( const int argc, char **argv ) {
    for ( int i = 1; i < argc; ++i )
      this->tokens.push_back ( std::string ( argv[i] ) );
  }
  /// @author iain
  const std::string& getCmdOption ( const std::string &option ) const {
    std::vector<std::string>::const_iterator itr;
    itr =  std::find ( this->tokens.begin(), this->tokens.end(), option );
    if ( itr != this->tokens.end() && ++itr != this->tokens.end() ) {
      return *itr;
    }
    static const std::string empty_string ( "" );
    return empty_string;
  }
  /// @author iain
  bool cmdOptionExists ( const std::string &option ) const {
    return std::find ( this->tokens.begin(), this->tokens.end(), option ) != this->tokens.end();
  }
private:
  std::vector <std::string> tokens;
};

}

template <typename ConfiguratorType>
class ThrombAnalysis {
  typedef typename ConfiguratorType::RealType RealType;
  const qc::ScalarArray<typename ConfiguratorType::RealType, ConfiguratorType::Dim> &_dataImage;
  const RealType _maxIntensityValue;
  const typename ConfiguratorType::InitType _grid;
  qc::ScalarArray<int, ConfiguratorType::Dim> _labelArray;
  aol::Vector<int> _corePixels;
  aol::Vector<int> _tailPixels;
  int _numberOfLabels;
  aol::RandomAccessContainer<aol::Matrix22<int> > _boundingBox;
public:
  ThrombAnalysis ( const qc::ScalarArray<RealType, ConfiguratorType::Dim> &DataImage )
   : _dataImage ( DataImage ),
     _maxIntensityValue ( _dataImage.getMaxValue() ),
     _grid ( qc::GridSize<ConfiguratorType::Dim> ( _dataImage.getSize() ) ),
     _labelArray ( _grid ) {

    qc::BitArray<ConfiguratorType::Dim> mask ( _grid );
    mask.thresholdFrom ( _dataImage, _dataImage.getMinValue() );
    _numberOfLabels = qc::ConnectedComponentsLabeler::doLabel ( mask, _labelArray );

     cerr << "Found " << _numberOfLabels << " components." << endl;

    _corePixels.reallocate ( _numberOfLabels );
    _tailPixels.reallocate ( _numberOfLabels );
    _boundingBox.reallocate ( _numberOfLabels );
    for ( int i = 0; i < _numberOfLabels; ++i ) {
      _boundingBox[i][0][0] = _labelArray.getNumX();
      _boundingBox[i][0][1] = -1;
      _boundingBox[i][1][0] = _labelArray.getNumY();
      _boundingBox[i][1][1] = -1;
    }

    for ( int j = 0; j < _labelArray.getNumY(); ++j ) {
      for ( int i = 0; i < _labelArray.getNumX(); ++i ) {
        if ( _labelArray.get(i,j) > 0 ) {
          const int labelIndex = _labelArray.get(i,j) - 1;
          _boundingBox[labelIndex][0][0] = aol::Min ( _boundingBox[labelIndex][0][0], i );
          _boundingBox[labelIndex][0][1] = aol::Max ( _boundingBox[labelIndex][0][1], i );
          _boundingBox[labelIndex][1][0] = aol::Min ( _boundingBox[labelIndex][1][0], j );
          _boundingBox[labelIndex][1][1] = aol::Max ( _boundingBox[labelIndex][1][1], j );
          if ( aol::appeqAbsolute ( _dataImage.get(i,j), _maxIntensityValue ) )
            _corePixels[labelIndex]++;
          else
            _tailPixels[labelIndex]++;
        }
      }
    }
  }
  
  void run ( const bool saveComponentMasks = false ) const {
    aol::Vector<int> boundaryPixels ( _numberOfLabels );
    aol::Vector<int> brightComponents ( _numberOfLabels );
    aol::Vector<int> darkComponents ( _numberOfLabels );

    const int numX = _labelArray.getNumX();
    const int numY = _labelArray.getNumY();
    qc::BitArray<ConfiguratorType::Dim> boundaryMask ( _grid );
    for ( int j = 0; j < numY; ++j ) {
      for ( int i = 0; i < numX; ++i ) {
        const int label = _labelArray.get ( i, j );
        
        // Skip background pixels
        if ( label == 0 )
          continue;

        bool onBoundary = false;
        for ( int y = aol::Max ( j-1, 0 ); y <= aol::Min ( j+1, numY - 1 ); ++y ) {
          for ( int x = aol::Max ( i-1, 0 ); x <= aol::Min ( i+1, numX - 1 ); ++x ) {
            if ( label != _labelArray.get ( x, y ) ) {
              onBoundary = true;
              break;
            }
          }
          if ( onBoundary == true )
            break;
        }

        if ( onBoundary ) {
          boundaryPixels[ label - 1]++;
          boundaryMask.set ( i, j, true );
        }
      }
    }
    boundaryMask.save ( "boundary.pgm" );

    //aol::Vector<int> histo;
    //pixelVolume.createHistogramOfValues ( histo, 256 );
    //aol::plotHistogram<RType> ( histo, "histo", false, false );

    aol::RandomGenerator rndGen;
    aol::RandomAccessContainer<aol::Vec3<int> > colors ( _numberOfLabels + 1 );
    for ( int i = 1; i <= _numberOfLabels; ++i )
      for ( int j = 0; j < 3; ++j )
        colors[i][j] = rndGen.rUnsignedInt ( 256 );

    qc::MultiArray<unsigned char, ConfiguratorType::Dim, 3> bufArray ( _grid );
    for ( int i = 0; i < _labelArray.size(); ++i )
      for ( int j = 0; j < 3; ++j )
        bufArray[j][i] = colors[_labelArray[i]][j];

    bufArray.savePNG( "componentsCol.png" );

    bufArray.setZero();

    if ( saveComponentMasks && ( aol::directoryExists ( "components" ) == false ) )
      aol::makeDirectory ( "components" );

    for ( int i = 0; i < _numberOfLabels; ++i ) {
      const aol::Vec2<int> start = getLocalComponentBBStart(i);
      const aol::Vec2<int> stop = getLocalComponentBBStop(i);
      qc::BitArray<ConfiguratorType::Dim> localMask ( stop[0] - start[0] + 1, stop[1] - start[1] + 1 );
      qc::ScalarArray<int, ConfiguratorType::Dim> localLabelArray ( localMask.getSize() );

      getLocalComponentMask ( localMask, i, false, false );
      if ( saveComponentMasks )
        localMask.save ( aol::strprintf ( "components/componentBright-%03d.pgm", i ).c_str() );
      qc::cleanMask ( localMask, 2, false );
      brightComponents[i] = qc::ConnectedComponentsLabeler::doLabel ( localMask, localLabelArray );

      getLocalComponentMask ( localMask, i, false, true );
      if ( saveComponentMasks )
        localMask.save ( aol::strprintf ( "components/componentDark-%03d.pgm", i ).c_str() );
      qc::cleanMask ( localMask, 2, false );
      localLabelArray.setZero();
      darkComponents[i] = qc::ConnectedComponentsLabeler::doLabel ( localMask, localLabelArray );

      // Visualize the bounding boxes
      for ( int x = start[0]; x <= stop[0]; ++x )
        for ( int j = 0; j < 3; ++j ) {
          for ( int k = 0; k <= 1; ++k ) {
            const int y = _boundingBox[i][1][k];
            bufArray[j].set ( x, y, colors[i+1][j] );
          }
        }
    }
    bufArray.savePNG( "componentsBBCol.png" );

    _labelArray.saveTIFF ( "components.tif" );
    _labelArray.saveMetaImageFile ( "components", qc::PGM_UNSIGNED_SHORT_BINARY );

    {
      std::ofstream out ( "stat.csv" );
      if ( out.is_open() )
      {
        out << "corePixels;tailPixels;boundaryPixels;brightComponents;darkComponents\n";
        for ( int i = 0; i < _numberOfLabels; ++i ) {
          out << _corePixels[i] << ";" << _tailPixels[i] << ";" << boundaryPixels[i]
          << ";" << brightComponents[i] << ";" << darkComponents[i] << endl;
        }
        out.close();
      }
      else
        cerr << "Unable to open file \"stat.csv\" for writing." << endl;
    }
  }
  
  aol::Vec2<int> getLocalComponentBBStart ( const int CompNum ) const {
    return  aol::Vec2<int> ( _boundingBox[CompNum][0][0], _boundingBox[CompNum][1][0] );
  }

  aol::Vec2<int> getLocalComponentBBStop ( const int CompNum ) const {
    return  aol::Vec2<int> ( _boundingBox[CompNum][0][1], _boundingBox[CompNum][1][1] );
  }

  void getLocalComponentMask ( qc::BitArray<ConfiguratorType::Dim> &LocalMask,
                               const int CompNum,
                               const bool BackgroundValue,
                               const bool MaxValToFalse ) const {
    const aol::Vec2<int> start = getLocalComponentBBStart(CompNum);
    LocalMask.setZero();
    for ( int y = 0; y < LocalMask.getNumY(); ++y ) {
      for ( int x = 0; x < LocalMask.getNumX(); ++x ) {
        if ( _labelArray.get(x + start[0],y + start[1]) != CompNum + 1 )
          LocalMask.set ( x, y, BackgroundValue );
        else {
          const bool isMaxVal = aol::appeqAbsolute ( _dataImage.get(x + start[0],y + start[1]), _maxIntensityValue );
          LocalMask.set ( x, y, MaxValToFalse ? !isMaxVal : isMaxVal );
        }
      }
    }
  }
};

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main( int argc, char **argv ) {
  try {
    if ( argc < 2 ){
      cerr << "USAGE: " << argv[0] << " <input_file> [options]" << endl;
      return EXIT_FAILURE;
    }
    
    aol::CommandLineArgumentParser CLI ( argc, argv );

    cerr << "Loading file \"" << argv[1] << "\" ... ";
    qc::ScalarArray<RType, qc::QC_2D> input ( argv[1] );
    cerr << "done." << endl;
    ThrombAnalysis<ConfType> thrombAnalysis ( input );
    thrombAnalysis.run( CLI.cmdOptionExists ( "--saveComponentMasks" ) );
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
