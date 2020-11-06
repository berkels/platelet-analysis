/**
 * \file
 * \brief Segments an image into multiple piecewise constant regions with the Mumford Shah model.
 *
 * Segments an image into multiple piecewise constant regions with the Mumford Shah model.
 * The minimizer is approximated by the primal / dual algorithm propsed
 * by Antonin Chambolle and Thomas Pock in
 * {A first-order primal-dual algorithm for convex problems with applications to imaging}.
 *
 * \author Mevenkamp
 */

#include <aol.h>
#include <parameterParser.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <vectorExtensions.h>
#include <solver.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <configurators.h>
#include <ChanVese.h>
#include <generator.h>
#include <segmentation.h>
#include <primalDualSegmentation.h>

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main( int argc, char **argv ) {

  try {
    aol::ParameterParser parser( argc, argv, "MultiPhaseMSSeg.par" );

    aol::StopWatch watch;
    watch.start();

    typedef RType RealType;
    typedef ConfType ConfiguratorType;

    const RealType gamma = parser.getDouble( "gamma" );

    qc::DefaultArraySaver<RealType, ConfiguratorType::Dim> saver;
    saver.initFromParser( parser, true );
    
    int imageDim = parser.getInt ( "imageDim" );
    
    if ( imageDim == 1 ) {
      qc::MultiArray<RealType, ConfiguratorType::Dim, 1> u0 ( parser.getString( "input-image" ).c_str() );
      if ( parser.hasVariable( "thresholdInputAt" ) )
        u0.threshold ( parser.getDouble ( "thresholdInputAt" ), 0, 1 );
      else
        u0 /= u0.getMaxValue();
      u0.setOverflowHandling( aol::CLIP_THEN_SCALE, 0, 1 );
      u0.savePNG( saver.createSaveName( "", ".png", -1, "input" ).c_str() );
      ConfiguratorType::InitType grid ( qc::GridSize<ConfType::Dim> ( u0.getSize() ) );
      
      const int numSegments = parser.getInt ( "numSegments" );

      im::PiecewiseConstantMultiPhaseMSSegmentor<ConfiguratorType, im::FirstOrderPrimalDualMultiPhaseMSSegmentor<ConfiguratorType> > segmentor ( grid, gamma, u0, true, false, "", numSegments );
      aol::VectorContainer<qc::ScalarArray<RealType, ConfiguratorType::Dim> > temp ( numSegments, qc::ScalarArray<RealType, ConfiguratorType::Dim> ( grid ) );

      segmentor.setMaxIterations ( parser.getInt( "numSteps") );
      segmentor.setStopEpsilon ( parser.getDouble ( "epsilon" ) );
      if ( parser.hasVariable( "outerIterations"  ) )
        segmentor.setOuterIterations( parser.getInt ( "outerIterations" ) );
      
      if ( parser.hasVariable( "initialGrayValues" ) ) {
        aol::MultiVector<RealType> initialGrayValues;
        parser.getRealMultiVec( "initialGrayValues", initialGrayValues );
        segmentor.getMeanValuesReference() = initialGrayValues;
        std::cerr << initialGrayValues << std::endl;
      }

      aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > pDual ( numSegments, qc::MultiArray<RealType, ConfiguratorType::Dim> ( grid ) );
      if ( parser.hasVariable ( "adjustGrayValues" ) && !parser.checkAndGetBool ( "adjustGrayValues" ) )
        segmentor.segment ( temp, &pDual );
      else
        segmentor.segmentAndAdjustGrayValues ( temp, &pDual );

      qc::ScalarArray<int, qc::QC_2D> hardSegmentation ( grid );
      segmentor.getHardSegmentation ( hardSegmentation, temp );
      hardSegmentation.setOverflowHandlingToCurrentValueRange ( );
      hardSegmentation.savePNG ( saver.createSaveName( "", ".png", -1, "labels" ).c_str() );
      
      qc::MultiArray<RealType, ConfiguratorType::Dim, 1> meanValuesImage ( grid );
      segmentor.getMeanValuesImage ( meanValuesImage, hardSegmentation );
      meanValuesImage.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
      meanValuesImage.savePNG ( saver.createSaveName( "", ".png", -1, "segmentation" ).c_str() );
    } else if ( imageDim == 3 ) {
      qc::MultiArray<RealType, ConfiguratorType::Dim, 3> u0 ( parser.getString( "input-image" ).c_str() );
      if ( parser.hasVariable( "thresholdInputAt" ) )
        u0.threshold ( parser.getDouble ( "thresholdInputAt" ), 0, 1 );
      else
        u0 /= u0.getMaxValue();
      u0.setOverflowHandling( aol::CLIP_THEN_SCALE, 0, 1 );
      u0.savePNG( saver.createSaveName( "", ".png", -1, "input" ).c_str() );
      ConfiguratorType::InitType grid ( qc::GridSize<ConfType::Dim> ( u0.getSize() ) );
      
      const int numSegments = parser.getInt ( "numSegments" );
      
      im::PiecewiseConstantMultiPhaseMSSegmentor<ConfiguratorType, im::FirstOrderPrimalDualMultiPhaseMSSegmentor<ConfiguratorType> > segmentor ( grid, gamma, u0, true, false, "", numSegments );
      aol::VectorContainer<qc::ScalarArray<RealType, ConfiguratorType::Dim> > temp ( numSegments, qc::ScalarArray<RealType, ConfiguratorType::Dim> ( grid ) );
      
      segmentor.setMaxIterations ( parser.getInt( "numSteps") );
      segmentor.setStopEpsilon ( parser.getDouble ( "epsilon" ) );
      if ( parser.hasVariable( "outerIterations"  ) )
        segmentor.setOuterIterations( parser.getInt ( "outerIterations" ) );
      
      if ( parser.hasVariable( "initialGrayValues" ) ) {
        aol::MultiVector<RealType> initialGrayValues;
        parser.getRealMultiVec( "initialGrayValues", initialGrayValues );
        segmentor.getMeanValuesReference() = initialGrayValues;
        std::cerr << initialGrayValues << std::endl;
      }
      
      aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > pDual ( numSegments, qc::MultiArray<RealType, ConfiguratorType::Dim> ( grid ) );
      if ( parser.hasVariable ( "adjustGrayValues" ) && !parser.checkAndGetBool ( "adjustGrayValues" ) )
        segmentor.segment ( temp, &pDual );
      else
        segmentor.segmentAndAdjustGrayValues ( temp, &pDual );
      
      qc::ScalarArray<int, qc::QC_2D> hardSegmentation ( grid );
      segmentor.getHardSegmentation ( hardSegmentation, temp );
      hardSegmentation.setOverflowHandlingToCurrentValueRange ( );
      hardSegmentation.savePNG ( saver.createSaveName( "", ".png", -1, "labels" ).c_str() );
      
      qc::MultiArray<RealType, ConfiguratorType::Dim, 3> meanValuesImage ( grid );
      segmentor.getMeanValuesImage ( meanValuesImage, hardSegmentation );
      meanValuesImage.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
      meanValuesImage.savePNG ( saver.createSaveName( "", ".png", -1, "segmentation" ).c_str() );
    }

    watch.stop();
    watch.printReport( cerr );
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
