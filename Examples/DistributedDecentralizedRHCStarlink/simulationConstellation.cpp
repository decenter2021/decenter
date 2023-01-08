// ---------- Headers ----------
// Standard headers
#include <ctime>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <fstream>

// Boost headers
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>

// Eigen header 
#include <Eigen/Core>
// Tudat headers
#include <tudat/simulation/simulation.h>
#include <tudat/astro/basic_astro.h>
#include <tudat/astro/gravitation.h>
#include <tudat/astro/propulsion.h>

// User defined libraries
#include "applicationOutput.h" // Minor changes in relation to the file available in tudat examples
#include "tudatThrustFeedbackMatlab.h"
#include "constellationParameters.h" 


// ---------- Main application ----------
int main(){

        // --------------------------------------
        // ---------- Using statements ----------
        // --------------------------------------
        using namespace tudat;
        using namespace tudat::simulation_setup;
        using namespace tudat::propagators;
        using namespace tudat::numerical_integrators;
        using namespace tudat::orbital_element_conversions;
        using namespace tudat::basic_mathematics;
        using namespace tudat::unit_conversions;
        using namespace tudat::interpolators;
        using namespace tudat::gravitation;
        using namespace tudat::numerical_integrators;

        // -----------------------------------------------------
        // ---------- Create environment and vehicle  ----------
        // -----------------------------------------------------

        // ---------- Define constellation ----------
        // Set number of satellites in constellation.
        
        const unsigned int numberOfPlanes = CONSTELLATION_N_PLANES;
        const unsigned int numberOfSatellitesPerPlane = CONSTELLATION_N_PER_PLANE;
        const unsigned int numberOfSatellites = numberOfPlanes*numberOfSatellitesPerPlane;

        // ---------- Set simulation epoch ----------
        const double simulationStartEpoch = EPOCH_START;
        const double simulationEndEpoch = EPOCH_END;

        // Load Spice kernels.
        spice_interface::loadStandardSpiceKernels( );
        
        // Create body object (array of strings)
        std::vector <std::string> bodiesToCreate;
        bodiesToCreate.push_back("Sun");
        bodiesToCreate.push_back("Earth");
        bodiesToCreate.push_back("Moon");
        bodiesToCreate.push_back("Mars");
        bodiesToCreate.push_back("Venus");
        bodiesToCreate.push_back("Jupiter");
        

        // Create settings for objects (key-value pairs)
        BodyListSettings bodySettings = getDefaultBodySettings(bodiesToCreate,simulationStartEpoch - ENVIRONMENT_TIME_BUFFER,simulationEndEpoch + ENVIRONMENT_TIME_BUFFER, "Earth", "J2000" );
        // Adjust settings (Ephemeris)

        for (unsigned int i = 0; i < bodiesToCreate.size( ); i++){
            bodySettings.at(bodiesToCreate.at(i))->ephemerisSettings->resetFrameOrientation("J2000");
            bodySettings.at(bodiesToCreate.at(i))->rotationModelSettings->resetOriginalFrame("J2000");
        }
        // Adjust atmosphere settings
        bodySettings.at("Earth")->atmosphereSettings = std::make_shared< AtmosphereSettings >( nrlmsise00 );

        // ---------- Create planets objects ----------
        // NamedBodyMap == std::unordered_map< std::string, std::shared_ptr< simulation_setup::Body > >
        SystemOfBodies bodyMap = createSystemOfBodies(bodySettings);

        // ---------- Create vehicle object ----------
        // Set accelerations for each satellite.
        std::string currentSatelliteName;
        for ( unsigned int i = 0; i < numberOfSatellites; i++ ){
                currentSatelliteName =  "sat" + boost::lexical_cast< std::string >( i );
                bodyMap.createEmptyBody(currentSatelliteName);
                bodyMap.at(currentSatelliteName) = std::make_shared< simulation_setup::Body >( );
        }

        //---------- Create aerodynamic coefficient interface settings ----------
        std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
                std::make_shared< ConstantAerodynamicCoefficientSettings >(
                SAT_Ad, SAT_Cd * Eigen::Vector3d::UnitX( ), 1, 1 );
        // Create and set aerodynamic coefficients object
        for ( unsigned int i = 0; i < numberOfSatellites; i++ ){
                currentSatelliteName =  "sat" + boost::lexical_cast< std::string >( i );
                bodyMap.at(currentSatelliteName)->setAerodynamicCoefficientInterface(createAerodynamicCoefficientInterface(
                aerodynamicCoefficientSettings,currentSatelliteName));
        }

        //---------- Create radiation pressure settings ----------
        std::vector< std::string > occultingBodies;
        occultingBodies.push_back("Earth");
        std::shared_ptr< RadiationPressureInterfaceSettings > satRadiationPressureSettings =
                std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", SAT_SRPA, SAT_Cr, occultingBodies);
        // Create and set radiation pressure settings
        for ( unsigned int i = 0; i < numberOfSatellites; i++ ){
                currentSatelliteName =  "sat" + boost::lexical_cast< std::string >( i );
                bodyMap.at(currentSatelliteName)->setRadiationPressureInterface("Sun", createRadiationPressureInterface(satRadiationPressureSettings, currentSatelliteName, bodyMap));
        }

        // Finalize body creation
        setGlobalFrameBodyEphemerides( bodyMap.getMap(), "Earth", "J2000" );
        
        
        
        // --------------------------------------------------------
        // ---------- Define constellation initial state ----------
        // --------------------------------------------------------
        /* // Set orbital parameters of Galileo constellation.
        const double semiMajorAxis = CONSTELLATION_SMA;
        const double eccentricity = CONSTELLATION_ECC;
        const double inclination = CONSTELLATION_INC;
        const double argumentOfPeriapsis = CONSTELLATION_AOP;
        const double longitudeOfAscendingNodeSpacing = 2.0*mathematical_constants::PI/numberOfPlanes;
        const double trueAnomalySpacing = 2.0*mathematical_constants::PI/numberOfSatellitesPerPlane;
        const double adjacentPlaneSpacing = CONSTELLATION_F*2.0*mathematical_constants::PI/numberOfSatellites;
        // Declare size of state.
        const unsigned int sizeOfState = 6;
        // Set Keplerian elements
        Eigen::MatrixXd initialConditionsInKeplerianElements(sizeOfState,numberOfSatellites);
        initialConditionsInKeplerianElements.row(0) = Eigen::MatrixXd::Constant(1,numberOfSatellites,semiMajorAxis);
        initialConditionsInKeplerianElements.row(1) = Eigen::MatrixXd::Constant(1,numberOfSatellites,eccentricity);
        initialConditionsInKeplerianElements.row(2) = Eigen::MatrixXd::Constant(1,numberOfSatellites,inclination);
        initialConditionsInKeplerianElements.row(3) = Eigen::MatrixXd::Constant(1,numberOfSatellites,argumentOfPeriapsis);
         // Set longitude of ascending node.
        for (unsigned int i = 0; i < numberOfPlanes; i++ ){
                initialConditionsInKeplerianElements.block(3,i*numberOfSatellitesPerPlane,1,numberOfSatellitesPerPlane) = Eigen::MatrixXd::Constant(1,numberOfSatellitesPerPlane,i*adjacentPlaneSpacing);
        }
        // Set longitude of ascending node.
        for (unsigned int i = 0; i < numberOfPlanes; i++ ){
                initialConditionsInKeplerianElements.block(4,i*numberOfSatellitesPerPlane,1,numberOfSatellitesPerPlane) = Eigen::MatrixXd::Constant(1,numberOfSatellitesPerPlane,i*longitudeOfAscendingNodeSpacing);
        }
        // Set true anomaly.
        Eigen::RowVectorXd trueAnomalySpacingIntegers(numberOfSatellitesPerPlane);
        for ( unsigned int i = 0; i < numberOfSatellitesPerPlane; i++ ){
                trueAnomalySpacingIntegers(i) =  i*1.0;
        }
        for ( unsigned int i = 0; i < numberOfPlanes; i++ ){
                initialConditionsInKeplerianElements.block(5,i*numberOfSatellitesPerPlane, 1, numberOfSatellitesPerPlane ) = Eigen::MatrixXd::Constant( 1, numberOfSatellitesPerPlane, trueAnomalySpacing ).array( )*trueAnomalySpacingIntegers.array( );
        }
        
        // Convert initial conditions to Cartesian elements.
        Eigen::MatrixXd initialConditions(sizeOfState, numberOfSatellites);
        double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        for ( unsigned int i = 0; i < numberOfSatellites; i++ ){
                Eigen::Vector6d initKepl = initialConditionsInKeplerianElements.col( i ).cast< double >();
                initialConditions.col( i ) = convertKeplerianToCartesianElements(initKepl, static_cast< double >(earthGravitationalParameter));
        }

        Eigen::Matrix<double,Eigen::Dynamic,1> systemIterationInitialState = Eigen::MatrixXd::Zero(sizeOfState * numberOfSatellites,1);
        for (unsigned int i = 0; i < numberOfSatellites; i++ ){
                systemIterationInitialState.segment( i * sizeOfState, sizeOfState ) = initialConditions.col( i );
        }*/
        
        // Declare size of state.
        const unsigned int sizeOfState = 6;
        // Load x0 from txt file
        //std::ofstream x0file_;
        std::string filepath(X0_FILE_PATH);
        
        // Open file
    	//x0file_.open(X0_FILE_PATH,std::ios::in);
    	std::ifstream x0file_(filepath);
    	// Initial state vectors
        Eigen::Matrix<double,Eigen::Dynamic,1> systemIterationInitialState = Eigen::MatrixXd::Zero(sizeOfState * numberOfSatellites,1);
        Eigen::Matrix<double,numberOfSatellites,1> systemIterationInitialBodyMasses;     
        double aux;
        for (unsigned int i = 0; i < numberOfSatellites; i++ ){
        	// Initial position + velocity in cartesian coordinates
        	for (unsigned int j = 0; j < sizeOfState; j++ ){
        		x0file_ >> aux;
        		systemIterationInitialState(i*sizeOfState+j,0) = aux;
        	}
        	// Initial mass
        	x0file_ >> aux;
        	systemIterationInitialBodyMasses(i,0) = aux;
        }
        // Close file
        x0file_.close();
        
        
        // -----------------------------------------------
        // ---------- Setup acceleration models ----------
        // -----------------------------------------------

        // SelectedAccelerationMap == std::map< std::string, std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > >
        SelectedAccelerationMap accelerationMap;
        std::vector <std::string> bodiesToPropagate;
        std::vector <std::string> centralBodies;

        // ----------- Thrust feedback
        // 4 evaluations per interval in RK4s
        std::shared_ptr <tudatThrustFeedbackMatlab> ttfm = std::make_shared <tudatThrustFeedbackMatlab>(SERVER_PORT,numberOfSatellites,EPOCH_CONTROL_UPDATE,5*EPOCH_CONTROL_UPDATE/EPOCH_SAMPLE,SAT_ISP*SAT_g0,simulationEndEpoch);

        // Create acceleration models
        // Set accelerations for each satellite.
        for ( unsigned int i = 0; i < numberOfSatellites; i++ ){
                currentSatelliteName = "sat" + boost::lexical_cast< std::string >( i );
                std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfCurrentSatellite;
                // Central gravity + Spherical Harmonics + Atmosferic drag 
                accelerationsOfCurrentSatellite["Earth"] = {sphericalHarmonicAcceleration( 24, 24 ),aerodynamicAcceleration()};
                // Solar radiation pressure
                accelerationsOfCurrentSatellite["Sun"] = {pointMassGravityAcceleration(),cannonBallRadiationPressureAcceleration( )};
                // Third body gravity
                accelerationsOfCurrentSatellite["Moon"] = {pointMassGravityAcceleration()};
                accelerationsOfCurrentSatellite["Mars"] = {pointMassGravityAcceleration()};
                accelerationsOfCurrentSatellite["Venus"] = {pointMassGravityAcceleration()};    
                accelerationsOfCurrentSatellite["Jupiter"] = {pointMassGravityAcceleration()}; 
                
                //Thrust - Custom function                              
                std::function <Eigen::Vector3d(const double)> thrustFeedback = std::bind(&tudatThrustFeedbackMatlab::thrustFeedbackWrapper, ttfm, std::placeholders::_1,i,bodyMap.getMap());
                accelerationsOfCurrentSatellite[currentSatelliteName] = {std::make_shared<ThrustAccelerationSettings>(	thrustFeedback,
                																										SAT_ISP,
                																										THRUST_FRAME,
                																										"Earth")};
                          
                // Push satellite acceleration map
                accelerationMap[currentSatelliteName] = accelerationsOfCurrentSatellite;
                // Define body to propagate
                bodiesToPropagate.push_back(currentSatelliteName);
                 // Define central body
                centralBodies.push_back("Earth");
        }
        
        // Create acceleration models and propagation settings
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(bodyMap,accelerationMap,bodiesToPropagate,centralBodies);

        // ------------------------------------------------
        // ---------- Setup propagation settings ----------
        // ------------------------------------------------
        // ---------- Translational propagator settings ----------      
        // Define propagation termination conditions
        std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
               std::make_shared< propagators::PropagationTimeTerminationSettings >(simulationEndEpoch);
        std::shared_ptr <TranslationalStatePropagatorSettings <double>> translationalPropagatorSettings = 
                std::make_shared <TranslationalStatePropagatorSettings<double>>(centralBodies,accelerationModelMap,bodiesToPropagate,systemIterationInitialState,simulationEndEpoch);
        std::map <double,Eigen::VectorXd> integrationResult;
        
        // ---------- Mass propagator settings ----------
        // Create mass rate models
        std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
        std::vector< std::string > bodiesWithMassToPropagate;
        //Eigen::Matrix<double,numberOfSatellites,1> systemIterationInitialBodyMasses;
        for ( unsigned int i = 0; i < numberOfSatellites; i++ ){
                currentSatelliteName = "sat" + boost::lexical_cast< std::string >( i );
                //Mass rate for 6 thrusters aligned with TWN frame axis - Custom function                             
                std::function <double(double)> massRateThrustFeedback= std::bind(&tudatThrustFeedbackMatlab::massRateThrustFeedbackWrapper, ttfm, std::placeholders::_1,i,bodyMap.getMap());
                massRateModels[currentSatelliteName] = createMassRateModel(currentSatelliteName, std::make_shared< CustomMassRateSettings >( massRateThrustFeedback ),bodyMap,accelerationModelMap);
                bodiesWithMassToPropagate.push_back(currentSatelliteName);
                //systemIterationInitialBodyMasses(i,0) = SAT_MASS*1.0;
        } 
        // Create settings for propagating the mass of the vehicle
        std::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings = std::make_shared< MassPropagatorSettings< double > >(bodiesWithMassToPropagate, massRateModels,systemIterationInitialBodyMasses,terminationSettings);
        
        // ---------- Global state propagator ----------
        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
        propagatorSettingsVector.push_back(massPropagatorSettings);
        propagatorSettingsVector.push_back(translationalPropagatorSettings);
               
        // ---------- Define dependent variables ----------
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
        std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
                std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );
        std::shared_ptr< PropagatorSettings< > > propagatorSettings =
                std::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsVector, terminationSettings, dependentVariablesToSave );

        // Initialialize integrator
        std::shared_ptr <IntegratorSettings <>> integratorSettings = std::make_shared <IntegratorSettings <>>(rungeKutta4,simulationStartEpoch,EPOCH_SAMPLE);
        
        // ------------------------------------------------
        // ---------- Propagate ---------------------------
        // ------------------------------------------------
        std::cout << "Started simulation.\n"; 
        // ---------- Declarion of loop variables ----------
        // Setup dynamics simulator
        SingleArcDynamicsSimulator <> dynamicsSimulator(bodyMap,integratorSettings,propagatorSettings,true);
        //std::map <double,Eigen::VectorXd> integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution();
        integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution();
        // Close feedback connection
        ttfm->closeConnection();
       
        std::cout << "Saving results.\n"; 
        // ---------- Save data ----------
        // Retrieve numerically integrated state for each satellite.
        std::vector< std::map< double, Eigen::VectorXd > > allSatellitesPropagationHistory;
        allSatellitesPropagationHistory.resize( numberOfSatellites );
        for(std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( ); stateIterator != integrationResult.end( ); stateIterator++ ){
                for( unsigned int i = 0; i < allSatellitesPropagationHistory.size( ); i++ ){
                        allSatellitesPropagationHistory[i][stateIterator->first] = Eigen::VectorXd::Zero(7);
                        allSatellitesPropagationHistory[i][stateIterator->first].segment(0,6) = stateIterator->second.segment(i*6,6);
                        allSatellitesPropagationHistory[i][stateIterator->first](6) = stateIterator->second(numberOfSatellites*6+i);
                }
        }
        for (unsigned int i = 0; i < numberOfSatellites; i++){
                // Set filename for output data.
                std::stringstream outputFilename;
                outputFilename << "stateSat" << i << ".dat";
                // Write satellite propagation history to file.
                writeDataMapToTextFile( allSatellitesPropagationHistory.at(i).begin(),allSatellitesPropagationHistory.at(i).end(),
                        outputFilename.str( ),
                        "",
                        std::numeric_limits< double >::digits10,
                        std::numeric_limits< double >::digits10,
                        "\t" );
        }
        std::cout << "Simulation complete.\n"; 
         
    return EXIT_SUCCESS;
}
