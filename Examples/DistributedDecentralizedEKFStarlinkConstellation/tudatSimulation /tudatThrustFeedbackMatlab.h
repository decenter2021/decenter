#ifndef TUDAT_THRUST_FEEDBACK_MATLAB_H
#define TUDAT_THRUST_FEEDBACK_MATLAB_H

// Standard headers
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <stdexcept>

#define DISCONNECTED_SOCKET -11
#define COM_SUCCESS 0
#define COM_ERROR -1
#define MAX_ATOMIC_LENGTH 1024*8

#include <Eigen/Core>
#include "applicationOutput.h"
#include "simulationParameters.h"

class tudatThrustFeedbackMatlab{
public:
    // ---------- Construct / Destructor ----------
    tudatThrustFeedbackMatlab(const unsigned int numberOfSatellites, const double epochControlUpdatePeriod, const int evaluationsPerUpdatePeriod, const double isp_times_g0, double simulationEndEpoch = -1){
        std::cout << "@Tudat: Creating MATLAB thrust feedback object for " << numberOfSatellites << " satellites." << std::endl;
        // Init number of satellites and dynamic matrices 
        this->numberOfSatellites_ = numberOfSatellites;
        this->simulationEndEpoch_ = simulationEndEpoch*1.0;
        // Init control update peiod
        this->epochControlUpdatePeriod_ = epochControlUpdatePeriod;
        this->endTimeInstant_ = (int) (simulationEndEpoch_/epochControlUpdatePeriod + 0.5);
        this->evaluationsPerUpdatePeriod_ = (int) (evaluationsPerUpdatePeriod*numberOfSatellites+0.5);
        
        this->isp_times_g0_ = isp_times_g0;
        this->thrust_ = Eigen::MatrixXd::Zero(3,numberOfSatellites);
        
        this->establishConnectionMatlabUDP();
        // Init files to write output actuation
        this->initWriteActuation();
    }

    void closeConnection(){
        this->terminateConnectionMatlabUDP();
    }

    ~tudatThrustFeedbackMatlab(){
        this->terminateConnectionMatlabUDP();
    }

    // ---------- Feedback wrapper functions ----------
    Eigen::Vector3d thrustFeedbackWrapper(const double time, const unsigned int satId, const std::unordered_map< std::string, std::shared_ptr< tudat::simulation_setup::Body > >& bodyMap){
        // Check need to update control action
        int time_instant = (int) time/epochControlUpdatePeriod_; 
        if (time_instant == lastComputationTimeInstant_ && currentThrustEvaluations_ < evaluationsPerUpdatePeriod_){
        	currentThrustEvaluations_++;
        	//std::cout << "currentThrustEvaluations_ " << currentThrustEvaluations_ << std::endl;
        	return thrust_.col(satId);
        }else if (time_instant == lastComputationTimeInstant_ + 1){
        	if (currentThrustEvaluations_ < evaluationsPerUpdatePeriod_){
        		currentThrustEvaluations_++;
        		return thrust_.col(satId);
        	}else if (currentThrustEvaluations_ == evaluationsPerUpdatePeriod_){
        		if (time_instant != endTimeInstant_){
        			this->updateActuation(time,bodyMap);
        			lastComputationTimeInstant_ = time_instant;
        		}
        		currentThrustEvaluations_ = 1;
        		return thrust_.col(satId);
        	}else{       		
        		throw std::invalid_argument("Lost at least one sample of the thrust time function.");
        	}
        }else if (lastComputationTimeInstant_ < 0){
        	this->updateActuation(time,bodyMap);
        	lastComputationTimeInstant_ = time_instant;
        	currentThrustEvaluations_ = 1;
        	return thrust_.col(satId);
        }
        throw std::invalid_argument("Requested thrust for an invalid time-instant.");
    }
    
    // ---------- Feedback wrapper functions ----------
    double massRateThrustFeedbackWrapper(const double time, const unsigned int satId, const std::unordered_map< std::string, std::shared_ptr< tudat::simulation_setup::Body > >& bodyMap){       
        return -(this->thrust_).col(satId).cwiseAbs().sum()/(this->isp_times_g0_);
    }
    
    // ---------- Public variables ----------
    unsigned int numberOfSatellites_;
    double epochControlUpdatePeriod_;

protected:
    // ---------- Private Matlab interface functions ----------
    void updateActuation(const double time, const std::unordered_map< std::string, std::shared_ptr< tudat::simulation_setup::Body > >& bodyMap){
        // Define state struct to sent to Matlab 
        //Eigen::Matrix<double,7,Eigen::Dynamic> state = Eigen::MatrixXd::Zero(7,this->numberOfSatellites_);
        Eigen::Matrix<double,Eigen::Dynamic,1> state = Eigen::MatrixXd::Zero(7*this->numberOfSatellites_+1,1);
        state(0,0) = time;
        this->getStateFromBodyMap(&state,bodyMap);
        // send double time - 8 bytes
        // send state stateStruct - state.data() - 7*nsats*8 bytes
        this->getFeedback(&state);

        if (simulationEndEpoch_>0){
            this->progressBar(time/simulationEndEpoch_);
        }
    }

    // Get satellite state in ECI frame from body map
    void getStateFromBodyMap(Eigen::Matrix<double,Eigen::Dynamic,1> * state, const std::unordered_map< std::string, std::shared_ptr< tudat::simulation_setup::Body > >& bodyMap){
        //Eigen::Matrix<double, 6,1> earthSSBJ2000 = bodyMap.at("Earth")->getState();
        std::string currentSatelliteName;
        for (unsigned int i = 0; i < this->numberOfSatellites_; i++){
            currentSatelliteName = "sat" + boost::lexical_cast< std::string >(i);
            state->block(i*7+1,0,6,1) =  (bodyMap.at(currentSatelliteName)->getState());//-earthSSBJ2000;
            (*state)(i*7+7,0) = bodyMap.at(currentSatelliteName)->getBodyMass();
        }
    }

    int establishConnectionMatlabUDP(){
        // ---------- Setup client connection ----------
        // Create client socket 
        std::cout << "@Tudat: Entered establish connection UDP.\n";
        this->clientSocket_ = DISCONNECTED_SOCKET;
        if ((this->clientSocket_ = socket(AF_INET, SOCK_DGRAM,0)) == DISCONNECTED_SOCKET) {
            perror("@Tudat: Socket creation failed");
            return COM_ERROR;
        }
        // Setup server address
        struct sockaddr_in serverAddr;
        socklen_t lenServerAddr;
        memset(&(this->serverAddr_), 0, sizeof(this->serverAddr_));
        
        // Fill server information
        this->serverAddr_.sin_family = AF_INET;
        this->serverAddr_.sin_port = htons(SERVER_PORT);
        this->serverAddr_.sin_addr.s_addr = INADDR_ANY;
        
        // Ping server
        Eigen::Matrix<double,Eigen::Dynamic,1> state = Eigen::MatrixXd::Random(1+7*this->numberOfSatellites_,1);
        const char * pingMsg =  (const char *) state.data();
        
        unsigned int numberOfBytesToSend = (7*this->numberOfSatellites_+1)*sizeof(double);
        unsigned int numberOfBytesSent = 0;
        unsigned int lengthPacket = 0;
        int status;
        while(numberOfBytesToSend > 0){
            lengthPacket = MAX_ATOMIC_LENGTH;
            if(lengthPacket > numberOfBytesToSend){
                lengthPacket = numberOfBytesToSend;
            }
            int status = sendto(this->clientSocket_,  (void *)(pingMsg+numberOfBytesSent),lengthPacket,
                MSG_CONFIRM, (const struct sockaddr *) &(this->serverAddr_), 
                sizeof(this->serverAddr_)); 
            numberOfBytesSent += lengthPacket;
            numberOfBytesToSend -= lengthPacket;
        }
        std::cout << "@Tudat: Sent ping.\n";       
        //std::cout << "@Tudat: Ended ping " << status << " " << EMSGSIZE <<  " \n";
        
        unsigned int numberOfBytesToReceive = 3*this->numberOfSatellites_*sizeof(double);
        char * byteStream = (char *) malloc((this->numberOfSatellites_)*3*sizeof(double));
        unsigned int numberOfBytesReceived = 0;
        while(numberOfBytesToReceive > 0){
            int n = recv(this->clientSocket_, (void *)(byteStream + numberOfBytesReceived), numberOfBytesToReceive, 
                MSG_WAITALL);
            numberOfBytesToReceive -= n;
            numberOfBytesReceived += n;
        }
        free(byteStream);
    
        std::cout << "@Tudat: MATLAB UDP connection established.\n";
        std::flush(std::cout);
        return COM_SUCCESS;
    }

    int terminateConnectionMatlabUDP(){
        std::cout << std::endl; // not to overwrite progress bar
        if(this->clientSocket_ == DISCONNECTED_SOCKET){
            std::cout << "Error terminating connection: Socket disconnected.\n";
            return COM_ERROR;
        }
       
        Eigen::Matrix<double,Eigen::Dynamic,1> state = Eigen::MatrixXd::Random(1+7*this->numberOfSatellites_,1);
        state(0,0) = -1;
        const char * pingMsg = (const char *) state.data();
        unsigned int numberOfBytesToSend = (7*this->numberOfSatellites_+1)*sizeof(double);
        unsigned int numberOfBytesSent = 0;
        unsigned int lengthPacket = 0;
        int status;
        while(numberOfBytesToSend > 0){
            lengthPacket = MAX_ATOMIC_LENGTH;
            if(lengthPacket > numberOfBytesToSend){
                lengthPacket = numberOfBytesToSend;
            }
            status = sendto(this->clientSocket_,  (void *)(pingMsg+numberOfBytesSent),lengthPacket,
                MSG_CONFIRM, (const struct sockaddr *) &(this->serverAddr_), 
                sizeof(this->serverAddr_)); 
            numberOfBytesSent += lengthPacket;
            numberOfBytesToSend -= lengthPacket;
        }
        
        // Close client socket 
        close(this->clientSocket_);
        this->clientSocket_ = DISCONNECTED_SOCKET;
        return COM_SUCCESS;
    }

    int getFeedback(Eigen::Matrix<double,Eigen::Dynamic,1> * state){
        if(this->clientSocket_ == DISCONNECTED_SOCKET){
            std::cout << "Error getting feedback: Socket disconnected.\n";
            return COM_ERROR;
        }

        const char * msgOut = (const char *) state->data();
        
        unsigned int numberOfBytesToSend = (7*this->numberOfSatellites_+1)*sizeof(double);
        unsigned int numberOfBytesSent = 0;
        unsigned int lengthPacket = 0;
        int status;
        while(numberOfBytesToSend > 0){
            lengthPacket = MAX_ATOMIC_LENGTH;
            if(lengthPacket > numberOfBytesToSend){
                lengthPacket = numberOfBytesToSend;
            }
            status = sendto(this->clientSocket_,  (void *)(msgOut+numberOfBytesSent),lengthPacket,
                MSG_CONFIRM, (const struct sockaddr *) &(this->serverAddr_), 
                sizeof(this->serverAddr_)); 
            numberOfBytesSent += lengthPacket;
            numberOfBytesToSend -= lengthPacket;
        }

        unsigned int numberOfBytesToReceive = 3*this->numberOfSatellites_*sizeof(double);
        char * byteStream = (char *) malloc((this->numberOfSatellites_)*3*sizeof(double));
        unsigned int numberOfBytesReceived = 0;
        while(numberOfBytesToReceive > 0){
            int n = recv(this->clientSocket_, (void *)(byteStream + numberOfBytesReceived), numberOfBytesToReceive, 
                MSG_WAITALL);
            numberOfBytesToReceive -= n;
            numberOfBytesReceived += n;
        }
        double * doubleStream = (double * ) byteStream;

        this->writeActuation((*state)(0,0),doubleStream);
        
        // Post-process actuation
        Eigen::Matrix<double,3,1> thrust; 
        for (unsigned int i = 0; i < this->numberOfSatellites_; i++){
            this->thrust_.col(i) = Eigen::Map<Eigen::Matrix<double,3,1>>(doubleStream+i*3);
        }
        free(byteStream);

        return COM_SUCCESS;
    }

    void progressBar(float progress){
        int barWidth = 70;
        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();
    }
    
    void initWriteActuation(){
    	for(unsigned int i = 0; i < this->numberOfSatellites_; i++){           
            std::stringstream outputFilename;
            // Check if output directory exists; create it if it doesn't.
            if ( !boost::filesystem::exists( "./ouput/" ) ){
                boost::filesystem::create_directories( "./output/" );
            }
            outputFilename << "./output/inputSat" << i << ".dat";
            // Open output file.
            std::ofstream outputFile_;
            outputFile_.open(outputFilename.str(),std::ios::out);
            outputFile_ << "";
            outputFile_.close();
        }
    }

    void writeActuation(const double time, double * u){ 
        for(unsigned int i = 0; i < this->numberOfSatellites_; i++){
            std::stringstream outputFilename;
            outputFilename << "./output/inputSat" << i << ".dat";
            // Open output file.
            std::ofstream outputFile_;
            outputFile_.open(outputFilename.str(), std::ios::app);
            outputFile_ << time << "\t";
            for(int j = 0; j<3;j++){
                outputFile_ << *(u+3*i+j); 
                if(j<2) outputFile_ << "\t";
            }
            outputFile_ << std::endl;
            outputFile_.close();
        }
    }   

    // ---------- Private matlab interface variables ----------
    // Init last computation time to invalid vaue to force an update
    int lastComputationTimeInstant_ = -2;
    int endTimeInstant_;
    double simulationEndEpoch_ = -1;
    Eigen::MatrixXd thrust_;
    int clientSocket_;
    struct sockaddr_in serverAddr_;
    socklen_t lenServerAddr_;
    int evaluationsPerUpdatePeriod_;
    int currentThrustEvaluations_ = 0;
    double isp_times_g0_;
};

#endif

