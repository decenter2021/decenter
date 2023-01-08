/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
*/
/*
 * Modified version of 'applicationOutput.h' available in 
 * https://github.com/Tudat/tudatExampleApplications 
 * */
 
#ifndef TUDAT_APPLICATIONOUTPUT_H
#define TUDAT_APPLICATIONOUTPUT_H

#include <iostream>


template< typename InputIterator >
void writeDataMapToTextFile(
                InputIterator iteratorDataMap, InputIterator last,
                const std::string& filename,
                const std::string& fileHeader,
                const int precisionOfKeyType, const int precisionOfValueType,
                const std::string& delimiter ){
                
	std::stringstream outputFilename;
	// Check if output directory exists; create it if it doesn't.
    if ( !boost::filesystem::exists( "./ouput/" ) ){
    	boost::filesystem::create_directories( "./output/" );
    }
    outputFilename << "./output/" << filename;
    // Open output file.
    // Write file header to file.
    std::ofstream outputFile_;
    outputFile_.open(outputFilename.str(),std::ios::out);

    // Write file header to file.
    outputFile_ << fileHeader;
 
    // Loop over map of propagation history.
    for ( ; iteratorDataMap != last; iteratorDataMap++ ){
                // Print map data to output file.
                outputFile_ << std::setprecision( precisionOfKeyType )
                     << std::left << std::setw( precisionOfKeyType + 1 )
                     << iteratorDataMap->first;
                tudat::input_output::writeValueToStream( outputFile_, iteratorDataMap->second, precisionOfValueType,
                             delimiter );
	}
        // Close output file.
        outputFile_.close( );
}


#endif // TUDAT_APPLICATIONOUTPUT_H
