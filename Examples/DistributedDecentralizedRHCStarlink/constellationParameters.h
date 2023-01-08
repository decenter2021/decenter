#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H 

// ---------- All parameters in SI units ----------
// ---------- Simulation pamake allrameters ------------
// Epochs
// EPOCH_START must be a multile of EPOCH_CONTROL_UPDATE
#define EPOCH_START 0.0
#define EPOCH_END 77000.0
#define EPOCH_SAMPLE 10.0
#define EPOCH_CONTROL_UPDATE 10.0
#define ENVIRONMENT_TIME_BUFFER 300.0

// Thrust frame
#define THRUST_FRAME tnw_thrust_frame
// #define THRUST_FRAME inertial_thrust_frame

// ---------- Constellation parameters ------------
// Atmospheric drag
#define SAT_Cd 2.2
#define SAT_Ad 24.0

// Thrust
#define SAT_MASS 260.0
#define SAT_Ct1 0.068
#define SAT_ISP 1640.0
#define SAT_g0 9.81

// SRP
#define SAT_Cr 1.2
#define SAT_SRPA 10.0

// Nominal constellation definition
#define CONSTELLATION_N_PLANES 72
#define CONSTELLATION_N_PER_PLANE 22
#define CONSTELLATION_F 17
#define CONSTELLATION_SMA 6921000.0
#define CONSTELLATION_INC convertDegreesToRadians(53.0)
#define CONSTELLATION_INC_DEG 53.0
#define CONSTELLATION_ECC 0.0
#define CONSTELLATION_AOP 0.0
#define X0_FILE_PATH "./data/x0/x0_tmp4_1584.txt"

// ---------- Server parameters ------------
#define SERVER_PORT 6038

#endif
