#include "gradient_optimize.h"
#include "transient_shaper.h"
#include <stdint.h>

#define M_PI 3.14159265358979323846f
#define GYRO_MEAS_ERROR M_PI * (10.0f / 180.0f)
#define GYRO_DRIFT_ERROR M_PI * (0.01f / 180.0f)
#define alpha 10.0f
#define beta sqrtf(3.0f / 4.0f) * GYRO_MEAS_ERROR
#define zeta sqrtf(3.0f / 4.0f) * GRYO_DRIFT_ERROR
#define MADGWICK_ITER_TIMES 5

// Madgwick fusion algo setup
float deltat = 0.0f;
uint64_t now = 0, lastUpdate = 0;

float g_qmag = 1.0f;
float gamma = 0.0f;

float ax=0.0f, ay=0.0f, az=0.0f, gx=0.0f, gy=0.0f, gz=0.0f, mx=0.0f, my=0.0f, mz=0.0f;
float q[4] = {1.0f, 0.0f, 0.0f, 0.0f};
float yaw_pitch_roll[3] = {0.0f, 0.0f, 0.0f};

//float g_bx = 0.0f, g_by = 0.0f, g_bz = 0.0f;

// Raw data pre-processing
float r_ax, r_ay, r_az, r_gx, r_gy, r_gz, r_mx, r_my, r_mz;

// Functions
inline float inv_sqrt(float x);
inline void normalize3d(float* x, float* y, float* z);
inline void normalize4d(float* x, float* y, float* z, float* w);

float inv_sqrt(float x) {
	float halfx = 0.5f * x;
	float y = x;
	long i = *(long*)&y;
	i = 0x5f3759df - (i>>1);
	y = *(float*)&i;
	y = y * (1.5f - (halfx * y * y));
	y = y * (1.5f - (halfx * y * y));
	return y;
}

void normalize3d(float* x, float* y, float* z)
{
	float norm_invert = inv_sqrt(pow(*x, 2) + pow(*y, 2) + pow(*z, 2));
	*x *= norm_invert;
	*y *= norm_invert;
	*z *= norm_invert;
}

void normalize4d(float* x, float* y, float* z, float* w)
{
	float norm_invert = inv_sqrt(pow(*x, 2) + pow(*y, 2) + pow(*z, 2) + pow(*w, 2));
	*x *= norm_invert;
	*y *= norm_invert;
	*z *= norm_invert;
	*w *= norm_invert;
}

void quaternion_to_YPR(void)
{
	const float a12 = 2.0f * (q[1] * q[2] + q[0] * q[3]);
	const float a22 = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];;
	const float a31 = 2.0f * (q[0] * q[1] + q[2] * q[3]);
	const float a32 = 2.0f * (q[1] * q[3] - q[0] * q[2]);
	const float a33 = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];;
	float pitch = -asinf(a32);
	float roll = atan2f(a31, a33);
	float yaw = atan2f(a12, a22);
	
	pitch /= (M_PI / 2.0f);
	roll /= M_PI;
	yaw /= M_PI;
	
	if (pitch > 0.54f) 
	{
        roll += 0.65f * (pitch - 0.54f) * 1.3333f;
        yaw += 0.6f * (pitch - 0.54f) * 1.3333f;
    }

	yaw_pitch_roll[0] = yaw;
	yaw_pitch_roll[1] = pitch;
	yaw_pitch_roll[2] = roll;
}

void motion_calculation_init(void)
{
	// a transient_shaper maybe used to smooth the g_qmag, which is the quaternion changing rate
	// transient_shaper_init(2.0f, 10.0f, 20.0f, -24.0f);
}

// the gradient from the original paper, where reference vecs are g = [0,0,1] and b = [bx, 0, bz]  
void madgwick_gradient(float* g1, float* g2, float* g3, float* g4, float q1, float q2, float q3, float q4)
{
	float _2q1 = 2.0f * q1;
	float _2q2 = 2.0f * q2;
	float _2q3 = 2.0f * q3;
	float _2q4 = 2.0f * q4;
	float q1q1 = q1 * q1;
	float q1q2 = q1 * q2;
	float q1q3 = q1 * q3;
	float q1q4 = q1 * q4;
	float q2q2 = q2 * q2;
	float q2q3 = q2 * q3;
	float q2q4 = q2 * q4;
	float q3q3 = q3 * q3;
	float q3q4 = q3 * q4;
	float q4q4 = q4 * q4;
	float q1mx = q1 * mx;
	float q1my = q1 * my;
	float q1mz = q1 * mz;
	float q2mx = q2 * mx;
	// Reference direction of Earth's magnetic field
	float hx = mx * q1q1 - 2.0f * q1my * q4 + 2.0f * q1mz * q3 + mx * q2q2 + _2q2 * my * q3 + _2q2 * mz * q4 - mx * q3q3 - mx * q4q4;
	float hy = 2.0f * q1mx * q4 + my * q1q1 - 2.0f * q1mz * q2 + 2.0f * q2mx * q3 - my * q2q2 + my * q3q3 + _2q3 * mz * q4 - my * q4q4;
	float bx = 0.5f * sqrtf(hx * hx + hy * hy);
	float bz = -q1mx * q3 + q1my * q2 + 0.5f * mz * q1q1 + q2mx * q4 - 0.5f * mz * q2q2 + q3 * my * q4 - 0.5f * mz * q3q3 + 0.5f * mz * q4q4;
	// calculate objective function
	float f1 = 2.0f * q2q4 - 2.0f * q1q3 - ax;
	float f2 = 2.0f * q1q2 + 2.0f * q3q4 - ay;
	float f3 = 1.0f - 2.0f * q2q2 - 2.0f * q3q3 - az;
	float f4 = 2.0f * bx * (0.5f - q3q3 - q4q4) + 2.0f * bz * (q2q4 - q1q3) - mx;
	float f5 = 2.0f * bx * (q2q3 - q1q4) + 2.0f * bz * (q1q2 + q3q4) - my;
	float f6 = 2.0f * bx * (q1q3 + q2q4) + 2.0f * bz * (0.5f - q2q2 - q3q3) - mz;
	// calculate gradient
	*g1 = 2.0f * (-q3 * f1 + q2 * f2 - bz * q3 * f4 + (-bx * q4 + bz * q2) * f5 + bx * q3 * f6);
	*g2 = 2.0f * (q4 * f1 + q1 * f2 - _2q2 * f3 + bz * q4 * f4 + (bx * q3 + bz * q1) * f5 + (bx * q4 - _2q2 * bz) * f6);
	*g3 = 2.0f * (-q1 * f1 + q4 * f2 - _2q3 * f3 + (-_2q3 * bx - bz * q1) * f4 + (bx * q2 + bz * q4) * f5 + (bx * q1 - _2q3* bz) * f6);
	*g4 = 2.0f * (q2 * f1 + q3 * f2 + (-_2q4 * bx + bz * q2) * f4 + (-bx * q1 + bz * q3) * f5 + bx * q2 * f6); 
}

// implements https://www.sciencedirect.com/science/article/pii/S0888327019303012
// where a new ref vec for magnetic field is generated via cross product of measurement from acc and gyro
void improved_madgwick_gradient(float* g1, float* g2, float* g3, float* g4, float q1, float q2, float q3, float q4)
{
	float q1q1 = q1 * q1;
	float q1q2 = q1 * q2;
	float q1q3 = q1 * q3;
	float q1q4 = q1 * q4;
	float q2q2 = q2 * q2;
	float q2q3 = q2 * q3;
	float q2q4 = q2 * q4;
	float q3q3 = q3 * q3;
	float q3q4 = q3 * q4;
	float q4q4 = q4 * q4;
	// new vector that will subsitute magnetic vec
	float vx = ay * mz - az * my;
	float vy = az * mx - ax * mz;
	float vz = ax * my - ay * mx;
	// normalize new vec
	normalize3d(&vx, &vy, &vz);
	// calculate objective function
	float f1 = 2.0f * (q2q4 - q1q3) - ax;
	float f2 = 2.0f * (q1q2 + q3q4) - ay;
	float f3 = 1.0f - 2.0f * q2q2 - 2.0f * q3q3 - az;
	float f4 = 2.0f * (q1q4 + q2q3) - vx;
	float f5 = q1q1 - q2q2 + q3q3 - q4q4 - vy;
	float f6 = 2.0f * (q3q4 - q1q2) - vz;
	// calculate gradient
	*g1 = 2.0f * (-q3 * f1 + q2 * f2 + q4 * f4 + q1 * f5 - q2 * f6);
	*g2 = 2.0f * (q4 * f1 + q1 * f2 - 2.0f * q2 * f3 + q3 * f4 - q2 * f5 - q1 * f6);
	*g3 = 2.0f * (-q1 * f1 + q4 * f2 - 2.0f * q3 * f3 + q2 * f4 + q3 * f5 + q4 * f6);
	*g4 = 2.0f * (q2 * f1 + q3 * f2 + q1 * f4 - q4 * f5 + q3 * f6);
}

void MadgwickQuaternionUpdate(void)
{
	// platform dependant delta_t
	now = app_timer_cnt_get();
	deltat = (float)app_timer_cnt_diff_compute(now, lastUpdate) / APP_TIMER_CLOCK_FREQ;
	lastUpdate = now;
	
	float q1 = q[0], q2 = q[1], q3 = q[2], q4 = q[3];   // short name local variable for readability
	float s1, s2, s3, s4;

    // Normalise accelerometer measurement
	normalize3d(&ax, &ay, &az);

    // Normalise magnetometer measurement
	normalize3d(&mx, &my, &mz);
	
//	madgwick_gradient(&s1, &s2, &s3, &s4, q1, q2, q3, q4);
	improved_madgwick_gradient(&s1, &s2, &s3, &s4, q1, q2, q3, q4);

	// normalise step magnitude
	normalize4d(&s1, &s2, &s3, &s4);

//  // compensate gryo drift, optional since good sensor's bias is really small
//	float g_err_x = 2 * (q1 * s2 - q2 * s1 - q3 * s4 + q4 * s3);
//	float g_err_y = 2 * (q1 * s3 + q2 * s4 - q3 * s1 - q4 * s2);
//	float g_err_z = 2 * (q1 * s4 - q2 * s3 + q3 * s2 - q4 * s1);
//	
//	g_bx += g_err_x * deltat * zeta;
//	g_by += g_err_y * deltat * zeta;
//	g_bz += g_err_y * deltat * zeta;
//	
//	gx -= g_bx;
//	gy -= g_by;
//	gz -= g_bz;
	
	// quaternion estimation based on gryo integration
	float q_wdot1 = -q2 * gx - q3 * gy - q4 * gz;
	float q_wdot2 = q1 * gx + q3 * gz - q4 * gy;
	float q_wdot3 = q1 * gy - q2 * gz + q4 * gx;
	float q_wdot4 = q1 * gz + q2 * gy - q3 * gx;
	
	float q_g1 = 0.5f * q_wdot1 * deltat;
	float q_g2 = 0.5f * q_wdot2 * deltat;
	float q_g3 = 0.5f * q_wdot3 * deltat;
	float q_g4 = 0.5f * q_wdot4 * deltat;
	
	// calcuate convergence rate
    g_qmag = 0.5f * sqrtf(q_wdot1 * q_wdot1 + q_wdot2 * q_wdot2 + q_wdot3 * q_wdot3 + q_wdot4 * q_wdot4);
	float mue = alpha * g_qmag * deltat;

    //  vanilla decent
	float q_d1 = - mue * s1;
	float q_d2 = - mue * s2;
	float q_d3 = - mue * s3;
	float q_d4 = - mue * s4;
	
	// calcuate fution ratio
	gamma = beta / (g_qmag * alpha + beta);
	// fuse two estimation
	q1 += gamma * q_d1 + (1 - gamma) * q_g1;
	q2 += gamma * q_d2 + (1 - gamma) * q_g2;
	q3 += gamma * q_d3 + (1 - gamma) * q_g3;
	q4 += gamma * q_d4 + (1 - gamma) * q_g4;
    
	// normalise quaternion
	normalize4d(&q1, &q2, &q3, &q4);
	q[0] = q1;
	q[1] = q2;
	q[2] = q3;
	q[3] = q4;
}
// raw sensor data as input, output calculated quaternion
void process_raw_motion_data(accel_values_t* acc_ptr, gyro_values_t* gyro_ptr, magn_values_t* magn_ptr, QUAT_DATA* quat_dest_ptr)
{	
	int16_t motion_data[9];
	memcpy(&motion_data[0], &acc_ptr->x, 2);
	memcpy(&motion_data[1], &acc_ptr->y, 2);
	memcpy(&motion_data[2], &acc_ptr->z, 2);
	memcpy(&motion_data[3], &gyro_ptr->x, 2);
	memcpy(&motion_data[4], &gyro_ptr->y, 2);
	memcpy(&motion_data[5], &gyro_ptr->z, 2);
	memcpy(&motion_data[6], &magn_ptr->x, 2);
	memcpy(&motion_data[7], &magn_ptr->y, 2);
	memcpy(&motion_data[8], &magn_ptr->z, 2);
	
	r_ax = motion_data[0] * 2.0 / 32768.0;
	r_ay = motion_data[1] * 2.0 / 32768.0;
	r_az = motion_data[2] * 2.0 / 32768.0;
	r_gx = motion_data[3] * 500.0/32768.0 ;
    r_gy = motion_data[4] * 500.0/32768.0 ;
    r_gz = motion_data[5] * 500.0/32768.0;
    r_mx = motion_data[6] * 10.*4912./8190.;
    r_my = motion_data[7] * 10.*4912./8190.;
    r_mz = motion_data[8] * 10.*4912./8190.;
	
	// update based on your sensor registration
	ax = -r_ay;
	ay = -r_ax;
	az = r_az;
	gx = r_gy * M_PI / 180.0f;
	gy = r_gx * M_PI / 180.0f;
	gz = -r_gz * M_PI / 180.0f;
	mx = r_mx;
	my = r_my;
	mz = r_mz;

	for (uint8_t i = 0; i < MADGWICK_ITER_TIMES; i++) MadgwickQuaternionUpdate();

	quat_dest_ptr->w = q[0];
	quat_dest_ptr->x = q[1];
	quat_dest_ptr->y = q[2];
	quat_dest_ptr->z = q[3];
}
