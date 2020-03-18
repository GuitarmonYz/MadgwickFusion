#include <math.h>

struct RMSFollower
{
	float y_pre;
	float alpha_attack;
	float alpha_release;
} rms_follower_slow, rms_follower_fast;

float threshold_;
float sensitivity_;
float sample_rate_ = 10000.0f;

inline float rms_follower_process(struct RMSFollower* follower_ptr, float x_in);
inline float origin_to_poly_skew(float new_val, float skew_factor, float skew_mag, float min, float max);
inline void quadratic_solver(float a, float b, float c, float* res);
inline float origin_to_linear_skew (float sourceValue, float sourceRangeMin, float sourceRangeMax, float targetRangeMin, float targetRangeMax);

void transient_shaper_init(float attack, float attack_ratio, float release, float threshold)
{
	threshold_ = threshold;
	
	rms_follower_slow.y_pre = 1.0f;
	rms_follower_slow.alpha_attack = exp(-1e3/(attack * sample_rate_));
	rms_follower_slow.alpha_release = exp(-1e3/(release * sample_rate_));
	
	rms_follower_fast.y_pre = 1.0f;
	rms_follower_fast.alpha_attack = exp(-1e3/(attack * sample_rate_ * attack_ratio / 100.0f));
	rms_follower_fast.alpha_release = exp(-1e3/(release * sample_rate_));
	
	sensitivity_ = 4.0f;
}

float transient_shaper_get_gain(float sample_in)
{
	float y_slow = rms_follower_process(&rms_follower_slow, sample_in);
	float y_fast = rms_follower_process(&rms_follower_fast, sample_in);
	
    float diff = y_fast - y_slow;
    float diff_db = 20.0f * log10(diff);
    float gain = 1.0f;
    if (diff_db > threshold_)
    {
    //    gain = origin_to_poly_skew(diff, 0.0f, 1.0f, 0.0f, 1.0f);
       gain = origin_to_linear_skew(gain, 0.0f, 1.0f, 1.0f, sensitivity_);
    }
    return gain;
}

//Below are helper functions
float rms_follower_process(struct RMSFollower* follower_ptr, float x_in)
{
	x_in = fabs(x_in);
	float y = 0.0f;
	if (x_in > follower_ptr->y_pre)
	{
		y = sqrt(pow(follower_ptr->y_pre, 2) * follower_ptr->alpha_attack + pow(x_in, 2) * (1.0f - follower_ptr->alpha_attack));
	}
	else
	{
		y = sqrt(pow(follower_ptr->y_pre, 2) * follower_ptr->alpha_release + pow(x_in, 2) * (1.0f - follower_ptr->alpha_release));
	}
	follower_ptr->y_pre = y;
	return y;
}

float origin_to_poly_skew(float new_val, float skew_factor, float skew_mag, float min, float max)
{
	float new_val_01 = origin_to_linear_skew(new_val, min, max, 0.0f, 1.0f);
	float skew_point_x = skew_factor;
	float skew_point_y = skew_mag;
	float a = 1 - 2 * skew_point_x;
	float b = 2 * skew_point_x;
	float c = -new_val_01;
	float res[2] = {0, 0};
	quadratic_solver(a, b ,c, res);
	float t = res[0];
	float y_01 = 2 * (1 - t) * t * skew_point_y + t * t;
	return origin_to_linear_skew(y_01, 0.0f, 1.0f, min, max);
}

void quadratic_solver(float a, float b, float c, float* res)
{
    if (a < (float)1e-4 && a > (float)-1e-4)
    {
        res[0] = -c/b;
        res[1] = -c/b;
        return;
    }
    res[0] = (-b + sqrtf(b*b - 4*a*c)) / (2*a);
    res[1] = (-b - sqrtf(b*b - 4*a*c)) / (2*a);
}

float origin_to_linear_skew (float src_val, float src_min, float src_max, float tgt_min, float tgt_max)
{
    return tgt_min + ((tgt_max - tgt_min) * (src_val - src_min)) / (src_max - src_min);
}
