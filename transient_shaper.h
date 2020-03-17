#include "transient_shaper.h"
#include <math.h>

struct ARFollower
{
	float y_pre;
	float alpha_attack;
	float alpha_release;
}arfollower_slow, arfollower_fast;

float threshold_;
float sensitivity_;
float sample_rate_ = 10000.0f;

//Helper function declarations
float arfollower_process(struct ARFollower* follower_ptr, float x_in);
float originToPolySkew(float new_val, float skew_factor, float skew_mag, float min, float max);
void QuadraticSolver(float a, float b, float c, float* res);
float jmap (float sourceValue, float sourceRangeMin, float sourceRangeMax, float targetRangeMin, float targetRangeMax);

void transient_shaper_init(float attack, float attack_ratio, float release, float threshold)
{
	threshold_ = threshold;
	
	arfollower_slow.y_pre = 1.0f;
	arfollower_slow.alpha_attack = exp(-1e3/(attack * sample_rate_));
	arfollower_slow.alpha_release = exp(-1e3/(release * sample_rate_));
	
	arfollower_fast.y_pre = 1.0f;
	arfollower_fast.alpha_attack = exp(-1e3/(attack * sample_rate_ * attack_ratio / 100.0f));
	arfollower_fast.alpha_release = exp(-1e3/(release * sample_rate_));
	
	sensitivity_ = 4.0f;
}

float transient_shaper_get_gain(float sample_in)
{
	float y_slow = arfollower_process(&arfollower_slow, sample_in);
	float y_fast = arfollower_process(&arfollower_fast, sample_in);
	
    float diff = y_fast - y_slow;
    float diff_db = 20.0f * log10(diff);
    float gain = 1.0f;
    if (diff_db > threshold_)
    {
//        gain = originToPolySkew(diff, 0.0f, 1.0f, 0.0f, 1.0f);
//        gain = jmap(gain, 0.0f, 1.0f, 1.0f, sensitivity_);
		gain = sensitivity_;
    }
    return gain;
}

//Helper functions
float arfollower_process(struct ARFollower* follower_ptr, float x_in)
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

float originToPolySkew(float new_val, float skew_factor, float skew_mag, float min, float max)
{
	float new_val_01 = jmap(new_val, min, max, 0.0f, 1.0f);
	float skew_point_x = skew_factor;
	float skew_point_y = skew_mag;
	float a = 1 - 2 * skew_point_x;
	float b = 2 * skew_point_x;
	float c = -new_val_01;
	float res[2] = {0, 0};
	QuadraticSolver(a, b ,c, res);
	float t = res[0];
	float y_01 = 2 * (1 - t) * t * skew_point_y + t * t;
	return jmap(y_01, 0.0f, 1.0f, min, max);
}

void QuadraticSolver(float a, float b, float c, float* res)
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

float jmap (float sourceValue, float sourceRangeMin, float sourceRangeMax, float targetRangeMin, float targetRangeMax)
{
    return targetRangeMin + ((targetRangeMax - targetRangeMin) * (sourceValue - sourceRangeMin)) / (sourceRangeMax - sourceRangeMin);
}

