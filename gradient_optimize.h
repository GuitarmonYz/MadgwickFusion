// implemented several extention to venilla gradient decent, may have better convergence rate with proper param setting 
#include <math.h>

struct LowpassParams
{
	float lp_alpha; 
	float lp_beta;
	float lp_y;
};

float lowpassUpdate(struct LowpassParams* p, float x) 
{
	p->lp_y = p->lp_alpha * x + p->lp_beta * p->lp_y;
	return p->lp_y;
}

float lowpassUpdate2(struct LowpassParams* p, float x)
{
	p->lp_y = p->lp_alpha * x * x + p->lp_beta * p->lp_y;
	return p->lp_y;
}

struct AdamParams
{
	float elp;
	struct LowpassParams m1_lp, m2_lp, m3_lp, m4_lp,  
						 v1_lp, v2_lp, v3_lp, v4_lp;
};

void adamUpdate(float* g1, float* g2, float* g3, float* g4, struct AdamParams* p, float s1, float s2, float s3, float s4, float mue)
{
	lowpassUpdate(&(p->m1_lp), s1);
	lowpassUpdate(&(p->m2_lp), s2);
	lowpassUpdate(&(p->m3_lp), s3);
	lowpassUpdate(&(p->m4_lp), s4);
	lowpassUpdate(&(p->v1_lp), s1);
	lowpassUpdate(&(p->v2_lp), s2);
	lowpassUpdate(&(p->v3_lp), s3);
	lowpassUpdate(&(p->v4_lp), s4);
	*g1 = - mue * p->m1_lp.lp_y / sqrt(p->v1_lp.lp_y + p->elp);
	*g2 = - mue * p->m1_lp.lp_y / sqrt(p->v1_lp.lp_y + p->elp);
	*g3 = - mue * p->m1_lp.lp_y / sqrt(p->v1_lp.lp_y + p->elp);
	*g4 = - mue * p->m1_lp.lp_y / sqrt(p->v1_lp.lp_y + p->elp);
}

struct AdaDeltaParams
{
	float elp;
	struct LowpassParams g1_lp, g2_lp, g3_lp, g4_lp,
						 s1_lp, s2_lp, s3_lp, s4_lp;
};

void adaDeltaUpdate(float* g1, float* g2, float* g3, float* g4, struct AdaDeltaParams* p, float s1, float s2, float s3, float s4, float mue)
{	
	lowpassUpdate(&(p->g1_lp), s1);
	lowpassUpdate(&(p->g2_lp), s2);
	lowpassUpdate(&(p->g3_lp), s3);
	lowpassUpdate(&(p->g4_lp), s4);
	
	float delta_q1 = sqrtf(p->s1_lp.lp_y + p->elp) / sqrtf(p->g1_lp.lp_y + p->elp) * s1;
	float delta_q2 = sqrtf(p->s2_lp.lp_y + p->elp) / sqrtf(p->g1_lp.lp_y + p->elp) * s2;
	float delta_q3 = sqrtf(p->s3_lp.lp_y + p->elp) / sqrtf(p->g1_lp.lp_y + p->elp) * s3;
	float delta_q4 = sqrtf(p->s4_lp.lp_y + p->elp) / sqrtf(p->g1_lp.lp_y + p->elp) * s4;
	
	lowpassUpdate(&(p->s1_lp), delta_q1);
	lowpassUpdate(&(p->s2_lp), delta_q2);
	lowpassUpdate(&(p->s3_lp), delta_q3);
	lowpassUpdate(&(p->s4_lp), delta_q4);
	
	*g1 = - mue * p->s1_lp.lp_y;
	*g2 = - mue * p->s2_lp.lp_y;
	*g3 = - mue * p->s3_lp.lp_y;
	*g4 = - mue * p->s4_lp.lp_y;
}
