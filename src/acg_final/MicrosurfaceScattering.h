#ifndef _MICROSURFACESCATTERING_
#define _MICROSURFACESCATTERING_

#include "pbrt.h"
#include "geometry.h"
using namespace pbrt;


/************* MICROSURFACE HEIGHT DISTRIBUTION *************/

/* API */
class MicrosurfaceHeight
{
public:
	// height PDF	
	virtual float P1(const float h) const=0; 
	// height CDF	
	virtual float C1(const float h) const=0; 
	// inverse of the height CDF
	virtual float invC1(const float U) const=0; 
};

/* Uniform height distribution in [-1, 1] */
class MicrosurfaceHeightUniform : public MicrosurfaceHeight
{
public:	
	// height PDF	
	virtual float P1(const float h) const; 
	// height CDF	
	virtual float C1(const float h) const; 
	// inverse of the height CDF
	virtual float invC1(const float U) const; 
};

/* Gaussian height distribution N(0,1) */
class MicrosurfaceHeightGaussian : public MicrosurfaceHeight
{
public:	
	// height PDF	
	virtual float P1(const float h) const; 
	// height CDF	
	virtual float C1(const float h) const; 
	// inverse of the height CDF
	virtual float invC1(const float U) const; 
};


/************* MICROSURFACE SLOPE DISTRIBUTION *************/

/* API */
class MicrosurfaceSlope
{
public:
	MicrosurfaceSlope()
	{}

public:
	// roughness
	// projected roughness in wi
	float alpha_i(const Vector3f& wi, const float m_alpha_x, const float m_alpha_y) const; 

public:
	// distribution of normals (NDF)	
	float D(const Vector3f& wm, const float m_alpha_x, const float m_alpha_y) const; 
	// distribution of visible normals (VNDF)
	float D_wi(const Vector3f& wi, const Vector3f& wm, const float m_alpha_x, const float m_alpha_y) const; 
	// sample the VNDF
	Vector3f sampleD_wi(const Vector3f& wi, const float U1, const float U2, const float m_alpha_x, const float m_alpha_y) const;

public:
	// distribution of slopes
	virtual float P22(const float slope_x, const float slope_y, const float m_alpha_x, const float m_alpha_y) const=0; 
	// Smith's Lambda function
	virtual float Lambda(const Vector3f& wi, const float alpha_u, const float alpha_v) const=0;
	// projected area towards incident direction
	virtual float projectedArea(const Vector3f& wi, const float alpha_u, const float alpha_v) const=0;
	// sample the distribution of visible slopes with alpha=1.0
	virtual Vector2f sampleP22_11(const float theta_i, const float U1, const float U2) const=0;
};

/* Beckmann slope distribution */
class MicrosurfaceSlopeBeckmann : public MicrosurfaceSlope
{
public:
	MicrosurfaceSlopeBeckmann()
	{}

	// distribution of slopes
	virtual float P22(const float slope_x, const float slope_y, const float m_alpha_x, const float m_alpha_y) const; 
	// Smith's Lambda function
	virtual float Lambda(const Vector3f& wi, const float alpha_u, const float alpha_v) const;
	// projected area towards incident direction
	virtual float projectedArea(const Vector3f& wi, const float alpha_u, const float alpha_v) const;
	// sample the distribution of visible slopes with alpha=1.0
	virtual Vector2f sampleP22_11(const float theta_i, const float U1, const float U2) const;
};

/* GGX slope distribution */
class MicrosurfaceSlopeGGX : public MicrosurfaceSlope
{
public:
	MicrosurfaceSlopeGGX()
	{}

	// distribution of slopes
	virtual float P22(const float slope_x, const float slope_y, const float m_alpha_x, const float m_alpha_y) const; 
	// Smith's Lambda function
	virtual float Lambda(const Vector3f& wi, const float alpha_u, const float alpha_v) const;
	// projected area towards incident direction
	virtual float projectedArea(const Vector3f& wi, const float alpha_u, const float alpha_v) const;
	// sample the distribution of visible slopes with alpha=1.0
	virtual Vector2f sampleP22_11(const float theta_i, const float U1, const float U2) const;
};


/************* MICROSURFACE *************/

/* API */
class Microsurface 
{
public:
	// height distribution
	const MicrosurfaceHeight* m_microsurfaceheight; 
	// slope distribution
	const MicrosurfaceSlope* m_microsurfaceslope; 

public:

	Microsurface(const bool height_uniform, // uniform or Gaussian height distribution
				const bool slope_beckmann) : // Beckmann or GGX slope distribution  
		m_microsurfaceheight((height_uniform) ? 
		  static_cast<MicrosurfaceHeight*>(new MicrosurfaceHeightUniform) 
		: static_cast<MicrosurfaceHeight*>(new MicrosurfaceHeightGaussian)),
		m_microsurfaceslope((slope_beckmann) ? 
		  static_cast<MicrosurfaceSlope*>(new MicrosurfaceSlopeBeckmann()) 
		: static_cast<MicrosurfaceSlope*>(new MicrosurfaceSlopeGGX()))
	{}
		
	~Microsurface() 
	{
		delete m_microsurfaceheight; 
		delete m_microsurfaceslope;
	}

	// evaluate BSDF with a random walk (stochastic but unbiased)
	// scatteringOrder=0 --> contribution from all scattering events
	// scatteringOrder=1 --> contribution from 1st bounce only
	// scatteringOrder=2 --> contribution from 2nd bounce only, etc..
	virtual float eval(const Vector3f& wi, const Vector3f& wo, const int scatteringOrder=0) const; 

	// sample BSDF with a random walk
	// scatteringOrder is set to the number of bounces computed for this sample
	virtual Vector3f sample(const Vector3f& wi, int& scatteringOrder) const;
	Vector3f sample(const Vector3f& wi) const {int scatteringOrder; return sample(wi, scatteringOrder);}

public:
	// masking function
	float G_1(const Vector3f& wi) const;
	// masking function at height h0
	float G_1(const Vector3f& wi, const float h0) const;
	// sample height in outgoing direction
	float sampleHeight(const Vector3f& wo, const float h0, const float U) const;

public:
	// evaluate local phase function 
	virtual float evalPhaseFunction(const Vector3f& wi, const Vector3f& wo) const=0;
	// sample local phase function
	virtual Vector3f samplePhaseFunction(const Vector3f& wi) const=0; 

	// evaluate BSDF limited to single scattering 
	// this is in average equivalent to eval(wi, wo, 1);
	virtual float evalSingleScattering(const Vector3f& wi, const Vector3f& wo) const=0;

public:
	void refresh_alpha(float alphau, float alphav);

protected:
	float m_alphau, m_alphav;
};

/* Microsurface made of conductor material */
class MicrosurfaceConductor : public Microsurface
{
public:
	MicrosurfaceConductor(const bool height_uniform, // uniform or Gaussian
				 const bool slope_beckmann // Beckmann or GGX
)
		: Microsurface(height_uniform, slope_beckmann)
	{}

public:
	// evaluate local phase function 
	virtual float evalPhaseFunction(const Vector3f& wi, const Vector3f& wo) const;
	// sample local phase function
	virtual Vector3f samplePhaseFunction(const Vector3f& wi) const; 

	// evaluate BSDF limited to single scattering 
	// this is in average equivalent to eval(wi, wo, 1);
	virtual float evalSingleScattering(const Vector3f& wi, const Vector3f& wo) const; 
};

/* Microsurface made of conductor material */
class MicrosurfaceDielectric : public Microsurface
{
public:
	const float m_eta;
public:
	MicrosurfaceDielectric(const bool height_uniform, // uniform or Gaussian
				 const bool slope_beckmann, // Beckmann or GGX
				 const float eta = 1.5f)
		: Microsurface(height_uniform, slope_beckmann),
		m_eta(eta)
	{}

	// evaluate BSDF with a random walk (stochastic but unbiased)
	// scatteringOrder=0 --> contribution from all scattering events
	// scatteringOrder=1 --> contribution from 1st bounce only
	// scatteringOrder=2 --> contribution from 2nd bounce only, etc..
	virtual float eval(const Vector3f& wi, const Vector3f& wo, const int scatteringOrder=0) const; 

	// sample final BSDF with a random walk
	// scatteringOrder is set to the number of bounces computed for this sample
	virtual Vector3f sample(const Vector3f& wi, int& scatteringOrder) const;

public:
	// evaluate local phase function 
	virtual float evalPhaseFunction(const Vector3f& wi, const Vector3f& wo) const;
	float evalPhaseFunction(const Vector3f& wi, const Vector3f& wo, const bool wi_outside, const bool wo_outside) const;
	// sample local phase function
	virtual Vector3f samplePhaseFunction(const Vector3f& wi) const; 
	Vector3f samplePhaseFunction(const Vector3f& wi, const bool wi_outside, bool& wo_outside) const; 

	// evaluate BSDF limited to single scattering 
	// this is in average equivalent to eval(wi, wo, 1);
	virtual float evalSingleScattering(const Vector3f& wi, const Vector3f& wo) const; 

protected:
	float Fresnel(const Vector3f& wi, const Vector3f& wm, const float eta) const;
	Vector3f refract(const Vector3f &wi, const Vector3f &wm, const float eta) const;
};

/* Microsurface made of conductor material */
class MicrosurfaceDiffuse : public Microsurface
{
public:
	MicrosurfaceDiffuse(const bool height_uniform, // uniform or Gaussian
				 const bool slope_beckmann // Beckmann or GGX
)
		: Microsurface(height_uniform, slope_beckmann)
	{}

public:
	// evaluate local phase function 
	virtual float evalPhaseFunction(const Vector3f& wi, const Vector3f& wo) const;
	// sample local phase function
	virtual Vector3f samplePhaseFunction(const Vector3f& wi) const; 

	// evaluate BSDF limited to single scattering 
	// this is in average equivalent to eval(wi, wo, 1);
	virtual float evalSingleScattering(const Vector3f& wi, const Vector3f& wo) const; 
};


#endif

