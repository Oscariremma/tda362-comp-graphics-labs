#include "material.h"
#include "sampling.h"
#include "labhelper.h"

using namespace labhelper;

namespace pathtracer
{
WiSample sampleHemisphereCosine(const vec3& wo, const vec3& n)
{
	mat3 tbn = tangentSpace(n);
	vec3 sample = cosineSampleHemisphere();
	WiSample r;
	r.wi = tbn * sample;
	if(dot(r.wi, n) > 0.0f)
		r.pdf = max(0.0f, dot(r.wi, n)) / M_PI;
	return r;
}

///////////////////////////////////////////////////////////////////////////
// A Lambertian (diffuse) material
///////////////////////////////////////////////////////////////////////////
vec3 Diffuse::f(const vec3& wi, const vec3& wo, const vec3& n) const
{
	vec3 N = n;
	if (dot(wo, N) < 0.0f) N = -N;

	if(dot(wi, N) <= 0.0f)
		return vec3(0.0f);
	if(!sameHemisphere(wi, wo, N))
		return vec3(0.0f);
	return (1.0f / M_PI) * color;
}

WiSample Diffuse::sample_wi(const vec3& wo, const vec3& n) const
{
	vec3 N = n;
	if (dot(wo, N) < 0.0f) N = -N;
	
	WiSample r = sampleHemisphereCosine(wo, N);
	r.f = f(r.wi, wo, n);
	return r;
}

vec3 MicrofacetBRDF::f(const vec3& wi, const vec3& wo, const vec3& n) const
{
	vec3 N = n;
	if(dot(wo, N) < 0.0f) N = -N;

	float cosThetaI = dot(N, wi);
	float cosThetaO = dot(N, wo);

	if(cosThetaI <= 0.0f || cosThetaO <= 0.0f)
		return vec3(0.0f);

	vec3 wh = normalize(wi + wo);
	float cosThetaH = max(0.0f, dot(N, wh));
	float woDotWh = max(0.0f, dot(wo, wh));
    
    // D(wh)
    float D = ((shininess + 2.0f) / (2.0f * M_PI)) * pow(cosThetaH, shininess);

    // G(wi, wo)
	float G = min(1.0f, min(2.0f * cosThetaH * cosThetaO / woDotWh,
	                        2.0f * cosThetaH * cosThetaI / woDotWh));
    
    // BRDF
    float brdf_val = (D * G) / (4.0f * cosThetaO * cosThetaI);
	return vec3(brdf_val);
}

WiSample MicrofacetBRDF::sample_wi(const vec3& wo, const vec3& n) const
{
	WiSample r;
	
	vec3 N = n;
	if(dot(wo, N) < 0.0f) N = -N;

	// Tangent space
	mat3 tbn = tangentSpace(N);
	vec3 tangent = tbn[0];
	vec3 bitangent = tbn[1];

	// Sample wh
	float phi = 2.0f * M_PI * randf();
	float cos_theta = pow(randf(), 1.0f / (shininess + 1.0f));
	float sin_theta = sqrt(max(0.0f, 1.0f - cos_theta * cos_theta));
	vec3 wh = normalize(sin_theta * cos(phi) * tangent + 
	                    sin_theta * sin(phi) * bitangent + 
	                    cos_theta * N);

	// Reflect wo around wh to get wi
	// wi = -wo + 2*dot(wo, wh)*wh
	r.wi = -wo + 2.0f * dot(wo, wh) * wh;
	r.wi = normalize(r.wi);

	// Calculate PDF
	// Calculate PDF
	float cosThetaH = max(0.0f, dot(N, wh));
	float pdf_wh = ((shininess + 1.0f) * pow(cosThetaH, shininess)) / (2.0f * M_PI);
	float woDotWh = max(0.0f, dot(wo, wh));
	r.pdf = pdf_wh / (4.0f * woDotWh);

	// Handle horizon check:
	// If sampled wi is below horizon, reflect it back up to preserve energy at grazing angles
	// (Common fix for microfacet darkening)
	if(dot(r.wi, N) <= 0.0f)
	{
		r.wi = reflect(r.wi, -N);
	}

	// Calculate f with the potentially new wi
	r.f = f(r.wi, wo, n);

	return r;
}


float BSDF::fresnel(const vec3& wi, const vec3& wo) const
{
	vec3 wh = normalize(wi + wo);
	float cosTheta = max(0.0f, dot(wh, wi));
	return R0 + (1.0f - R0) * pow(1.0f - cosTheta, 5.0f);
}


vec3 DielectricBSDF::f(const vec3& wi, const vec3& wo, const vec3& n) const
{
	float F;
	if (dot(wo, n) < 0.0f)
	{
		// Internal reflection
		// Recover IoR from R0
		float sqrtR0 = sqrt(clamp(R0, 0.0f, 0.99f));
		float ior = (1.0f + sqrtR0) / (1.0f - sqrtR0);
		
		float cosThetaI = dot(wo, -n);
		float sin2ThetaI = 1.0f - cosThetaI * cosThetaI;
		float sin2ThetaT = (ior * ior) * sin2ThetaI;
		
		if(sin2ThetaT > 1.0f)
		{
			// TIR
			F = 1.0f;
		}
		else
		{
			float cosThetaT = sqrt(1.0f - sin2ThetaT);
			F = R0 + (1.0f - R0) * pow(1.0f - cosThetaT, 5.0f);
		}
	}
	else
	{
		// External reflection
		F = fresnel(wi, wo);
	}

	return F * reflective_material->f(wi, wo, n) + (1.0f - F) * transmissive_material->f(wi, wo, n);
}

WiSample DielectricBSDF::sample_wi(const vec3& wo, const vec3& n) const
{
	WiSample r;
	
	float F;
	if (dot(wo, n) < 0.0f)
	{
		// Internal reflection
		// Recover IoR from R0
		float sqrtR0 = sqrt(clamp(R0, 0.0f, 0.99f));
		float ior = (1.0f + sqrtR0) / (1.0f - sqrtR0);
		
		float cosThetaI = dot(wo, -n);
		float sin2ThetaI = 1.0f - cosThetaI * cosThetaI;
		float sin2ThetaT = (ior * ior) * sin2ThetaI;
		
		if(sin2ThetaT > 1.0f)
		{
			// TIR
			F = 1.0f;
		}
		else
		{
			float cosThetaT = sqrt(1.0f - sin2ThetaT);
			F = R0 + (1.0f - R0) * pow(1.0f - cosThetaT, 5.0f);
		}
	}
	else
	{
		// External reflection sampling
		// Use ideal reflection vector to estimate F
		vec3 wi_reflect = reflect(-wo, n);
		F = fresnel(wi_reflect, wo);
	}

	if(randf() < F)
	{
		// Sample BRDF (Reflection)
		r = reflective_material->sample_wi(wo, n);
		r.pdf *= F;
		r.f *= F;
	}
	else
	{
		// Sample BTDF (Traismission/Diffuse)
		r = transmissive_material->sample_wi(wo, n);
		r.pdf *= (1.0f - F);
		r.f *= (1.0f - F);
	}
	
	return r;
}

vec3 MetalBSDF::f(const vec3& wi, const vec3& wo, const vec3& n) const
{
	float F = fresnel(wi, wo);
	return F * reflective_material->f(wi, wo, n) * color;
}

WiSample MetalBSDF::sample_wi(const vec3& wo, const vec3& n) const
{
	WiSample r = reflective_material->sample_wi(wo, n);
	float F = fresnel(r.wi, wo);
	r.f *= F * color;
	return r;
}


vec3 BSDFLinearBlend::f(const vec3& wi, const vec3& wo, const vec3& n) const
{
	return w * bsdf0->f(wi, wo, n) + (1.0f - w) * bsdf1->f(wi, wo, n);
}

WiSample BSDFLinearBlend::sample_wi(const vec3& wo, const vec3& n) const
{
	if(randf() < w)
	{
		return bsdf0->sample_wi(wo, n);
	}
	else
	{
		return bsdf1->sample_wi(wo, n);
	}
}



///////////////////////////////////////////////////////////////////////////
// A perfect specular refraction.
///////////////////////////////////////////////////////////////////////////
vec3 GlassBTDF::f(const vec3& wi, const vec3& wo, const vec3& n) const
{
	if(sameHemisphere(wi, wo, n))
	{
		return vec3(0.0f);
	}
	else
	{
		// Calculate ideal refraction vector wt
		float eta;
		vec3 N;
		if(dot(wo, n) > 0.0f)
		{
			N = n;
			eta = 1.0f / ior;
		}
		else
		{
			N = -n;
			eta = ior;
		}

		float w = dot(wo, N) * eta;
		float k = 1.0f + (w - eta) * (w + eta);
		
		if(k < 0.0f)
		{
			// TIR: No transmission
			return vec3(0.0f);
		}
		
		k = sqrt(k);
		vec3 wt = normalize(-eta * wo + (w - k) * N);
		
		// Evaluate Blinn-Phong-like lobe centered at wt
		// Using a high exponent to simulate sharp refraction of the point light
		float cosAlpha = max(0.0f, dot(wi, wt));
		float shininess = 10000.0f; // Very sharp to mimicking delta
		float normalization = (shininess + 2.0f) / (2.0f * M_PI);
		
		return vec3(1.0f) * normalization * pow(cosAlpha, shininess);
	}
}

WiSample GlassBTDF::sample_wi(const vec3& wo, const vec3& n) const
{
	WiSample r;

	float eta;
	glm::vec3 N;
	if(dot(wo, n) > 0.0f)
	{
		N = n;
		eta = 1.0f / ior;
	}
	else
	{
		N = -n;
		eta = ior;
	}

	// Alternatively:
	// d = dot(wo, N)
	// k = d * d (1 - eta*eta)
	// wi = normalize(-eta * wo + (d * eta - sqrt(k)) * N)

	// or

	// d = dot(n, wo)
	// k = 1 - eta*eta * (1 - d * d)
	// wi = - eta * wo + ( eta * d - sqrt(k) ) * N

	float w = dot(wo, N) * eta;
	float k = 1.0f + (w - eta) * (w + eta);
	if(k < 0.0f)
	{
		// Total internal reflection
		r.wi = reflect(-wo, n);
	}
	else
	{
		k = sqrt(k);
		r.wi = normalize(-eta * wo + (w - k) * N);
	}
	r.pdf = abs(dot(r.wi, n));
	r.f = vec3(1.0f, 1.0f, 1.0f);

	return r;
}

vec3 BTDFLinearBlend::f(const vec3& wi, const vec3& wo, const vec3& n) const
{
	return w * btdf0->f(wi, wo, n) + (1.0f - w) * btdf1->f(wi, wo, n);
}

WiSample BTDFLinearBlend::sample_wi(const vec3& wo, const vec3& n) const
{
	if(randf() < w)
	{
		WiSample r = btdf0->sample_wi(wo, n);
		return r;
	}
	else
	{
		WiSample r = btdf1->sample_wi(wo, n);
		return r;
	}
}


} // namespace pathtracer
