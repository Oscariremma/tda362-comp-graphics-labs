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

	const float epsilon = 0.001f;

	const vec3 wh = normalize(wi + wo);
	const float NdotWh = max(epsilon, dot(n, wh));
	const float NdotWi = max(epsilon, dot(n, wi));
	const float NdotWo = max(epsilon, dot(n, wo));
	const float WoDotWh = max(epsilon, dot(wo, wh));

	const float D = (shininess + 2.0f) / (2.0f * M_PI) * pow(NdotWh, shininess);
	const float G = min(1.0f,
	                       min((2.0f * NdotWh * NdotWo) / WoDotWh,
	                           (2.0f * NdotWh * NdotWi) / WoDotWh));

	return vec3(D * G / (4.0f * NdotWi * NdotWo));
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
	const float F = fresnel(wi, wo);
	const vec3 reflective_part = reflective_material->f(wi, wo, n);
	const vec3 transmissive_part = transmissive_material->f(wi, wo, n);
	return F * reflective_part + (1.0f - F) * transmissive_part;
}

WiSample DielectricBSDF::sample_wi(const vec3& wo, const vec3& n) const
{
	WiSample r;
	vec3 N = n;
	if(dot(wo, N) < 0.0f)
		N = -N;

	if(randf() < 0.5)
	{
		r = reflective_material->sample_wi(wo, n);
		r.pdf *= 0.5f;
		float F = fresnel(r.wi, wo);
		r.f *= F;
	}
	else
	{
		r = transmissive_material->sample_wi(wo, n);
		r.pdf *= 0.5f;
		// Use macro-surface Fresnel approximation for transmission weighting
		// since we don't have a valid reflection vector for this lobe in the helper
		vec3 wi_reflect = reflect(-wo, N);
		float F = fresnel(wi_reflect, wo);
		r.f *= (1 - F);
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
		return vec3(0);
	}
	else
	{
		return vec3(1);
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
