#include "Pathtracer.h"
#include <memory>
#include <iostream>
#include <map>
#include <algorithm>
#include "material.h"
#include "embree.h"
#include "sampling.h"
#include "labhelper.h"

using namespace std;
using namespace glm;
using namespace labhelper;

namespace pathtracer
{
///////////////////////////////////////////////////////////////////////////////
// Global variables
///////////////////////////////////////////////////////////////////////////////
Settings settings;
Environment environment;
Image rendered_image;
PointLight point_light;
std::vector<DiscLight> disc_lights;

///////////////////////////////////////////////////////////////////////////
// Restart rendering of image
///////////////////////////////////////////////////////////////////////////
void restart()
{
	// No need to clear image,
	rendered_image.number_of_samples = 0;
}

int getSampleCount()
{
	return std::max(rendered_image.number_of_samples - 1, 0);
}

///////////////////////////////////////////////////////////////////////////
// On window resize, window size is passed in, actual size of pathtraced
// image may be smaller (if we're subsampling for speed)
///////////////////////////////////////////////////////////////////////////
void resize(int w, int h)
{
	rendered_image.width = w / settings.subsampling;
	rendered_image.height = h / settings.subsampling;
	rendered_image.data.resize(rendered_image.width * rendered_image.height);
	restart();
}

///////////////////////////////////////////////////////////////////////////
/// Return the radiance from a certain direction wi from the environment
/// map.
///////////////////////////////////////////////////////////////////////////
vec3 Lenvironment(const vec3& wi)
{
	const float theta = acos(std::max(-1.0f, std::min(1.0f, wi.y)));
	float phi = atan(wi.z, wi.x);
	if(phi < 0.0f)
		phi = phi + 2.0f * M_PI;
	vec2 lookup = vec2(phi / (2.0 * M_PI), 1 - theta / M_PI);
	return environment.multiplier * environment.map.sample(lookup.x, lookup.y);
}

///////////////////////////////////////////////////////////////////////////
/// Calculate the radiance going from one point (r.hitPosition()) in one
/// direction (-r.d), through path tracing.
///////////////////////////////////////////////////////////////////////////
vec3 Li(Ray& primary_ray)
{
	vec3 L = vec3(0.0f);
	vec3 path_throughput = vec3(1.0f);
	Ray current_ray = primary_ray;

	for(int bounces = 0; bounces <= settings.max_bounces; bounces++)
	{
		///////////////////////////////////////////////////////////////////
		// Get the intersection information from the ray
		///////////////////////////////////////////////////////////////////
		Intersection hit = getIntersection(current_ray);

		///////////////////////////////////////////////////////////////////
		// Create a Material tree for evaluating brdfs and calculating
		// sample directions.
		///////////////////////////////////////////////////////////////////

		Diffuse diffuse(hit.material->m_color);
		GlassBTDF glass(hit.material->m_ior);
		BTDFLinearBlend transparency_blend(hit.material->m_transparency, &glass, &diffuse);
		MicrofacetBRDF microfacet(hit.material->m_shininess);
		
		// Calculate R0 from IoR for dielectric interface (air-dielectric)
		// R0 = ((n-1)/(n+1))^2
		float R0 = pow((hit.material->m_ior - 1.0f) / (hit.material->m_ior + 1.0f), 2.0f);
		DielectricBSDF dielectric(&microfacet, &transparency_blend, R0);
		
		MetalBSDF metal(&microfacet, hit.material->m_color, hit.material->m_fresnel);
		BSDFLinearBlend metal_blend(hit.material->m_metalness, &metal, &dielectric);
		BSDF& mat = metal_blend;

		///////////////////////////////////////////////////////////////////
		// Calculate Direct Illumination from light.
		///////////////////////////////////////////////////////////////////
		{
			const float distance_to_light = length(point_light.position - hit.position);
			const float falloff_factor = 1.0f / (distance_to_light * distance_to_light);
			vec3 Li = point_light.intensity_multiplier * point_light.color * falloff_factor;
			vec3 wi = normalize(point_light.position - hit.position);

			Ray shadow_ray(hit.position + hit.geometry_normal * EPSILON, wi);
			shadow_ray.tfar = distance_to_light;
			
			// If shadow ray hits something, we need to check if it's transparent
			// Iterate through transparent surfaces
			vec3 shadow_throughput(1.0f);
			bool visible = true;
			
			// Limit shadow transparency depth
			for(int shadow_bounces = 0; shadow_bounces < 10; ++shadow_bounces)
			{
				if(!intersect(shadow_ray))
				{
					// Reached light (or at least no occluder within distance)
					break; 
				}

				Intersection shadow_hit = getIntersection(shadow_ray);
				
				// Check if occluder is beyond light source?
				// intersect() respects tfar, so we don't need manual distance check strictly if set correctly.
				// However, intersect updates tfar to hit distance.
				
				if(shadow_hit.material->m_transparency > 0.0f)
				{
					// Transparent object
					shadow_throughput *= shadow_hit.material->m_color * shadow_hit.material->m_transparency;
					
					// Continue ray
					// Bias slightly forward from hit
					shadow_ray.o = shadow_hit.position + shadow_ray.d * EPSILON;
					shadow_ray.tfar = distance_to_light - length(shadow_hit.position - hit.position); // Approximate remaining distance
                    // Actually clearer to just reset tfar relative to new origin or just trust tfar logic if we track distance.
                    // Let's just set tfar to a large value or remaining distance.
                    float dist_so_far = length(shadow_hit.position - hit.position); 
                    // Better: Point light is at known position.
                    shadow_ray.tfar = length(point_light.position - shadow_ray.o);
				}
				else
				{
					// Opaque occluder
					visible = false;
					break;
				}
			}

			if(visible)
			{
				L += path_throughput * mat.f(wi, hit.wo, hit.shading_normal) * Li * shadow_throughput * std::abs(dot(wi, hit.shading_normal));
			}
		}

		// Add emitted radiance from intersection
		L += path_throughput * hit.material->m_emission;

		// Sample an incoming direction (and the brdf and pdf for that direction)
		WiSample sample = mat.sample_wi(hit.wo, hit.shading_normal);

		// If the pdf is too close to zero, it means that the current path is extremely
		// unlikely to exist, so we break to avoid numerical instability
		if(sample.pdf < EPSILON) return L;

		float cosineterm = std::abs(dot(sample.wi, hit.shading_normal));

		path_throughput = path_throughput * (sample.f * cosineterm) / sample.pdf;

		// If pathThroughput is zero there is no need to continue, as no more light comes from this path
		if(path_throughput == vec3(0.0f)) return L;

		// Create next ray on path (existing instance can't be reused)
		vec3 offset = hit.geometry_normal * EPSILON;
		if(dot(sample.wi, hit.geometry_normal) < 0.0f) offset = -offset;
		current_ray = Ray(hit.position + offset, sample.wi); // Bias the ray

		// Intersect the new ray and if there is no intersection just
		// add environment contribution and finish
		if(!intersect(current_ray))
		{
			L += path_throughput * Lenvironment(current_ray.d);
			return L;
		}

		// Otherwise, reiterate for the new intersection
	}
	return L;
}

///////////////////////////////////////////////////////////////////////////
/// Used to homogenize points transformed with projection matrices
///////////////////////////////////////////////////////////////////////////
inline static glm::vec3 homogenize(const glm::vec4& p)
{
	return glm::vec3(p * (1.f / p.w));
}

///////////////////////////////////////////////////////////////////////////
/// Trace one path per pixel and accumulate the result in an image
///////////////////////////////////////////////////////////////////////////
void tracePaths(const glm::mat4& V, const glm::mat4& P)
{
	// Stop here if we have as many samples as we want
	if((int(rendered_image.number_of_samples) > settings.max_paths_per_pixel)
	   && (settings.max_paths_per_pixel != 0))
	{
		return;
	}
	vec3 camera_pos = vec3(glm::inverse(V) * vec4(0.0f, 0.0f, 0.0f, 1.0f));
	// Trace one path per pixel (the omp parallel stuf magically distributes the
	// pathtracing on all cores of your CPU).
	int num_rays = 0;
	vector<vec4> local_image(rendered_image.width * rendered_image.height, vec4(0.0f));

#pragma omp parallel for
	for(int y = 0; y < rendered_image.height; y++)
	{
		for(int x = 0; x < rendered_image.width; x++)
		{
			vec3 color;
			Ray primaryRay;
			primaryRay.o = camera_pos;
			// Create a ray that starts in the camera position and points toward
			// the current pixel on a virtual screen.
			vec2 screenCoord = vec2((float(x) + randf()) / float(rendered_image.width),
			                        (float(y) + randf()) / float(rendered_image.height));
			// Calculate direction
			vec4 viewCoord = vec4(screenCoord.x * 2.0f - 1.0f, screenCoord.y * 2.0f - 1.0f, 1.0f, 1.0f);
			vec3 p = homogenize(inverse(P * V) * viewCoord);
			primaryRay.d = normalize(p - camera_pos);
			// Intersect ray with scene
			if(intersect(primaryRay))
			{
				// If it hit something, evaluate the radiance from that point
				color = Li(primaryRay);
			}
			else
			{
				// Otherwise evaluate environment
				color = Lenvironment(primaryRay.d);
			}
			// Accumulate the obtained radiance to the pixels color
			float n = float(rendered_image.number_of_samples);
			rendered_image.data[y * rendered_image.width + x] =
			    rendered_image.data[y * rendered_image.width + x] * (n / (n + 1.0f))
			    + (1.0f / (n + 1.0f)) * color;
		}
	}
	rendered_image.number_of_samples += 1;
}
}; // namespace pathtracer
