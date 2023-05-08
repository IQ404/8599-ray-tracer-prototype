/*****************************************************************//**
 * \file   Camera.h
 * \brief  The class for the camera
 * 
 * \author Xiaoyang Liu
 * \date   May 2023
 *********************************************************************/

#ifndef CAMERA_H
#define CAMERA_H

#include "RayTracingToolbox.h"

class Camera
{
	Point3D origin;			// where camera locates.
	Vector3D horizontal;	// for calculating the left-to-right offset of the endpoint on the viewport
	Vector3D vertical;		// for calculating the bottom-to-top offset of the endpoint on the viewport
	Point3D bottom_left;	// the bottom-left point on the viewpoint

public:

	// Constructors:

	Camera(const double vertical_fov, const double aspect_ratio)
	// the fov should be provided in degree
	{
		double theta = degrees_to_radians(vertical_fov);
		double half_height = std::tan(theta / 2.0);
		double viewport_height = half_height * 2.0;
		double viewport_width = viewport_height * aspect_ratio;		// viewport has the same aspect ratio as the image if the pixels on the display is square shaped.
		double focal_length = 1.0;		// this is the distance from the camera to the viewport (projection plane) which, currently, is the z = -1 plane.
		origin = Point3D{ 0.0,0.0,0.0 };	// where camera locates.
		horizontal = Vector3D{ viewport_width, 0.0,0.0 };		// for calculating the left-to-right offset of the endpoint on the viewport
		vertical = Vector3D{ 0.0,viewport_height,0.0 };		// for calculating the bottom-to-top offset of the endpoint on the viewport
		bottom_left = origin - Vector3D{ 0.0,0.0,focal_length } - (horizontal / 2.0) - (vertical / 2.0);		// the bottom-left point on the viewpoint
	}

	// Methods:

	Ray extract_ray(double horizontal_offset_factor, double vertical_offset_factor)
	{
		return Ray{ origin, bottom_left + horizontal_offset_factor * horizontal + vertical_offset_factor * vertical - origin };
	}
};

#endif // !CAMERA_H
