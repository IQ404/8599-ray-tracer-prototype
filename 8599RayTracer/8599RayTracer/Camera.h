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

	Camera(const Point3D& look_from, const Point3D& look_at, const Vector3D& up_direction, const double vertical_fov, const double aspect_ratio)
	// the fov should be provided in degree.
	// For the current implementation, up_direction can not be aligned on w (the inverse view direction).
	{
		double theta = degrees_to_radians(vertical_fov);
		double half_height = std::tan(theta / 2.0);
		double viewport_height = half_height * 2.0;
		double viewport_width = viewport_height * aspect_ratio;		// viewport has the same aspect ratio as the image if the pixels on the display is square shaped.
		
		Vector3D w = unit_vector(look_from - look_at);
		Vector3D u = unit_vector(cross(up_direction, w));
		Vector3D v = cross(w, u);

		double focal_length = 1.0;		// this is the distance from the camera to the viewport (projection plane) which, currently, is the z = -1 plane.
		origin = look_from;	// where camera locates.
		horizontal = viewport_width * u;		// for calculating the left-to-right offset of the endpoint on the viewport
		vertical = viewport_height * v;		// for calculating the bottom-to-top offset of the endpoint on the viewport
		bottom_left = origin - focal_length * w - (horizontal / 2.0) - (vertical / 2.0);		// the bottom-left point on the viewpoint
	}

	// Methods:

	Ray extract_ray(double horizontal_offset_factor, double vertical_offset_factor)
	{
		return Ray{ origin, bottom_left + horizontal_offset_factor * horizontal + vertical_offset_factor * vertical - origin };
	}
};

#endif // !CAMERA_H
