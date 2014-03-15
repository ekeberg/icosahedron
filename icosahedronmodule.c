#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>

PyDoc_STRVAR(icosahedron__doc__, "icosahedron(array_side, radius)\n\nGenerate an icosahedron. Radius is given in pixels and is defined as the distance to the corners.");
static PyObject *icosahedron(PyObject *self, PyObject *args, PyObject *kwargs)
{
  int image_side = 0.;
  double radius = 0.;

  if (!PyArg_ParseTuple(args, "id", &image_side, &radius)) {
    return NULL;
  }

  if (image_side <= 0) {
    PyErr_SetString(PyExc_ValueError, "Image side must be > 0.");
    return NULL;
  }

  if (radius <= 0.) {
    PyErr_SetString(PyExc_ValueError, "Radius must be > 0.");
  }

  int out_dim[] = {image_side, image_side, image_side};
  PyObject *out_array = (PyObject *)PyArray_FromDims(3, out_dim, NPY_FLOAT64);
  double *out = PyArray_DATA(out_array);

  double edge_thickness = 1./image_side;
  double half_edge_thickness = edge_thickness/2.;
  const double phi = (1.+sqrt(5.))/2.; //golden ratio

  double corner1[] = {0., 1., phi};
  double corner2[] = {1., phi, 0.};
  double corner3[] = {phi, 0., 1.};

  double original_radius = sqrt(1. + pow(phi, 2));
  double size_scaling = radius/original_radius;

  corner1[0] *= size_scaling; corner1[1] *= size_scaling; corner1[2] *= size_scaling;
  corner2[0] *= size_scaling; corner2[1] *= size_scaling; corner2[2] *= size_scaling;
  corner3[0] *= size_scaling; corner3[1] *= size_scaling; corner3[2] *= size_scaling;
  
  double center_z = (corner1[0]+corner2[0]+corner3[0])/3.;
  double center_y = (corner1[1]+corner2[1]+corner3[1])/3.;
  double center_x = (corner1[2]+corner2[2]+corner3[2])/3.;

  double normal1_z = (corner1[0] + corner2[0])/2. - center_z;
  double normal1_y = (corner1[1] + corner2[1])/2. - center_y;
  double normal1_x = (corner1[2] + corner2[2])/2. - center_x;

  double normal2_z = (corner2[0] + corner3[0])/2. - center_z;
  double normal2_y = (corner2[1] + corner3[1])/2. - center_y;
  double normal2_x = (corner2[2] + corner3[2])/2. - center_x;

  double normal3_z = (corner3[0] + corner1[0])/2. - center_z;
  double normal3_y = (corner3[1] + corner1[1])/2. - center_y;
  double normal3_x = (corner3[2] + corner1[2])/2. - center_x;

  double edge_distance = sqrt(pow(normal1_z, 2) + pow(normal1_y, 2) + pow(normal1_x, 2));
  normal1_z /= edge_distance; normal1_y /= edge_distance; normal1_x /= edge_distance;
  normal2_z /= edge_distance; normal2_y /= edge_distance; normal2_x /= edge_distance;
  normal3_z /= edge_distance; normal3_y /= edge_distance; normal3_x /= edge_distance;

  double face_normal_3[] = {phi/3., 0., (2.*phi+1.)/3.};
  double face_normal_2[] = {(2.*phi+1.)/3., phi/3., 0.};
  double face_normal_1[] = {0., (2.*phi+1.)/3., phi/3.};
  double face_normal_center[] = {1., 1., 1.};

  double face_distance = sqrt(pow(center_z, 2) + pow(center_y, 2) + pow(center_x, 2))/size_scaling;
  face_normal_1[0] /= face_distance; face_normal_1[1] /= face_distance; face_normal_1[2] /= face_distance;
  face_normal_2[0] /= face_distance; face_normal_2[1] /= face_distance; face_normal_2[2] /= face_distance;
  face_normal_3[0] /= face_distance; face_normal_3[1] /= face_distance; face_normal_3[2] /= face_distance;
  face_normal_center[0] /= sqrt(3.); face_normal_center[1] /= sqrt(3.);
  face_normal_center[2] /= sqrt(3.);
  
  face_distance = sqrt(pow(center_z, 2) + pow(center_y, 2) + pow(center_x, 2));

  double x, y, z;
  double projected_x, projected_y, projected_z;
  double scalar_product, distance;
  double image_side_float = (double) image_side;
  for (int x_pixel = 0; x_pixel < image_side; x_pixel++) {
    x = fabs(((double)x_pixel - image_side_float/2. + 0.5));///image_side_float*4.);
    for (int y_pixel = 0; y_pixel < image_side; y_pixel++) {
      y = fabs(((double)y_pixel - image_side_float/2. + 0.5));///image_side_float*4.);
      for (int z_pixel = 0; z_pixel < image_side; z_pixel++) {
	z = fabs(((double)z_pixel - image_side_float/2. + 0.5));///image_side_float*4.);
	scalar_product = x*face_normal_center[2] + y*face_normal_center[1] + z*face_normal_center[0];
	projected_x = x * face_distance/scalar_product;
	projected_y = y * face_distance/scalar_product;
	projected_z = z * face_distance/scalar_product;

	if ((projected_x-center_x)*normal1_x + (projected_y-center_y)*normal1_y +
	    (projected_z-center_z)*normal1_z > edge_distance) {
	  distance = x*face_normal_1[2] + y*face_normal_1[1] + z*face_normal_1[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.;
	  } else if (distance < face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.;
	  } else {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.-(face_distance-half_edge_thickness-distance);
	  }
	    
	} else if ((projected_x-center_x)*normal2_x + (projected_y-center_y)*normal2_y +
		   (projected_z-center_z)*normal2_z > edge_distance) {
	  distance = x*face_normal_2[2] + y*face_normal_2[1] + z*face_normal_2[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.;
	  } else if (distance < face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.;
	  } else {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.-(face_distance-half_edge_thickness-distance);
	  }

	} else if ((projected_x-center_x)*normal3_x + (projected_y-center_y)*normal3_y +
		   (projected_z-center_z)*normal3_z > edge_distance) {
	  distance = x*face_normal_3[2] + y*face_normal_3[1] + z*face_normal_3[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.;
	  } else if (distance < face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.;
	  } else {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.-(face_distance-half_edge_thickness-distance);
	  }

	} else {
	  distance = x*face_normal_center[2] + y*face_normal_center[1] + z*face_normal_center[0];
	  if (distance > face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 0.;
	  } else if (distance < face_distance + half_edge_thickness) {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.;
	  } else {
	    out[x_pixel*image_side*image_side + y_pixel*image_side + z_pixel] = 1.-(face_distance-half_edge_thickness-distance);
	  }
	}
      }
    }
  }
  return out_array;
}

static PyMethodDef IcosahedronMethods[] = {
  {"icosahedron", (PyCFunction)icosahedron, METH_VARARGS, icosahedron__doc__},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initicosahedron(void)
{
  import_array();
  PyObject *m = Py_InitModule3("icosahedron", IcosahedronMethods, "Create icosahedron density map");
  if (m == NULL)
    return;
}
