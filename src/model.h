#ifndef __MODEL_H__
#define __MODEL_H__
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include "internalfield.h"
#include "libcon2020.h"

namespace jupitermag {

void ModelField(double p0, double p1, double p2,
				const char *internal, const char *external,
				bool CartIn, bool CartOut,
				double *B0, double *B1, double *B2);

void ModelFieldArray(	int n, double *p0, double *p1, double *p2,
					const char *internal, const char *external,
					bool CartIn, bool CartOut,
					double *B0, double *B1, double *B2);

} /* namespace jupitermag */

extern "C" {
	void ModelField(double p0, double p1, double p2, 
					const char *internal, const char *external, 
					bool CartIn, bool CartOut,
					double *B0, double *B1, double *B2);

	void ModelFieldArray(	int n, double *p0, double *p1, double *p2, 
							const char *internal, const char *external, 
							bool CartIn, bool CartOut,
							double *B0, double *B1, double *B2);
}

#endif
