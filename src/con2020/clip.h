#ifndef __CLIP_H__
#define __CLIP_H__
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#endif



template <typename T> T clip(T x, T mn, T mx) {
	return std::min(mx,std::max(x,mn));
}
