---
title:  "Triangle Rasterization"
date:   2022-12-13 00:00:00 +0200
categories: graphics
tags: graphics software-rendering
---

For some reason I needed a software rasterizer to generate images like this offline:

![vg-renderer colormap](/assets/img/colormap_vg.png)

---

The above image has been generated using [vg-renderer](https://github.com/jdryg/vg-renderer)'s indexedTriList() function using per-vertex colors, inside my custom immediate mode UI. I know this is not a proper way to render a colormap but this is what I have at the moment so I had to replicate this offline.

After a 5-minute search I found [tinyrenderer](https://github.com/ssloy/tinyrenderer) and in the [2nd lesson](https://github.com/ssloy/tinyrenderer/wiki/Lesson-2:-Triangle-rasterization-and-back-face-culling) they describe how to rasterize a triangle. In the 5th and final attempt they use barycentric coordinates to determine whether a pixel is inside the triangle or not. So I thought it would be easy to adapt in order to smoothly interpolate the 3 vertex colors inside the triangle.

```c
typedef struct swr_context
{
	uint32_t* m_FrameBuffer;
	uint32_t m_Width;
	uint32_t m_Height;
} swr_context;

static inline int32_t swr_mini(int32_t a, int32_t b)
{
	return a < b ? a : b;
}

static inline int32_t swr_maxi(int32_t a, int32_t b)
{
	return a > b ? a : b;
}

static inline int32_t swr_min3i(int32_t a, int32_t b, int32_t c)
{
	return swr_mini(a, swr_mini(b, c));
}

static inline int32_t swr_max3i(int32_t a, int32_t b, int32_t c)
{
	return swr_maxi(a, swr_maxi(b, c));
}

static bool swr_calcBarycentricCoords(int32_t x0, int32_t y0, int32_t x1, int32_t y1, int32_t x2, int32_t y2, int32_t x, int32_t y, float* bc)
{
	const float dx20 = (float)(x2 - x0);
	const float dx10 = (float)(x1 - x0);
	const float dy20 = (float)(y2 - y0);
	const float dy10 = (float)(y1 - y0);

	const float uz = dx20 * dy10 - dx10 * dy20;
	if (fabsf(uz) < 1.0f) {
		return false;
	}

	const float dx0p = (float)(x0 - x);
	const float dy0p = (float)(y0 - y);
	const float ux = dx10 * dy0p - dx0p * dy10;
	const float uy = dx0p * dy20 - dx20 * dy0p;

	bc[0] = 1.0f - ((ux + uy) / uz);
	bc[1] = uy / uz;
	bc[2] = ux / uz;

	return bc[0] >= 0.0f && bc[1] >= 0.0f && bc[2] >= 0.0f;
}

static void swrDrawTriangle_1(swr_context* ctx, int32_t x0, int32_t y0, int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color0, uint32_t color1, uint32_t color2)
{
	const int32_t bboxMinX = swr_maxi(swr_min3i(x0, x1, x2), 0);
	const int32_t bboxMinY = swr_maxi(swr_min3i(y0, y1, y2), 0);
	const int32_t bboxMaxX = swr_mini(swr_max3i(x0, x1, x2), (int32_t)ctx->m_Width - 1);
	const int32_t bboxMaxY = swr_mini(swr_max3i(y0, y1, y2), (int32_t)ctx->m_Height - 1);

	const uint32_t c0r = (color0 & SWR_COLOR_RED_Msk) >> SWR_COLOR_RED_Pos;
	const uint32_t c0g = (color0 & SWR_COLOR_GREEN_Msk) >> SWR_COLOR_GREEN_Pos;
	const uint32_t c0b = (color0 & SWR_COLOR_BLUE_Msk) >> SWR_COLOR_BLUE_Pos;
	const uint32_t c0a = (color0 & SWR_COLOR_ALPHA_Msk) >> SWR_COLOR_ALPHA_Pos;
	const uint32_t c1r = (color1 & SWR_COLOR_RED_Msk) >> SWR_COLOR_RED_Pos;
	const uint32_t c1g = (color1 & SWR_COLOR_GREEN_Msk) >> SWR_COLOR_GREEN_Pos;
	const uint32_t c1b = (color1 & SWR_COLOR_BLUE_Msk) >> SWR_COLOR_BLUE_Pos;
	const uint32_t c1a = (color1 & SWR_COLOR_ALPHA_Msk) >> SWR_COLOR_ALPHA_Pos;
	const uint32_t c2r = (color2 & SWR_COLOR_RED_Msk) >> SWR_COLOR_RED_Pos;
	const uint32_t c2g = (color2 & SWR_COLOR_GREEN_Msk) >> SWR_COLOR_GREEN_Pos;
	const uint32_t c2b = (color2 & SWR_COLOR_BLUE_Msk) >> SWR_COLOR_BLUE_Pos;
	const uint32_t c2a = (color2 & SWR_COLOR_ALPHA_Msk) >> SWR_COLOR_ALPHA_Pos;

	for (int32_t y = bboxMinY; y <= bboxMaxY; ++y) {
		const uint32_t y_width = (uint32_t)y * ctx->m_Width;

		for (int32_t x = bboxMinX; x <= bboxMaxX; ++x) {
			float barycentricCoords[3];
			if (!swr_calcBarycentricCoords(x0, y0, x1, y1, x2, y2, x, y, &barycentricCoords[0])) {
				continue;
			}

			const uint32_t cr = (uint32_t)(c0r * barycentricCoords[0] + c1r * barycentricCoords[1] + c2r * barycentricCoords[2]);
			const uint32_t cg = (uint32_t)(c0g * barycentricCoords[0] + c1g * barycentricCoords[1] + c2g * barycentricCoords[2]);
			const uint32_t cb = (uint32_t)(c0b * barycentricCoords[0] + c1b * barycentricCoords[1] + c2b * barycentricCoords[2]);
			const uint32_t ca = (uint32_t)(c0a * barycentricCoords[0] + c1a * barycentricCoords[1] + c2a * barycentricCoords[2]);

			swrDrawPixel(ctx, x, y, SWR_COLOR(cr, cg, cb, ca));
		}
	}
}
```

`swrDrawTriangle_1()` takes around 9ms to render the colormap below on a 1024x1024 canvas.
The mesh consists of 115 vertices and 198 triangles. Note that the colors at each vertex are not the same as the original image above. They are random.

![Software Rasterized colormap](/assets/img/colormap_sw.png)

My machine is an i5-12400F with 32GB of DDR4 dual channel RAM. Code has been compiled with MSVC 2022 (17.4.2) in x64 Release mode.

### 1st attempt

My first optimization attempt was to SIMDify the per-pixel color calculations.

```c
static void swrDrawTriangle_2(swr_context* ctx, int32_t x0, int32_t y0, int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color0, uint32_t color1, uint32_t color2)
{
	const int32_t bboxMinX = swr_maxi(swr_min3i(x0, x1, x2), 0);
	const int32_t bboxMinY = swr_maxi(swr_min3i(y0, y1, y2), 0);
	const int32_t bboxMaxX = swr_mini(swr_max3i(x0, x1, x2), (int32_t)ctx->m_Width - 1);
	const int32_t bboxMaxY = swr_mini(swr_max3i(y0, y1, y2), (int32_t)ctx->m_Height - 1);

	const __m128i xmm_zero = _mm_setzero_si128();
	const __m128 xmm_c0 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color0), xmm_zero), xmm_zero));
	const __m128 xmm_c1 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color1), xmm_zero), xmm_zero));
	const __m128 xmm_c2 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color2), xmm_zero), xmm_zero));

	for (int32_t y = bboxMinY; y <= bboxMaxY; ++y) {
		const uint32_t y_width = (uint32_t)y * ctx->m_Width;

		for (int32_t x = bboxMinX; x <= bboxMaxX; ++x) {
			float barycentricCoords[3];
			if (!swr_calcBarycentricCoords(x0, y0, x1, y1, x2, y2, x, y, &barycentricCoords[0])) {
				continue;
			}

			__m128 xmm_c0_scaled = _mm_mul_ps(xmm_c0, _mm_load1_ps(&barycentricCoords[0]));
			__m128 xmm_c1_scaled = _mm_mul_ps(xmm_c1, _mm_load1_ps(&barycentricCoords[1]));
			__m128 xmm_c2_scaled = _mm_mul_ps(xmm_c2, _mm_load1_ps(&barycentricCoords[2]));
			__m128 xmm_c = _mm_add_ps(_mm_add_ps(xmm_c0_scaled, xmm_c1_scaled), xmm_c2_scaled);
			__m128i imm_c = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c), xmm_zero), xmm_zero);
			_mm_storeu_si32(&ctx->m_FrameBuffer[x + y_width], imm_c);
		}
	}
}
```

The code before the loop expands the 3 RGBA colors into XMM registers.
The code inside the loop scales each color with the corresponding barycentric coordinate, adds them all together and repacks them into a single RGBA uint32_t.

`swrDrawTriangle_2()` takes around 4.9ms for the same colormap.

### 2nd attempt

My second optimization attempt was to inline `swr_calcBarycentricCoords()` and try to move some of the pixel-independent code outside the inner loop.

```c
static void swrDrawTriangle_3(swr_context* ctx, int32_t x0, int32_t y0, int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color0, uint32_t color1, uint32_t color2)
{
	const int32_t dx20 = x2 - x0;
	const int32_t dx10 = x1 - x0;
	const int32_t dy20 = y2 - y0;
	const int32_t dy10 = y1 - y0;
	const int32_t iuz = dx20 * dy10 - dx10 * dy20;
	if (abs(iuz) < 1) {
		return;
	}

	const int32_t bboxMinX = swr_maxi(swr_min3i(x0, x1, x2), 0);
	const int32_t bboxMinY = swr_maxi(swr_min3i(y0, y1, y2), 0);
	const int32_t bboxMaxX = swr_mini(swr_max3i(x0, x1, x2), (int32_t)ctx->m_Width - 1);
	const int32_t bboxMaxY = swr_mini(swr_max3i(y0, y1, y2), (int32_t)ctx->m_Height - 1);

	const __m128i xmm_zero = _mm_setzero_si128();
	const __m128 xmm_c0 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color0), xmm_zero), xmm_zero));
	const __m128 xmm_c1 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color1), xmm_zero), xmm_zero));
	const __m128 xmm_c2 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color2), xmm_zero), xmm_zero));

	const float inv_uz = 1.0f / (float)iuz;

	for (int32_t y = bboxMinY; y <= bboxMaxY; ++y) {
		const int32_t dy0p = y0 - y;
		const int32_t dx10_dy0p = dx10 * dy0p;
		const int32_t dx20_dy0p = dx20 * dy0p;

		uint32_t* framebufferRow = &ctx->m_FrameBuffer[y * ctx->m_Width];
		for (int32_t x = bboxMinX; x <= bboxMaxX; ++x) {
			const int32_t dx0p = x0 - x;
			const int32_t iux = dx10_dy0p - dx0p * dy10;
			const int32_t iuy = dx0p * dy20 - dx20_dy0p;

			const float bcx = (float)iux * inv_uz;
			const float bcy = (float)iuy * inv_uz;
			const float bcz = 1.0f - (bcx + bcy);
			if (bcz < 0.0f || bcy < 0.0f || bcx < 0.0f) {
				continue;
			}

			const __m128 xmm_c0_scaled = _mm_mul_ps(xmm_c0, _mm_load1_ps(&bcz));
			const __m128 xmm_c1_scaled = _mm_mul_ps(xmm_c1, _mm_load1_ps(&bcy));
			const __m128 xmm_c2_scaled = _mm_mul_ps(xmm_c2, _mm_load1_ps(&bcx));
			const __m128 xmm_c = _mm_add_ps(_mm_add_ps(xmm_c0_scaled, xmm_c1_scaled), xmm_c2_scaled);
			const __m128i imm_c = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c), xmm_zero), xmm_zero);
			_mm_storeu_si32(&framebufferRow[x], imm_c);
		}
	}
}
```

First, the check for the degenerate triangle case is performed at the beginning of the function in integer coordinates. The `fabsf(uz) < 1.0f` test ends up being `iuz == 0`.
The second change was to move y-dependent variables into the y loop and leave only the final x-dependent calculations in the inner loop.

`swrDrawTriangle_3()` takes around 2.9ms for the same colormap.

### 3rd attempt

My third optimization attempt was to get rid of the conditional inside the x loop. Since we are only touching pixels inside the triangle's bounding box it's guaranteed that there will be a non-empty x-span for each y. So I tried to calculate the range of the x span for each y.

Each barycentric coordinate has an equation in the form of: 
```
bc = (A + x * B) / signed_area;
```
which should always be greater than or equal to 0.

Depending on the sign of the numerator and the denominator I ended up with the following truth table for the limit affected by the corresponding combination:

`signed_area`  | `A` | `B` | limit
:---:|:---:|:--:|-----
 `> 0` | `> 0` | `> 0` | `none`
 `> 0` | `> 0` | `< 0` | `xmax`
 `> 0` | `> 0` | `== 0` | `none`
 `> 0` | `< 0` | `> 0` | `xmin`
 `> 0` | `< 0` | `< 0` | `N/A`
 `> 0` | `< 0` | `== 0` | `N/A`
 `> 0` | `== 0` | `> 0` | `none`
 `> 0` | `== 0` | `< 0` | `xmax = 0`
 `> 0` | `== 0` | `== 0` | `none`
 `< 0` | `> 0` | `> 0` | `N/A`
 `< 0` | `> 0` | `< 0` | `xmin`
 `< 0` | `> 0` | `== 0` | `N/A`
 `< 0` | `< 0` | `> 0` | `xmax`
 `< 0` | `< 0` | `< 0` | `none`
 `< 0` | `< 0` | `== 0` | `none`
 `< 0` | `== 0` | `> 0` | `xmax = 0`
 `< 0` | `== 0` | `< 0` | `none`
 `< 0` | `== 0` | `== 0` | `none`

- `none` means that no limit is affected because the equation will always be greater than or equal to 0.
- `N/A` means that the specified combination should not come up due to the above guarantee (not empty x-span for each y)
- `xmin` means that the equation is initially (`x == bboxMinX`) negative and there is an `x` at which it changes sign and becomes positive.
- `xmax` means that the equation is initially (`x == bboxMinX`) positive and there is an `x` at which it changes sign and becomes negative.
- `xmax=0` means that the only value of `x` at which the above equation is greater than or equal to 0 is at `xmin == xmax == 0`.

After applying the above truth table for all 3 barycentric coordinates, rearranging the code and cleanup I ended up with the following function:

```c
static void swrDrawTriangle_4(swr_context* ctx, int32_t x0, int32_t y0, int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color0, uint32_t color1, uint32_t color2)
{
	const int32_t dx20 = x2 - x0;
	const int32_t dx10 = x1 - x0;
	const int32_t dy20 = y2 - y0;
	const int32_t dy10 = y1 - y0;
	const int32_t iarea = dx20 * dy10 - dx10 * dy20;
	if (iarea == 0) {
		return;
	}

	const int32_t bboxMinX = swr_maxi(swr_min3i(x0, x1, x2), 0);
	const int32_t bboxMinY = swr_maxi(swr_min3i(y0, y1, y2), 0);
	const int32_t bboxMaxX = swr_mini(swr_max3i(x0, x1, x2), (int32_t)ctx->m_Width - 1);
	const int32_t bboxMaxY = swr_mini(swr_max3i(y0, y1, y2), (int32_t)ctx->m_Height - 1);
	const int32_t bboxWidth = bboxMaxX - bboxMinX;
	const int32_t bboxHeight = bboxMaxY - bboxMinY;

	const __m128i xmm_zero = _mm_setzero_si128();
	const __m128 xmm_c0 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color0), xmm_zero), xmm_zero));
	const __m128 xmm_c1 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color1), xmm_zero), xmm_zero));
	const __m128 xmm_c2 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color2), xmm_zero), xmm_zero));

	const float inv_area = 1.0f / (float)iarea;

	const int32_t dy01 = -dy10;
	const int32_t dx01 = -dx10;
	const int32_t dy01_dy20 = dy01 + dy20;

	const bool iarea_dy01_samesign = (iarea > 0 && dy01 > 0) || (iarea < 0 && dy01 < 0);
	const bool iarea_dy20_samesign = (iarea > 0 && dy20 > 0) || (iarea < 0 && dy20 < 0);
	const bool iarea_dy01dy20_samesign = (iarea > 0 && dy01_dy20 > 0) || (iarea < 0 && dy01_dy20 < 0);

	int32_t ivx = ((x0 - bboxMinX) * dy01 - (y0 - bboxMinY) * dx01);
	int32_t ivy = ((x0 - bboxMinX) * dy20 - (y0 - bboxMinY) * dx20);
	int32_t ivz = iarea - ivx - ivy;

	uint32_t* framebufferRow = &ctx->m_FrameBuffer[bboxMinX + bboxMinY * ctx->m_Width];
	for (int32_t iy = 0; iy <= bboxHeight; ++iy) {
		int32_t ixmin = 0;
		int32_t ixmax = (uint32_t)bboxWidth;

		// Calculate ixmin and ixmax
		if (iarea_dy01_samesign) {
			ixmax = swr_mini(ixmax, (int32_t)floorf((float)ivx / (float)dy01));
		} else if (ivx != 0) {
			ixmin = swr_maxi(ixmin, (int32_t)ceilf((float)ivx / (float)dy01));
		}

		if (iarea_dy20_samesign) {
			ixmax = swr_mini(ixmax, (int32_t)floorf((float)ivy / (float)dy20));
		} else if (ivy != 0) {
			ixmin = swr_maxi(ixmin, (int32_t)ceilf((float)ivy / (float)dy20));
		}

		if ((iarea > 0 && dy01_dy20 < 0 && ivz >= 0) || (iarea < 0 && dy01_dy20 > 0 && ivz <= 0)) {
			ixmax = swr_mini(ixmax, -(int32_t)ceilf((float)ivz / (float)dy01_dy20));
		} else if ((iarea > 0 && dy01_dy20 > 0 && ivz < 0) || (iarea < 0 && dy01_dy20 < 0 && ivz > 0)) {
			ixmin = swr_maxi(ixmin, -(int32_t)floorf((float)ivz / (float)dy01_dy20));
		}

		int32_t iux = ivx - ixmin * dy01;
		int32_t iuy = ivy - ixmin * dy20;
		int32_t iuz = ivz + ixmin * dy01_dy20;
		for (int32_t ix = ixmin; ix <= ixmax; ++ix) {
			const float bcx = (float)iux * inv_area;
			const float bcy = (float)iuy * inv_area;
			const float bcz = (float)iuz * inv_area;
			assert(bcx >= 0.0f && bcy >= 0.0f && bcz >= 0.0f);

			const __m128 xmm_c0_scaled = _mm_mul_ps(xmm_c0, _mm_load1_ps(&bcz));
			const __m128 xmm_c1_scaled = _mm_mul_ps(xmm_c1, _mm_load1_ps(&bcy));
			const __m128 xmm_c2_scaled = _mm_mul_ps(xmm_c2, _mm_load1_ps(&bcx));
			const __m128 xmm_c = _mm_add_ps(_mm_add_ps(xmm_c0_scaled, xmm_c1_scaled), xmm_c2_scaled);
			const __m128i imm_c = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c), xmm_zero), xmm_zero);
			_mm_storeu_si32(&framebufferRow[ix], imm_c);

			iux -= dy01;
			iuy -= dy20;
			iuz += dy01_dy20;
		}

		ivx += dx01;
		ivy += dx20;
		ivz -= dx01 + dx20;
		framebufferRow += ctx->m_Width;
	}
}
```

It might seem like a big jump from the previous version but it's nothing more than moving constant calculations outside of inner loops and introducing a couple more variables (admittedly poorly named) for making things a bit more clear and symmetric.

`swrDrawTriangle_4()` takes around 1.7ms for the same colormap.

### 4th and final attempt

The last attempt was to avoid the cases were the signed area is negative by checking it at the beginning of the function and swapping points #1 and #2 before continuing with the rest of the code. 

Also, I SIMDified the `iu` counters. Unfortunately, I couldn't do the same of the `iv` counters because the code that uses them is scalar and requires reloading them on a stack-allocated array in every iteration (see the `#if 0` blocks in the code below).

```c
static void swrDrawTriangle(swr_context* ctx, int32_t x0, int32_t y0, int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color0, uint32_t color1, uint32_t color2)
{
	int32_t iarea = (x2 - x0) * (y1 - y0) - (x1 - x0) * (y2 - y0);
	if (iarea == 0) {
		// Degenerate triangle with 0 area.
		return;
	} else if (iarea < 0) {
		// Swap (x1, y1) <-> (x2, y2)
		{ int32_t tmp = x1; x1 = x2; x2 = tmp; }
		{ int32_t tmp = y1; y1 = y2; y2 = tmp; }
		{ uint32_t tmp = color1; color1 = color2; color2 = tmp; }
		iarea = -iarea;
	}

	const int32_t dx20 = x2 - x0;
	const int32_t dx10 = x1 - x0;
	const int32_t dy20 = y2 - y0;
	const int32_t dy10 = y1 - y0;

	const int32_t bboxMinX = swr_maxi(swr_min3i(x0, x1, x2), 0);
	const int32_t bboxMinY = swr_maxi(swr_min3i(y0, y1, y2), 0);
	const int32_t bboxMaxX = swr_mini(swr_max3i(x0, x1, x2), (int32_t)ctx->m_Width - 1);
	const int32_t bboxMaxY = swr_mini(swr_max3i(y0, y1, y2), (int32_t)ctx->m_Height - 1);
	const int32_t bboxWidth = bboxMaxX - bboxMinX;
	const int32_t bboxHeight = bboxMaxY - bboxMinY;

	const __m128i xmm_zero = _mm_setzero_si128();
	const __m128 xmm_c0 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color0), xmm_zero), xmm_zero));
	const __m128 xmm_c1 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color1), xmm_zero), xmm_zero));
	const __m128 xmm_c2 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color2), xmm_zero), xmm_zero));

	const int32_t dy01 = -dy10;
	const int32_t dx01 = -dx10;
	const int32_t dy01_dy20 = dy01 + dy20;

	const __m128 xmm_inv_area = _mm_set1_ps(1.0f / (float)iarea);

	const float inv_dy01 = 1.0f / (float)dy01;
	const float inv_dy20 = 1.0f / (float)dy20;
	const float inv_dy01_dy20 = 1.0f / (float)dy01_dy20;

	const int32_t dx0min = x0 - bboxMinX;
	const int32_t dy0min = y0 - bboxMinY;

#if 0
	const __m128i imm_div = _mm_set_epi32(0, -(dx01 + dx20), dx20, dx01);

	const int32_t iv0 = dx0min * dy01 - dy0min * dx01;
	const int32_t iv1 = dx0min * dy20 - dy0min * dx20;
	const int32_t iv2 = iarea - iv0 - iv1;
	__m128i imm_iv = _mm_set_epi32(0, iv2, iv1, iv0);
#else
	int32_t iv[4] = { 0 };
	iv[0] = dx0min * dy01 - dy0min * dx01;
	iv[1] = dx0min * dy20 - dy0min * dx20;
	iv[2] = iarea - iv[0] - iv[1];

	const int32_t div[4] = {
		dx01,
		dx20,
		-(dx01 + dx20),
		0
	};
#endif

	const __m128i imm_diu = _mm_set_epi32(0, dy01_dy20, -dy20, -dy01);

	uint32_t* framebufferRow = &ctx->m_FrameBuffer[bboxMinX + bboxMinY * ctx->m_Width];
	for (int32_t iy = 0; iy <= bboxHeight; ++iy) {
		int32_t ixmin = 0;
		int32_t ixmax = (uint32_t)bboxWidth;

		// Calculate ixmin and ixmax
		{
#if 0
			int32_t iv[4];
			_mm_storeu_si128((__m128i*)&iv[0], imm_iv);
#endif

			if (dy01 > 0) {
				ixmax = swr_mini(ixmax, (int32_t)floorf((float)iv[0] * inv_dy01));
			} else if (iv[0] != 0) {
				ixmin = swr_maxi(ixmin, (int32_t)ceilf((float)iv[0] * inv_dy01));
			}

			if (dy20 > 0) {
				ixmax = swr_mini(ixmax, (int32_t)floorf((float)iv[1] * inv_dy20));
			} else if (iv[1] != 0) {
				ixmin = swr_maxi(ixmin, (int32_t)ceilf((float)iv[1] * inv_dy20));
			}

			if (dy01_dy20 < 0 && iv[2] >= 0) {
				ixmax = swr_mini(ixmax, -(int32_t)ceilf((float)iv[2] * inv_dy01_dy20));
			} else if (dy01_dy20 > 0 && iv[2] < 0) {
				ixmin = swr_maxi(ixmin, -(int32_t)floorf((float)iv[2] * inv_dy01_dy20));
			}
		}

#if 0
		// TODO: Avoid _mm_set_epi32 (requires 32-bit integer multiply)
		__m128i imm_iu = _mm_add_epi32(imm_iv, _mm_set_epi32(0, ixmin * dy01_dy20, -ixmin * dy20, -ixmin * dy01));
#else
		__m128i imm_iu = _mm_set_epi32(0, iv[2] + ixmin * dy01_dy20, iv[1] - ixmin * dy20, iv[0] - ixmin * dy01);
#endif

		for (int32_t ix = ixmin; ix <= ixmax; ++ix) {
			const __m128 xmm_bc = _mm_mul_ps(_mm_cvtepi32_ps(imm_iu), xmm_inv_area);

			const __m128 xmm_bc0 = _mm_shuffle_ps(xmm_bc, xmm_bc, _MM_SHUFFLE(0, 0, 0, 0));
			const __m128 xmm_bc1 = _mm_shuffle_ps(xmm_bc, xmm_bc, _MM_SHUFFLE(1, 1, 1, 1));
			const __m128 xmm_bc2 = _mm_shuffle_ps(xmm_bc, xmm_bc, _MM_SHUFFLE(2, 2, 2, 2));

			const __m128 xmm_c0_scaled = _mm_mul_ps(xmm_c0, xmm_bc2);
			const __m128 xmm_c1_scaled = _mm_mul_ps(xmm_c1, xmm_bc1);
			const __m128 xmm_c2_scaled = _mm_mul_ps(xmm_c2, xmm_bc0);
			const __m128 xmm_c = _mm_add_ps(_mm_add_ps(xmm_c0_scaled, xmm_c1_scaled), xmm_c2_scaled);
			const __m128i imm_c = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c), xmm_zero), xmm_zero);
			_mm_storeu_si32(&framebufferRow[ix], imm_c);

			imm_iu = _mm_add_epi32(imm_iu, imm_diu);
		}

#if 0
		imm_iv = _mm_add_epi32(imm_iv, imm_div);
#else
		iv[0] += div[0];
		iv[1] += div[1];
		iv[2] += div[2];
#endif
		framebufferRow += ctx->m_Width;
	}
}
```

`swrDrawTriangle()` takes around 1.3ms for the same colormap.

That's all for now. I'll edit if I manage to squeeze some more perf out of this.