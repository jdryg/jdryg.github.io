---
title:  "Triangle Rasterization - Part 2: More optimizations"
date:   2022-12-17 00:00:00 +0200
categories: graphics
tags: graphics software-rendering
---

**WARNING**: Long post ahead with a lot of (repetitive) code. I wish I knew how to make foldable code blocks. 

I had some time to work on the triangle rasterization function again. I also found a great series of posts from Fabian Giesen on [Optimizing Software Occlusion Culling](https://fgiesen.wordpress.com/2013/02/17/optimizing-sw-occlusion-culling-index/) where he talks about triangle rasterization and barycentric coordinates in parts 6 - 8 (tbh, those are the only parts I read). After reading the articles I realized my code looked very similar to the first code snippet from part 8, albeit a lot more obscured. That's what happens when you don't fully understand the math and just look at the code :)

---

Anyway, the idea is the same:

```
for each row of pixels covered by the triangle
	calculate the barycentric coords at the start of the row
	for each row pixel
		if all barycentric coords are >= 0
			 draw the pixel
		move on to the next pixel on the row
```

The only thing missing from my last attempt from the previous post is the conditional in the inner loop for the sign of all the barycentric coords. In my case, since I calculated the `xmin` and `xmax` before entering the inner loop, this test is not needed. It is guaranteed that all pixels touched on each row will be inside the triangle. 

Here is the last attempt from the previous post (I'll call it `swrDrawTriangle_5()` to distinguish it from the newer versions I'll present below.

```c
static void swrDrawTriangle_5(swr_context* ctx, int32_t x0, int32_t y0, int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color0, uint32_t color1, uint32_t color2)
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

	const __m128i imm_div = _mm_set_epi32(0, -(dx01 + dx20), dx20, dx01);

	const int32_t iv0 = dx0min * dy01 - dy0min * dx01;
	const int32_t iv1 = dx0min * dy20 - dy0min * dx20;
	const int32_t iv2 = iarea - iv0 - iv1;
	__m128i imm_iv = _mm_set_epi32(0, iv2, iv1, iv0);

	const __m128i imm_diu = _mm_set_epi32(0, dy01_dy20, -dy20, -dy01);

	uint32_t* framebufferRow = &ctx->m_FrameBuffer[bboxMinX + bboxMinY * ctx->m_Width];
	for (int32_t iy = 0; iy <= bboxHeight; ++iy) {
		int32_t ixmin = 0;
		int32_t ixmax = (uint32_t)bboxWidth;

		// Calculate ixmin and ixmax
		{
			int32_t iv[4];
			_mm_storeu_si128((__m128i*)&iv[0], imm_iv);

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

		__m128i imm_iu = _mm_add_epi32(imm_iv, _mm_set_epi32(0, ixmin * dy01_dy20, -ixmin * dy20, -ixmin * dy01));

		for (int32_t ix = ixmin; ix <= ixmax; ++ix) {
			assert(_mm_movemask_ps(_mm_castsi128_ps(imm_iu)) == 0);

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

		imm_iv = _mm_add_epi32(imm_iv, imm_div);
		framebufferRow += ctx->m_Width;
	}
}
```

Before moving forward, I decided to add a 128-sample moving average for the frame time in my demo app to have a better understanding of the performance of various versions of the function. Here is a table showing the numbers for the 5 versions from the previous post:

Function | Average | Min | Max | Description
--|:-:|:-:|:-: |---
`swrDrawTriangle_1()` | 8.50ms | 8.43ms | 9.49ms | Original C implementation
`swrDrawTriangle_2()` | 4.76ms | 4.71ms | 4.98ms | SIMD color calculations
`swrDrawTriangle_3()` | 2.71ms | 2.65ms | 3.01ms | Moved constant calcs out of inner loops
`swrDrawTriangle_4()` | 1.59ms | 1.53ms | 1.81ms | Calculate `xmin` and `xmax` for which all coordinates will be greater than or equal to 0.
`swrDrawTriangle_5()` | 1.34ms | 1.28ms | 1.46ms | SIMD counters

### Ver. 6: Processing multiple pixels at once

The next step in Mr. Giesen's 8th post was to show how processing multiple pixels at once can be implemented. You might wonder why I haven't taken a step back and try to rewrite my code to be cleaner (based on the code from that post). The reason was that I was lazy. So instead, I decided to try and implement multiple pixel processing in my version of the code.

```c
static void swrDrawTriangle_6(swr_context* ctx, int32_t x0, int32_t y0, int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color0, uint32_t color1, uint32_t color2)
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

	const int32_t bboxMinX = swr_maxi(swr_min3i(x0, x1, x2), 0);
	const int32_t bboxMinY = swr_maxi(swr_min3i(y0, y1, y2), 0);
	const int32_t bboxMaxX = swr_mini(swr_max3i(x0, x1, x2), (int32_t)ctx->m_Width - 1);
	const int32_t bboxMaxY = swr_mini(swr_max3i(y0, y1, y2), (int32_t)ctx->m_Height - 1);
	const int32_t bboxWidth = bboxMaxX - bboxMinX;
	const int32_t bboxHeight = bboxMaxY - bboxMinY;

	const __m128i imm_zero = _mm_setzero_si128();
	const __m128 xmm_c0 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color0), imm_zero), imm_zero));
	const __m128 xmm_c1 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color1), imm_zero), imm_zero));
	const __m128 xmm_c2 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color2), imm_zero), imm_zero));
	const __m128 xmm_c2_c0 = _mm_sub_ps(xmm_c2, xmm_c0);
	const __m128 xmm_c1_c0 = _mm_sub_ps(xmm_c1, xmm_c0);

	const int32_t dy01 = y0 - y1;
	const int32_t dx01 = x0 - x1;
	const int32_t dx20 = x2 - x0;
	const int32_t dy20 = y2 - y0;
	const int32_t dy01_dy20 = dy01 + dy20;

	const __m128 xmm_zero = _mm_setzero_ps();
	const __m128 xmm_inv_area = _mm_set1_ps(1.0f / (float)iarea);

	// Barycentric coordinate deltas for the X direction
	const __m128i imm_x_duvw_ = _mm_set_epi32(0, dy01_dy20, -dy20, -dy01);
	const __m128 xmm_x_duvw_1 = _mm_mul_ps(_mm_cvtepi32_ps(imm_x_duvw_), xmm_inv_area);
	const __m128 xmm_x_duvw_2 = _mm_add_ps(xmm_x_duvw_1, xmm_x_duvw_1);
	const __m128 xmm_x_duvw_3 = _mm_add_ps(xmm_x_duvw_1, xmm_x_duvw_2);
	const __m128 xmm_x_duvw_4 = _mm_add_ps(xmm_x_duvw_2, xmm_x_duvw_2);

	// UV deltas for the 1st and 2nd pixel
	const __m128 xmm_x_duv0_duv1 = _mm_shuffle_ps(xmm_zero, xmm_x_duvw_1, _MM_SHUFFLE(1, 0, 1, 0));

	// UV deltas for the 3rd and 4th pixel
	const __m128 xmm_x_duv2_duv3 = _mm_shuffle_ps(xmm_x_duvw_2, xmm_x_duvw_3, _MM_SHUFFLE(1, 0, 1, 0));

	// UV deltas for the next set of pixels
	const __m128 xmm_x_duv4_duv4 = _mm_shuffle_ps(xmm_x_duvw_4, xmm_x_duvw_4, _MM_SHUFFLE(1, 0, 1, 0));

	// Barycentric coordinate deltas for the Y direction
	const __m128i imm_y_duvw_ = _mm_set_epi32(0, -(dx01 + dx20), dx20, dx01);

	// Calculate unnormalized barycentric coordinates of the bounding box min.
	const int32_t bboxMin_u = (x0 - bboxMinX) * dy01 - (y0 - bboxMinY) * dx01;
	const int32_t bboxMin_v = (x0 - bboxMinX) * dy20 - (y0 - bboxMinY) * dx20;
	const int32_t bboxMin_w = iarea - bboxMin_u - bboxMin_v;
	__m128i imm_row_uvw_ = _mm_set_epi32(0, bboxMin_w, bboxMin_v, bboxMin_u);

	// 
	const float inv_dy01_dy20 = 1.0f / (float)dy01_dy20;
	const float inv_dy20 = 1.0f / (float)dy20;
	const float inv_dy01 = 1.0f / (float)dy01;

	uint32_t* framebufferRow = &ctx->m_FrameBuffer[bboxMinX + bboxMinY * ctx->m_Width];
	for (int32_t iy = 0; iy <= bboxHeight; ++iy) {
		int32_t ixmin = 0;
		int32_t ixmax = (uint32_t)bboxWidth;

		// Calculate ixmin and ixmax
		{
			int32_t row_uvw_[4];
			_mm_storeu_si128((__m128i*)& row_uvw_[0], imm_row_uvw_);

			if (dy01 > 0) {
				ixmax = swr_mini(ixmax, (int32_t)floorf((float)row_uvw_[0] * inv_dy01));
			} else if (row_uvw_[0] != 0) {
				ixmin = swr_maxi(ixmin, (int32_t)ceilf((float)row_uvw_[0] * inv_dy01));
			}

			if (dy20 > 0) {
				ixmax = swr_mini(ixmax, (int32_t)floorf((float)row_uvw_[1] * inv_dy20));
			} else if (row_uvw_[1] != 0) {
				ixmin = swr_maxi(ixmin, (int32_t)ceilf((float)row_uvw_[1] * inv_dy20));
			}

			if (dy01_dy20 < 0 && row_uvw_[2] >= 0) {
				ixmax = swr_mini(ixmax, -(int32_t)ceilf((float)row_uvw_[2] * inv_dy01_dy20));
			} else if (dy01_dy20 > 0 && row_uvw_[2] < 0) {
				ixmin = swr_maxi(ixmin, -(int32_t)floorf((float)row_uvw_[2] * inv_dy01_dy20));
			}
		}

		if (ixmin <= ixmax) {
			// Calculate normalized barycentric coordinates at ixmin of the current row of pixels.
			const __m128i imm_p0uvw_ = _mm_add_epi32(imm_row_uvw_, _mm_mullo_epi32_SSE2(_mm_set1_epi32(ixmin), imm_x_duvw_));
			const __m128 xmm_p0uvw_ = _mm_mul_ps(_mm_cvtepi32_ps(imm_p0uvw_), xmm_inv_area);
			__m128 xmm_p0uvuv = _mm_shuffle_ps(xmm_p0uvw_, xmm_p0uvw_, _MM_SHUFFLE(1, 0, 1, 0));

			uint32_t* frameBuffer = &framebufferRow[ixmin];
			const uint32_t numPixels = (uint32_t)((ixmax - ixmin) + 1);
			const uint32_t numIter = numPixels >> 2; // 4 pixels per iteration
			for (uint32_t iIter = 0; iIter < numIter; ++iIter) {
				// Calculate barycentric coordinates for the 4 pixels.
				const __m128 xmm_p0uv_p1uv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv0_duv1); // Barycentric coordinates of 1st and 2nd pixels
				const __m128 xmm_p2uv_p3uv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv2_duv3); // Barycentric coordinates of 3rd and 4th pixels

				// Extract barycentric coordinates for each pixel
				const __m128 xmm_p0u = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_p0v = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(1, 1, 1, 1));
				const __m128 xmm_p1u = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(2, 2, 2, 2));
				const __m128 xmm_p1v = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(3, 3, 3, 3));
				const __m128 xmm_p2u = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_p2v = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(1, 1, 1, 1));
				const __m128 xmm_p3u = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(2, 2, 2, 2));
				const __m128 xmm_p3v = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(3, 3, 3, 3));

				// Calculate color of each pixel
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128 xmm_c_p1 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p1u), _mm_mul_ps(xmm_c1_c0, xmm_p1v)));
				const __m128 xmm_c_p2 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p2u), _mm_mul_ps(xmm_c1_c0, xmm_p2v)));
				const __m128 xmm_c_p3 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p3u), _mm_mul_ps(xmm_c1_c0, xmm_p3v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				const __m128i imm_c_p1 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p1), imm_zero), imm_zero);
				const __m128i imm_c_p2 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p2), imm_zero), imm_zero);
				const __m128i imm_c_p3 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p3), imm_zero), imm_zero);

				// Pack the 4 colors into a XMM registers and store into framebuffer
				const __m128i imm_c_p0011 = _mm_shuffle_si128(imm_c_p0, imm_c_p1, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128i imm_c_p2233 = _mm_shuffle_si128(imm_c_p2, imm_c_p3, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128i imm_c_p0123 = _mm_shuffle_si128(imm_c_p0011, imm_c_p2233, _MM_SHUFFLE(2, 0, 2, 0));
				_mm_storeu_si128((__m128i*)frameBuffer, imm_c_p0123);

				// Move on to the next set of pixels
				xmm_p0uvuv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv4_duv4);
				frameBuffer += 4;
			}

			// Handle the remainder of pixels for this row
			const uint32_t rem = numPixels & 3;
			switch (rem) {
			case 3: {
				const __m128 xmm_p0u = _mm_shuffle_ps(xmm_p0uvuv, xmm_p0uvuv, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_p0v = _mm_shuffle_ps(xmm_p0uvuv, xmm_p0uvuv, _MM_SHUFFLE(1, 1, 1, 1));
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_c_p0);
				xmm_p0uvuv = _mm_add_ps(xmm_p0uvuv, xmm_x_duvw_1);
				frameBuffer++;
			} // fallthrough
			case 2: {
				const __m128 xmm_p0u = _mm_shuffle_ps(xmm_p0uvuv, xmm_p0uvuv, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_p0v = _mm_shuffle_ps(xmm_p0uvuv, xmm_p0uvuv, _MM_SHUFFLE(1, 1, 1, 1));
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_c_p0);
				xmm_p0uvuv = _mm_add_ps(xmm_p0uvuv, xmm_x_duvw_1);
				frameBuffer++;
			} // fallthrough
			case 1: {
				const __m128 xmm_p0u = _mm_shuffle_ps(xmm_p0uvuv, xmm_p0uvuv, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_p0v = _mm_shuffle_ps(xmm_p0uvuv, xmm_p0uvuv, _MM_SHUFFLE(1, 1, 1, 1));
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_c_p0);
				xmm_p0uvuv = _mm_add_ps(xmm_p0uvuv, xmm_x_duvw_1);
			} // fallthrough
			case 0:
			default:
				break;
			}
		}

		// Move on to the next row of pixels.
		imm_row_uvw_ = _mm_add_epi32(imm_row_uvw_, imm_y_duvw_);
		framebufferRow += ctx->m_Width;
	}
}
```

This is the first version of the function which doesn't use the 3rd barycentric coordinate. Color calculations have been rearranged so that only the first two barycentric coords are needed. So instead of:

```
pixelColor = vertexColor[0] * bc[0] 
           + vertexColor[1] * bc[1] 
           + vertexColor[2] * bc[2];
```

the code calculates:

```
pixelColor = vertexColor[2] 
           + (vertexColor[0] - vertexColor[2]) * bc[0] 
           + (vertexColor[1] - vertexColor[2]) * bc[1];
```

Having to track only two barycentric coords in the inner loop means I can calculate all four pixels with 2 SIMD additions and the start of the next set of pixels with 1 addition.

Handling of the remainder pixels for each row is performed one pixel at a time with a fallthrough switch statement.

One last thing to note about this version is that the "guarantee" I mentioned several times in the past about the non-empty x-span for each row of pixels doesn't seem to be true. Since I am now iterating in sets of 4 pixels at a time, I realized that in some cases `xmax` ends up being less than `xmin`. So a test before entering the inner loop is needed to avoid accessing invalid pixels.

`swrDrawTriangle_6()` stats:
- Average: 1.17ms
- Min: 1.10ms
- Max: 1.23ms

### Ver. 7: Rearranging the inner loop

Version 6 has too many shuffles in the inner loop in order to extract the 2 barycentric coords of each of the 4 pixels. Even though the [codegen for the inner loop](https://gist.github.com/jdryg/9d2e36a2540c73a41dd4d4a6531a6e05#file-swrdrawtriangle_6_inner_loop-asm) looks clean and the compiler can keep all required quantities in XMM regs without touching the stack, it's still a lot of work.

My next attempt was to move the shuffles out of the inner loop and instead of incrementing the 1 packed counter for the start of each 4-pixel set (`xmm_p0uvuv`) I tried incrementing each of the 8 barycentric coords separately.

```c
static void swrDrawTriangle_7(swr_context* ctx, int32_t x0, int32_t y0, int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color0, uint32_t color1, uint32_t color2)
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

	const int32_t bboxMinX = swr_maxi(swr_min3i(x0, x1, x2), 0);
	const int32_t bboxMinY = swr_maxi(swr_min3i(y0, y1, y2), 0);
	const int32_t bboxMaxX = swr_mini(swr_max3i(x0, x1, x2), (int32_t)ctx->m_Width - 1);
	const int32_t bboxMaxY = swr_mini(swr_max3i(y0, y1, y2), (int32_t)ctx->m_Height - 1);
	const int32_t bboxWidth = bboxMaxX - bboxMinX;
	const int32_t bboxHeight = bboxMaxY - bboxMinY;

	const __m128i imm_zero = _mm_setzero_si128();
	const __m128 xmm_c0 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color0), imm_zero), imm_zero));
	const __m128 xmm_c1 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color1), imm_zero), imm_zero));
	const __m128 xmm_c2 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color2), imm_zero), imm_zero));
	const __m128 xmm_c2_c0 = _mm_sub_ps(xmm_c2, xmm_c0);
	const __m128 xmm_c1_c0 = _mm_sub_ps(xmm_c1, xmm_c0);

	const int32_t dy01 = y0 - y1;
	const int32_t dx01 = x0 - x1;
	const int32_t dx20 = x2 - x0;
	const int32_t dy20 = y2 - y0;
	const int32_t dy01_dy20 = dy01 + dy20;

	const __m128 xmm_zero = _mm_setzero_ps();
	const __m128 xmm_inv_area = _mm_set1_ps(1.0f / (float)iarea);

	// Barycentric coordinate deltas for the X direction
	const __m128i imm_x_duvw_ = _mm_set_epi32(0, dy01_dy20, -dy20, -dy01);
	const __m128 xmm_x_duvw_1 = _mm_mul_ps(_mm_cvtepi32_ps(imm_x_duvw_), xmm_inv_area);
	const __m128 xmm_x_duvw_2 = _mm_add_ps(xmm_x_duvw_1, xmm_x_duvw_1);
	const __m128 xmm_x_duvw_3 = _mm_add_ps(xmm_x_duvw_1, xmm_x_duvw_2);
	const __m128 xmm_x_duvw_4 = _mm_add_ps(xmm_x_duvw_2, xmm_x_duvw_2);

	// UV deltas for the 1st and 2nd pixel
	const __m128 xmm_x_duv0_duv1 = _mm_shuffle_ps(xmm_zero, xmm_x_duvw_1, _MM_SHUFFLE(1, 0, 1, 0));

	// UV deltas for the 3rd and 4th pixel
	const __m128 xmm_x_duv2_duv3 = _mm_shuffle_ps(xmm_x_duvw_2, xmm_x_duvw_3, _MM_SHUFFLE(1, 0, 1, 0));

	const __m128 xmm_x_du = _mm_shuffle_ps(xmm_x_duvw_1, xmm_x_duvw_1, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_x_dv = _mm_shuffle_ps(xmm_x_duvw_1, xmm_x_duvw_1, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_x_du4 = _mm_shuffle_ps(xmm_x_duvw_4, xmm_x_duvw_4, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_x_dv4 = _mm_shuffle_ps(xmm_x_duvw_4, xmm_x_duvw_4, _MM_SHUFFLE(1, 1, 1, 1));

	// Barycentric coordinate deltas for the Y direction
	const __m128i imm_y_duvw_ = _mm_set_epi32(0, -(dx01 + dx20), dx20, dx01);

	// Calculate unnormalized barycentric coordinates of the bounding box min.
	const int32_t bboxMin_u = (x0 - bboxMinX) * dy01 - (y0 - bboxMinY) * dx01;
	const int32_t bboxMin_v = (x0 - bboxMinX) * dy20 - (y0 - bboxMinY) * dx20;
	const int32_t bboxMin_w = iarea - bboxMin_u - bboxMin_v;
	__m128i imm_row_uvw_ = _mm_set_epi32(0, bboxMin_w, bboxMin_v, bboxMin_u);

	// 
	const float inv_dy01_dy20 = 1.0f / (float)dy01_dy20;
	const float inv_dy20 = 1.0f / (float)dy20;
	const float inv_dy01 = 1.0f / (float)dy01;

	uint32_t* framebufferRow = &ctx->m_FrameBuffer[bboxMinX + bboxMinY * ctx->m_Width];
	for (int32_t iy = 0; iy <= bboxHeight; ++iy) {
		int32_t ixmin = 0;
		int32_t ixmax = (uint32_t)bboxWidth;

		// Calculate ixmin and ixmax
		{
			int32_t row_uvw_[4];
			_mm_storeu_si128((__m128i*)& row_uvw_[0], imm_row_uvw_);

			if (dy01 > 0) {
				ixmax = swr_mini(ixmax, (int32_t)floorf((float)row_uvw_[0] * inv_dy01));
			} else if (row_uvw_[0] != 0) {
				ixmin = swr_maxi(ixmin, (int32_t)ceilf((float)row_uvw_[0] * inv_dy01));
			}

			if (dy20 > 0) {
				ixmax = swr_mini(ixmax, (int32_t)floorf((float)row_uvw_[1] * inv_dy20));
			} else if (row_uvw_[1] != 0) {
				ixmin = swr_maxi(ixmin, (int32_t)ceilf((float)row_uvw_[1] * inv_dy20));
			}

			if (dy01_dy20 < 0 && row_uvw_[2] >= 0) {
				ixmax = swr_mini(ixmax, -(int32_t)ceilf((float)row_uvw_[2] * inv_dy01_dy20));
			} else if (dy01_dy20 > 0 && row_uvw_[2] < 0) {
				ixmin = swr_maxi(ixmin, -(int32_t)floorf((float)row_uvw_[2] * inv_dy01_dy20));
			}
		}

		if (ixmin <= ixmax) {
			// Calculate normalized barycentric coordinates at ixmin of the current row of pixels.
			const __m128i imm_p0uvw_ = _mm_add_epi32(imm_row_uvw_, _mm_mullo_epi32_SSE2(_mm_set1_epi32(ixmin), imm_x_duvw_));
			const __m128 xmm_p0uvw_ = _mm_mul_ps(_mm_cvtepi32_ps(imm_p0uvw_), xmm_inv_area);
			const __m128 xmm_p0uvuv = _mm_shuffle_ps(xmm_p0uvw_, xmm_p0uvw_, _MM_SHUFFLE(1, 0, 1, 0));

			// Calculate barycentric coordinates for the 4 pixels.
			const __m128 xmm_p0uv_p1uv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv0_duv1); // Barycentric coordinates of 1st and 2nd pixels
			const __m128 xmm_p2uv_p3uv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv2_duv3); // Barycentric coordinates of 3rd and 4th pixels

			// Extract barycentric coordinates for each pixel
			__m128 xmm_p0u = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(0, 0, 0, 0));
			__m128 xmm_p0v = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(1, 1, 1, 1));
			__m128 xmm_p1u = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(2, 2, 2, 2));
			__m128 xmm_p1v = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(3, 3, 3, 3));
			__m128 xmm_p2u = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(0, 0, 0, 0));
			__m128 xmm_p2v = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(1, 1, 1, 1));
			__m128 xmm_p3u = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(2, 2, 2, 2));
			__m128 xmm_p3v = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(3, 3, 3, 3));

			uint32_t* frameBuffer = &framebufferRow[ixmin];
			const uint32_t numPixels = (uint32_t)((ixmax - ixmin) + 1);
			const uint32_t numIter = numPixels >> 2; // 4 pixels per iteration
			for (uint32_t iIter = 0; iIter < numIter; ++iIter) {
				// Calculate color of each pixel
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128 xmm_c_p1 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p1u), _mm_mul_ps(xmm_c1_c0, xmm_p1v)));
				const __m128 xmm_c_p2 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p2u), _mm_mul_ps(xmm_c1_c0, xmm_p2v)));
				const __m128 xmm_c_p3 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p3u), _mm_mul_ps(xmm_c1_c0, xmm_p3v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				const __m128i imm_c_p1 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p1), imm_zero), imm_zero);
				const __m128i imm_c_p2 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p2), imm_zero), imm_zero);
				const __m128i imm_c_p3 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p3), imm_zero), imm_zero);

				// Pack the 4 colors into a XMM registers and store into framebuffer
				const __m128i imm_c_p0011 = _mm_shuffle_si128(imm_c_p0, imm_c_p1, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128i imm_c_p2233 = _mm_shuffle_si128(imm_c_p2, imm_c_p3, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128i imm_c_p0123 = _mm_shuffle_si128(imm_c_p0011, imm_c_p2233, _MM_SHUFFLE(2, 0, 2, 0));
				_mm_storeu_si128((__m128i*)frameBuffer, imm_c_p0123);

				// Move on to the next set of pixels
				xmm_p0u = _mm_add_ps(xmm_p0u, xmm_x_du4);
				xmm_p1u = _mm_add_ps(xmm_p1u, xmm_x_du4);
				xmm_p2u = _mm_add_ps(xmm_p2u, xmm_x_du4);
				xmm_p3u = _mm_add_ps(xmm_p3u, xmm_x_du4);
				xmm_p0v = _mm_add_ps(xmm_p0v, xmm_x_dv4);
				xmm_p1v = _mm_add_ps(xmm_p1v, xmm_x_dv4);
				xmm_p2v = _mm_add_ps(xmm_p2v, xmm_x_dv4);
				xmm_p3v = _mm_add_ps(xmm_p3v, xmm_x_dv4);
				frameBuffer += 4;
			}

			// Handle the remainder of pixels for this row
			const uint32_t rem = numPixels & 3;
			switch (rem) {
			case 3: {
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_c_p0);
				xmm_p0u = _mm_add_ps(xmm_p0u, xmm_x_du);
				xmm_p0v = _mm_add_ps(xmm_p0v, xmm_x_dv);
				frameBuffer++;
			} // fallthrough
			case 2: {
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_c_p0);
				xmm_p0u = _mm_add_ps(xmm_p0u, xmm_x_du);
				xmm_p0v = _mm_add_ps(xmm_p0v, xmm_x_dv);
				frameBuffer++;
			} // fallthrough
			case 1: {
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_c_p0);
				xmm_p0u = _mm_add_ps(xmm_p0u, xmm_x_du);
				xmm_p0v = _mm_add_ps(xmm_p0v, xmm_x_dv);
			} // fallthrough
			case 0:
			default:
				break;
			}
		}

		// Move on to the next row of pixels.
		imm_row_uvw_ = _mm_add_epi32(imm_row_uvw_, imm_y_duvw_);
		framebufferRow += ctx->m_Width;
	}
}
```

Unfortunately, the code requires too many XMM regs so the compiler needs to use the stack for storing some values (see [codegen](https://gist.github.com/jdryg/9d2e36a2540c73a41dd4d4a6531a6e05#file-swrdrawtriangle_7_inner_loop-asm)). Stil, it's faster than v6.

`swrDrawTriangle_7()` stats:
- Average: 1.11ms
- Min: 1.06ms
- Max: 1.21ms

### Ver. 8: Replacing floorf() and ceilf()

Taking a step out of the inner loop, it bothered me that there were calls to `floorf()` and `ceilf()` in the outer loop. Instead of conditionally calling those CRT functions, I thought of calculating all 3 floors and ceils at once out of the conditionals, using SIMD. Unfortunately, `_mm_round_ps()` is available on SSE4.1 which was something I wanted to avoid.

After a bit of googling I found [this article](http://dss.stephanierct.com/DevBlog/?p=8) which shows SSE2 versions of packed `_mm_floor_ps()` and `_mm_ceil_ps()`.

```c
// http://dss.stephanierct.com/DevBlog/?p=8
static const float xmm_ones[] = { 1.0f, 1.0f, 1.0f, 1.0f };

static inline __m128 _mm_floor_ps_SSE2(__m128 x)
{
#if 1
	__m128 j = _mm_load_ps(&xmm_ones[0]);
#else
	__m128i v0 = _mm_setzero_si128();
	__m128i v1 = _mm_cmpeq_epi32(v0, v0);
	__m128i ji = _mm_srli_epi32(v1, 25);
	__m128 j = _mm_castsi128_ps(_mm_slli_epi32(ji, 23)); //create vector 1.0f
#endif
	__m128i i = _mm_cvttps_epi32(x);
	__m128 fi = _mm_cvtepi32_ps(i);
	__m128 igx = _mm_cmpgt_ps(fi, x);
	j = _mm_and_ps(igx, j);
	return _mm_sub_ps(fi, j);
}

static inline __m128 _mm_ceil_ps_SSE2(__m128 x)
{
#if 1
	__m128 j = _mm_load_ps(&xmm_ones[0]);
#else
	__m128i v0 = _mm_setzero_si128();
	__m128i v1 = _mm_cmpeq_epi32(v0, v0);
	__m128i ji = _mm_srli_epi32(v1, 25);
	__m128 j = _mm_castsi128_ps(_mm_slli_epi32(ji, 23)); //create vector 1.0f
#endif
	__m128i i = _mm_cvttps_epi32(x);
	__m128 fi = _mm_cvtepi32_ps(i);
	__m128 igx = _mm_cmplt_ps(fi, x);
	j = _mm_and_ps(igx, j);
	return _mm_add_ps(fi, j);
}

static void swrDrawTriangle_8(swr_context* ctx, int32_t x0, int32_t y0, int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color0, uint32_t color1, uint32_t color2)
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

	const int32_t bboxMinX = swr_maxi(swr_min3i(x0, x1, x2), 0);
	const int32_t bboxMinY = swr_maxi(swr_min3i(y0, y1, y2), 0);
	const int32_t bboxMaxX = swr_mini(swr_max3i(x0, x1, x2), (int32_t)ctx->m_Width - 1);
	const int32_t bboxMaxY = swr_mini(swr_max3i(y0, y1, y2), (int32_t)ctx->m_Height - 1);
	const int32_t bboxWidth = bboxMaxX - bboxMinX;
	const int32_t bboxHeight = bboxMaxY - bboxMinY;

	const __m128i imm_zero = _mm_setzero_si128();
	const __m128 xmm_c0 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color0), imm_zero), imm_zero));
	const __m128 xmm_c1 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color1), imm_zero), imm_zero));
	const __m128 xmm_c2 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color2), imm_zero), imm_zero));
	const __m128 xmm_c2_c0 = _mm_sub_ps(xmm_c2, xmm_c0);
	const __m128 xmm_c1_c0 = _mm_sub_ps(xmm_c1, xmm_c0);

	const int32_t dy01 = y0 - y1;
	const int32_t dx01 = x0 - x1;
	const int32_t dx20 = x2 - x0;
	const int32_t dy20 = y2 - y0;
	const int32_t dy01_dy20 = dy01 + dy20;

	const __m128 xmm_zero = _mm_setzero_ps();
	const __m128 xmm_inv_area = _mm_set1_ps(1.0f / (float)iarea);

	// Barycentric coordinate deltas for the X direction
	const __m128i imm_x_duvw_ = _mm_set_epi32(0, dy01_dy20, -dy20, -dy01);
	const __m128 xmm_x_duvw_1 = _mm_mul_ps(_mm_cvtepi32_ps(imm_x_duvw_), xmm_inv_area);
	const __m128 xmm_x_duvw_2 = _mm_add_ps(xmm_x_duvw_1, xmm_x_duvw_1);
	const __m128 xmm_x_duvw_3 = _mm_add_ps(xmm_x_duvw_1, xmm_x_duvw_2);
	const __m128 xmm_x_duvw_4 = _mm_add_ps(xmm_x_duvw_2, xmm_x_duvw_2);

	// UV deltas for the 1st and 2nd pixel
	const __m128 xmm_x_duv0_duv1 = _mm_shuffle_ps(xmm_zero, xmm_x_duvw_1, _MM_SHUFFLE(1, 0, 1, 0));

	// UV deltas for the 3rd and 4th pixel
	const __m128 xmm_x_duv2_duv3 = _mm_shuffle_ps(xmm_x_duvw_2, xmm_x_duvw_3, _MM_SHUFFLE(1, 0, 1, 0));

	const __m128 xmm_x_du = _mm_shuffle_ps(xmm_x_duvw_1, xmm_x_duvw_1, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_x_dv = _mm_shuffle_ps(xmm_x_duvw_1, xmm_x_duvw_1, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_x_du4 = _mm_shuffle_ps(xmm_x_duvw_4, xmm_x_duvw_4, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_x_dv4 = _mm_shuffle_ps(xmm_x_duvw_4, xmm_x_duvw_4, _MM_SHUFFLE(1, 1, 1, 1));

	// Barycentric coordinate deltas for the Y direction
	const __m128i imm_y_duvw_ = _mm_set_epi32(0, -(dx01 + dx20), dx20, dx01);

	// Calculate unnormalized barycentric coordinates of the bounding box min.
	const int32_t bboxMin_u = (x0 - bboxMinX) * dy01 - (y0 - bboxMinY) * dx01;
	const int32_t bboxMin_v = (x0 - bboxMinX) * dy20 - (y0 - bboxMinY) * dx20;
	const int32_t bboxMin_w = iarea - bboxMin_u - bboxMin_v;
	__m128i imm_row_uvw_ = _mm_set_epi32(0, bboxMin_w, bboxMin_v, bboxMin_u);

	// 
	const __m128 xmm_row_uvw_scale = _mm_set_ps(0.0f, 1.0f / (float)dy01_dy20, 1.0f / (float)dy20, 1.0f / (float)dy01);

	uint32_t* framebufferRow = &ctx->m_FrameBuffer[bboxMinX + bboxMinY * ctx->m_Width];
	for (int32_t iy = 0; iy <= bboxHeight; ++iy) {
		int32_t ixmin = 0;
		int32_t ixmax = (uint32_t)bboxWidth;

		// Calculate ixmin and ixmax
		{
			int32_t row_uvw_[4];
			_mm_storeu_si128((__m128i*) & row_uvw_[0], imm_row_uvw_);

			const __m128 xmm_row_uvw_ = _mm_mul_ps(_mm_cvtepi32_ps(imm_row_uvw_), xmm_row_uvw_scale);
			const __m128i imm_row_uvw_floor = _mm_cvtps_epi32(_mm_floor_ps_SSE2(xmm_row_uvw_));
			const __m128i imm_row_uvw_ceil = _mm_cvtps_epi32(_mm_ceil_ps_SSE2(xmm_row_uvw_));

			int32_t row_uvw_floor[4];
			_mm_storeu_si128((__m128i*) & row_uvw_floor[0], imm_row_uvw_floor);

			int32_t row_uvw_ceil[4];
			_mm_storeu_si128((__m128i*) & row_uvw_ceil[0], imm_row_uvw_ceil);

			if (dy01 > 0) {
				ixmax = swr_mini(ixmax, row_uvw_floor[0]);
			} else if (row_uvw_[0] != 0) {
				ixmin = swr_maxi(ixmin, row_uvw_ceil[0]);
			}

			if (dy20 > 0) {
				ixmax = swr_mini(ixmax, row_uvw_floor[1]);
			} else if (row_uvw_[1] != 0) {
				ixmin = swr_maxi(ixmin, row_uvw_ceil[1]);
			}

			if (dy01_dy20 < 0 && row_uvw_[2] >= 0) {
				ixmax = swr_mini(ixmax, -row_uvw_ceil[2]);
			} else if (dy01_dy20 > 0 && row_uvw_[2] < 0) {
				ixmin = swr_maxi(ixmin, -row_uvw_floor[2]);
			}
		}

		if (ixmin <= ixmax) {
			// Calculate normalized barycentric coordinates at ixmin of the current row of pixels.
			const __m128i imm_p0uvw_ = _mm_add_epi32(imm_row_uvw_, _mm_mullo_epi32_SSE2(_mm_set1_epi32(ixmin), imm_x_duvw_));
			const __m128 xmm_p0uvw_ = _mm_mul_ps(_mm_cvtepi32_ps(imm_p0uvw_), xmm_inv_area);
			const __m128 xmm_p0uvuv = _mm_shuffle_ps(xmm_p0uvw_, xmm_p0uvw_, _MM_SHUFFLE(1, 0, 1, 0));

			// Calculate barycentric coordinates for the 4 pixels.
			const __m128 xmm_p0uv_p1uv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv0_duv1); // Barycentric coordinates of 1st and 2nd pixels
			const __m128 xmm_p2uv_p3uv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv2_duv3); // Barycentric coordinates of 3rd and 4th pixels

			// Extract barycentric coordinates for each pixel
			__m128 xmm_p0u = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(0, 0, 0, 0));
			__m128 xmm_p0v = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(1, 1, 1, 1));
			__m128 xmm_p1u = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(2, 2, 2, 2));
			__m128 xmm_p1v = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(3, 3, 3, 3));
			__m128 xmm_p2u = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(0, 0, 0, 0));
			__m128 xmm_p2v = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(1, 1, 1, 1));
			__m128 xmm_p3u = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(2, 2, 2, 2));
			__m128 xmm_p3v = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(3, 3, 3, 3));

			uint32_t* frameBuffer = &framebufferRow[ixmin];
			const uint32_t numPixels = (uint32_t)((ixmax - ixmin) + 1);
			const uint32_t numIter = numPixels >> 2; // 4 pixels per iteration
			for (uint32_t iIter = 0; iIter < numIter; ++iIter) {
				// Calculate color of each pixel
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128 xmm_c_p1 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p1u), _mm_mul_ps(xmm_c1_c0, xmm_p1v)));
				const __m128 xmm_c_p2 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p2u), _mm_mul_ps(xmm_c1_c0, xmm_p2v)));
				const __m128 xmm_c_p3 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p3u), _mm_mul_ps(xmm_c1_c0, xmm_p3v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				const __m128i imm_c_p1 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p1), imm_zero), imm_zero);
				const __m128i imm_c_p2 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p2), imm_zero), imm_zero);
				const __m128i imm_c_p3 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p3), imm_zero), imm_zero);

				// Pack the 4 colors into a XMM registers and store into framebuffer
				const __m128i imm_c_p0011 = _mm_shuffle_si128(imm_c_p0, imm_c_p1, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128i imm_c_p2233 = _mm_shuffle_si128(imm_c_p2, imm_c_p3, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128i imm_c_p0123 = _mm_shuffle_si128(imm_c_p0011, imm_c_p2233, _MM_SHUFFLE(2, 0, 2, 0));
				_mm_storeu_si128((__m128i*)frameBuffer, imm_c_p0123);

				// Move on to the next set of pixels
				xmm_p0u = _mm_add_ps(xmm_p0u, xmm_x_du4);
				xmm_p1u = _mm_add_ps(xmm_p1u, xmm_x_du4);
				xmm_p2u = _mm_add_ps(xmm_p2u, xmm_x_du4);
				xmm_p3u = _mm_add_ps(xmm_p3u, xmm_x_du4);
				xmm_p0v = _mm_add_ps(xmm_p0v, xmm_x_dv4);
				xmm_p1v = _mm_add_ps(xmm_p1v, xmm_x_dv4);
				xmm_p2v = _mm_add_ps(xmm_p2v, xmm_x_dv4);
				xmm_p3v = _mm_add_ps(xmm_p3v, xmm_x_dv4);
				frameBuffer += 4;
			}

			// Handle the remainder of pixels for this row
			const uint32_t rem = numPixels & 3;
			switch (rem) {
			case 3: {
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_c_p0);
				xmm_p0u = _mm_add_ps(xmm_p0u, xmm_x_du);
				xmm_p0v = _mm_add_ps(xmm_p0v, xmm_x_dv);
				frameBuffer++;
			} // fallthrough
			case 2: {
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_c_p0);
				xmm_p0u = _mm_add_ps(xmm_p0u, xmm_x_du);
				xmm_p0v = _mm_add_ps(xmm_p0v, xmm_x_dv);
				frameBuffer++;
			} // fallthrough
			case 1: {
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_c_p0);
				xmm_p0u = _mm_add_ps(xmm_p0u, xmm_x_du);
				xmm_p0v = _mm_add_ps(xmm_p0v, xmm_x_dv);
			} // fallthrough
			case 0:
			default:
				break;
			}
		}

		// Move on to the next row of pixels.
		imm_row_uvw_ = _mm_add_epi32(imm_row_uvw_, imm_y_duvw_);
		framebufferRow += ctx->m_Width;
	}
}
```

The rest of the code should be exactly the same as ver 7.

`swrDrawTriangle_8()` stats:
- Average: 1.03ms
- Min: 0.98ms
- Max: 1.13ms

### Ver. 9: Better RGBA32F to RGBA8 packing

Returning to the inner loop, I realized the code I used for packing the 4x RGBA32F colors into 4x RGBA8 was wasteful. There is a faster way to do it if you actually take a look into how `_mm_packs_epi32()` and `_mm_packus_epi16()` work.

```c
static void swrDrawTriangle_9(swr_context* ctx, int32_t x0, int32_t y0, int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color0, uint32_t color1, uint32_t color2)
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

	const int32_t bboxMinX = swr_maxi(swr_min3i(x0, x1, x2), 0);
	const int32_t bboxMinY = swr_maxi(swr_min3i(y0, y1, y2), 0);
	const int32_t bboxMaxX = swr_mini(swr_max3i(x0, x1, x2), (int32_t)ctx->m_Width - 1);
	const int32_t bboxMaxY = swr_mini(swr_max3i(y0, y1, y2), (int32_t)ctx->m_Height - 1);
	const int32_t bboxWidth = bboxMaxX - bboxMinX;
	const int32_t bboxHeight = bboxMaxY - bboxMinY;

	const __m128i imm_zero = _mm_setzero_si128();
	const __m128 xmm_c0 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color0), imm_zero), imm_zero));
	const __m128 xmm_c1 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color1), imm_zero), imm_zero));
	const __m128 xmm_c2 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color2), imm_zero), imm_zero));
	const __m128 xmm_c2_c0 = _mm_sub_ps(xmm_c2, xmm_c0);
	const __m128 xmm_c1_c0 = _mm_sub_ps(xmm_c1, xmm_c0);

	const int32_t dy01 = y0 - y1;
	const int32_t dx01 = x0 - x1;
	const int32_t dx20 = x2 - x0;
	const int32_t dy20 = y2 - y0;
	const int32_t dy01_dy20 = dy01 + dy20;

	const __m128 xmm_zero = _mm_setzero_ps();
	const __m128 xmm_inv_area = _mm_set1_ps(1.0f / (float)iarea);

	// Barycentric coordinate deltas for the X direction
	const __m128i imm_x_duvw_ = _mm_set_epi32(0, dy01_dy20, -dy20, -dy01);
	const __m128 xmm_x_duvw_1 = _mm_mul_ps(_mm_cvtepi32_ps(imm_x_duvw_), xmm_inv_area);
	const __m128 xmm_x_duvw_2 = _mm_add_ps(xmm_x_duvw_1, xmm_x_duvw_1);
	const __m128 xmm_x_duvw_3 = _mm_add_ps(xmm_x_duvw_1, xmm_x_duvw_2);
	const __m128 xmm_x_duvw_4 = _mm_add_ps(xmm_x_duvw_2, xmm_x_duvw_2);

	// UV deltas for the 1st and 2nd pixel
	const __m128 xmm_x_duv0_duv1 = _mm_shuffle_ps(xmm_zero, xmm_x_duvw_1, _MM_SHUFFLE(1, 0, 1, 0));

	// UV deltas for the 3rd and 4th pixel
	const __m128 xmm_x_duv2_duv3 = _mm_shuffle_ps(xmm_x_duvw_2, xmm_x_duvw_3, _MM_SHUFFLE(1, 0, 1, 0));

	const __m128 xmm_x_du = _mm_shuffle_ps(xmm_x_duvw_1, xmm_x_duvw_1, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_x_dv = _mm_shuffle_ps(xmm_x_duvw_1, xmm_x_duvw_1, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_x_du4 = _mm_shuffle_ps(xmm_x_duvw_4, xmm_x_duvw_4, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_x_dv4 = _mm_shuffle_ps(xmm_x_duvw_4, xmm_x_duvw_4, _MM_SHUFFLE(1, 1, 1, 1));

	// Barycentric coordinate deltas for the Y direction
	const __m128i imm_y_duvw_ = _mm_set_epi32(0, -(dx01 + dx20), dx20, dx01);

	// Calculate unnormalized barycentric coordinates of the bounding box min.
	const int32_t bboxMin_u = (x0 - bboxMinX) * dy01 - (y0 - bboxMinY) * dx01;
	const int32_t bboxMin_v = (x0 - bboxMinX) * dy20 - (y0 - bboxMinY) * dx20;
	const int32_t bboxMin_w = iarea - bboxMin_u - bboxMin_v;
	__m128i imm_row_uvw_ = _mm_set_epi32(0, bboxMin_w, bboxMin_v, bboxMin_u);

	// 
	const __m128 xmm_row_uvw_scale = _mm_set_ps(0.0f, 1.0f / (float)dy01_dy20, 1.0f / (float)dy20, 1.0f / (float)dy01);

	uint32_t* framebufferRow = &ctx->m_FrameBuffer[bboxMinX + bboxMinY * ctx->m_Width];
	for (int32_t iy = 0; iy <= bboxHeight; ++iy) {
		int32_t ixmin = 0;
		int32_t ixmax = (uint32_t)bboxWidth;

		// Calculate ixmin and ixmax
		{
			int32_t row_uvw_[4];
			_mm_storeu_si128((__m128i*) & row_uvw_[0], imm_row_uvw_);

			const __m128 xmm_row_uvw_ = _mm_mul_ps(_mm_cvtepi32_ps(imm_row_uvw_), xmm_row_uvw_scale);
			const __m128i imm_row_uvw_floor = _mm_cvtps_epi32(_mm_floor_ps_SSE2(xmm_row_uvw_));
			const __m128i imm_row_uvw_ceil = _mm_cvtps_epi32(_mm_ceil_ps_SSE2(xmm_row_uvw_));

			int32_t row_uvw_floor[4];
			_mm_storeu_si128((__m128i*) & row_uvw_floor[0], imm_row_uvw_floor);

			int32_t row_uvw_ceil[4];
			_mm_storeu_si128((__m128i*) & row_uvw_ceil[0], imm_row_uvw_ceil);

			if (dy01 > 0) {
				ixmax = swr_mini(ixmax, row_uvw_floor[0]);
			} else if (row_uvw_[0] != 0) {
				ixmin = swr_maxi(ixmin, row_uvw_ceil[0]);
			}

			if (dy20 > 0) {
				ixmax = swr_mini(ixmax, row_uvw_floor[1]);
			} else if (row_uvw_[1] != 0) {
				ixmin = swr_maxi(ixmin, row_uvw_ceil[1]);
			}

			if (dy01_dy20 < 0 && row_uvw_[2] >= 0) {
				ixmax = swr_mini(ixmax, -row_uvw_ceil[2]);
			} else if (dy01_dy20 > 0 && row_uvw_[2] < 0) {
				ixmin = swr_maxi(ixmin, -row_uvw_floor[2]);
			}
		}

		if (ixmin <= ixmax) {
			// Calculate normalized barycentric coordinates at ixmin of the current row of pixels.
			const __m128i imm_p0uvw_ = _mm_add_epi32(imm_row_uvw_, _mm_mullo_epi32_SSE2(_mm_set1_epi32(ixmin), imm_x_duvw_));
			const __m128 xmm_p0uvw_ = _mm_mul_ps(_mm_cvtepi32_ps(imm_p0uvw_), xmm_inv_area);
			const __m128 xmm_p0uvuv = _mm_shuffle_ps(xmm_p0uvw_, xmm_p0uvw_, _MM_SHUFFLE(1, 0, 1, 0));

			// Calculate barycentric coordinates for the 4 pixels.
			const __m128 xmm_p0uv_p1uv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv0_duv1); // Barycentric coordinates of 1st and 2nd pixels
			const __m128 xmm_p2uv_p3uv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv2_duv3); // Barycentric coordinates of 3rd and 4th pixels

			// Extract barycentric coordinates for each pixel
			__m128 xmm_p0u = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(0, 0, 0, 0));
			__m128 xmm_p0v = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(1, 1, 1, 1));
			__m128 xmm_p1u = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(2, 2, 2, 2));
			__m128 xmm_p1v = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p0uv_p1uv, _MM_SHUFFLE(3, 3, 3, 3));
			__m128 xmm_p2u = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(0, 0, 0, 0));
			__m128 xmm_p2v = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(1, 1, 1, 1));
			__m128 xmm_p3u = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(2, 2, 2, 2));
			__m128 xmm_p3v = _mm_shuffle_ps(xmm_p2uv_p3uv, xmm_p2uv_p3uv, _MM_SHUFFLE(3, 3, 3, 3));

			uint32_t* frameBuffer = &framebufferRow[ixmin];
			const uint32_t numPixels = (uint32_t)((ixmax - ixmin) + 1);
			const uint32_t numIter = numPixels >> 2; // 4 pixels per iteration
			for (uint32_t iIter = 0; iIter < numIter; ++iIter) {
				// Calculate the color of each pixel
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128 xmm_c_p1 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p1u), _mm_mul_ps(xmm_c1_c0, xmm_p1v)));
				const __m128 xmm_c_p2 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p2u), _mm_mul_ps(xmm_c1_c0, xmm_p2v)));
				const __m128 xmm_c_p3 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p3u), _mm_mul_ps(xmm_c1_c0, xmm_p3v)));

				// Convert to uint32_t
				const __m128i imm_rgba_p0_u32 = _mm_cvtps_epi32(xmm_c_p0);
				const __m128i imm_rgba_p1_u32 = _mm_cvtps_epi32(xmm_c_p1);
				const __m128i imm_rgba_p2_u32 = _mm_cvtps_epi32(xmm_c_p2);
				const __m128i imm_rgba_p3_u32 = _mm_cvtps_epi32(xmm_c_p3);

				// Pack into uint16_t
				const __m128i imm_rgba_p01_u16 = _mm_packs_epi32(imm_rgba_p0_u32, imm_rgba_p1_u32);
				const __m128i imm_rgba_p23_u16 = _mm_packs_epi32(imm_rgba_p2_u32, imm_rgba_p3_u32);

				// Pack into uint8_t
				const __m128i imm_rgba_p0123_u8 = _mm_packus_epi16(imm_rgba_p01_u16, imm_rgba_p23_u16);

				// Store
				_mm_storeu_si128((__m128i*)frameBuffer, imm_rgba_p0123_u8);

				// Move on to the next set of pixels
				xmm_p0u = _mm_add_ps(xmm_p0u, xmm_x_du4);
				xmm_p1u = _mm_add_ps(xmm_p1u, xmm_x_du4);
				xmm_p2u = _mm_add_ps(xmm_p2u, xmm_x_du4);
				xmm_p3u = _mm_add_ps(xmm_p3u, xmm_x_du4);
				xmm_p0v = _mm_add_ps(xmm_p0v, xmm_x_dv4);
				xmm_p1v = _mm_add_ps(xmm_p1v, xmm_x_dv4);
				xmm_p2v = _mm_add_ps(xmm_p2v, xmm_x_dv4);
				xmm_p3v = _mm_add_ps(xmm_p3v, xmm_x_dv4);
				frameBuffer += 4;
			}

			// Handle the remainder of pixels for this row
			const uint32_t rem = numPixels & 3;
			switch (rem) {
			case 3: {
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_c_p0);
				xmm_p0u = _mm_add_ps(xmm_p0u, xmm_x_du);
				xmm_p0v = _mm_add_ps(xmm_p0v, xmm_x_dv);
				frameBuffer++;
			} // fallthrough
			case 2: {
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_c_p0);
				xmm_p0u = _mm_add_ps(xmm_p0u, xmm_x_du);
				xmm_p0v = _mm_add_ps(xmm_p0v, xmm_x_dv);
				frameBuffer++;
			} // fallthrough
			case 1: {
				const __m128 xmm_c_p0 = _mm_add_ps(xmm_c0, _mm_add_ps(_mm_mul_ps(xmm_c2_c0, xmm_p0u), _mm_mul_ps(xmm_c1_c0, xmm_p0v)));
				const __m128i imm_c_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_c_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_c_p0);
				xmm_p0u = _mm_add_ps(xmm_p0u, xmm_x_du);
				xmm_p0v = _mm_add_ps(xmm_p0v, xmm_x_dv);
			} // fallthrough
			case 0:
			default:
				break;
			}
		}

		// Move on to the next row of pixels.
		imm_row_uvw_ = _mm_add_epi32(imm_row_uvw_, imm_y_duvw_);
		framebufferRow += ctx->m_Width;
	}
}
```

`swrDrawTriangle_9()` stats:
- Average: 0.90ms
- Min: 0.86ms
- Max: 1.00ms

Finally, submillisecond frame times. Yay!

### Ver. 10: Final SSE2 version (a regression)

One of the things I still didn't like about the code was that it treated colors as RGBA tuples. In other words, I won't be able to easily adapt the function if I want to interpolate other attributes over the triangle.

A more natural way of treating them would have been as independent R, G, B and A "varyings". This requires rearranging the inner loop and instead of calculating:

```
rgba_p0 = rgba0 + (rgba2 - rgba0) * u_p0 + (rgba1 - rgba0) * v_p0;
rgba_p1 = rgba0 + (rgba2 - rgba0) * u_p1 + (rgba1 - rgba0) * v_p1;
rgba_p2 = rgba0 + (rgba2 - rgba0) * u_p2 + (rgba1 - rgba0) * v_p2;
rgba_p3 = rgba0 + (rgba2 - rgba0) * u_p3 + (rgba1 - rgba0) * v_p3;
```

it calculates:

```
r_p0123 = r0 + (r2 - r0) * u_p0123 + (r1 - r0) * v_p0123;
g_p0123 = g0 + (g2 - g0) * u_p0123 + (g1 - g0) * v_p0123;
b_p0123 = b0 + (b2 - b0) * u_p0123 + (b1 - b0) * v_p0123;
a_p0123 = a0 + (b2 - a0) * u_p0123 + (a1 - a0) * v_p0123;
```

This way, a new varying attribute can be added by simply adding one more such equation:

```
var_p0123 = var0 + (var2 - var0) * u_p0123 + (var1 - var0) * v_p0123;
```

The only issue with this approach in the case of colors is that packing them into RGBA8 ends up being slower (no `_mm_shuffle_epi8()`, it has to be done with SSE2). [This answer](https://stackoverflow.com/a/24599209) from StackOverflow shows how it can be done.

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

	const int32_t bboxMinX = swr_maxi(swr_min3i(x0, x1, x2), 0);
	const int32_t bboxMinY = swr_maxi(swr_min3i(y0, y1, y2), 0);
	const int32_t bboxMaxX = swr_mini(swr_max3i(x0, x1, x2), (int32_t)ctx->m_Width - 1);
	const int32_t bboxMaxY = swr_mini(swr_max3i(y0, y1, y2), (int32_t)ctx->m_Height - 1);
	const int32_t bboxWidth = bboxMaxX - bboxMinX;
	const int32_t bboxHeight = bboxMaxY - bboxMinY;

	const __m128i imm_zero = _mm_setzero_si128();
	const __m128 xmm_rgba0 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color0), imm_zero), imm_zero));
	const __m128 xmm_rgba1 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color1), imm_zero), imm_zero));
	const __m128 xmm_rgba2 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color2), imm_zero), imm_zero));
	const __m128 xmm_drgba20 = _mm_sub_ps(xmm_rgba2, xmm_rgba0);
	const __m128 xmm_drgba10 = _mm_sub_ps(xmm_rgba1, xmm_rgba0);

	const __m128 xmm_r0 = _mm_shuffle_ps(xmm_rgba0, xmm_rgba0, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_g0 = _mm_shuffle_ps(xmm_rgba0, xmm_rgba0, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_b0 = _mm_shuffle_ps(xmm_rgba0, xmm_rgba0, _MM_SHUFFLE(2, 2, 2, 2));
	const __m128 xmm_a0 = _mm_shuffle_ps(xmm_rgba0, xmm_rgba0, _MM_SHUFFLE(3, 3, 3, 3));
	const __m128 xmm_dr20 = _mm_shuffle_ps(xmm_drgba20, xmm_drgba20, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_dg20 = _mm_shuffle_ps(xmm_drgba20, xmm_drgba20, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_db20 = _mm_shuffle_ps(xmm_drgba20, xmm_drgba20, _MM_SHUFFLE(2, 2, 2, 2));
	const __m128 xmm_da20 = _mm_shuffle_ps(xmm_drgba20, xmm_drgba20, _MM_SHUFFLE(3, 3, 3, 3));
	const __m128 xmm_dr10 = _mm_shuffle_ps(xmm_drgba10, xmm_drgba10, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_dg10 = _mm_shuffle_ps(xmm_drgba10, xmm_drgba10, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_db10 = _mm_shuffle_ps(xmm_drgba10, xmm_drgba10, _MM_SHUFFLE(2, 2, 2, 2));
	const __m128 xmm_da10 = _mm_shuffle_ps(xmm_drgba10, xmm_drgba10, _MM_SHUFFLE(3, 3, 3, 3));

	const int32_t dy01 = y0 - y1;
	const int32_t dx01 = x0 - x1;
	const int32_t dx20 = x2 - x0;
	const int32_t dy20 = y2 - y0;
	const int32_t dy01_dy20 = dy01 + dy20;

	const __m128 xmm_zero = _mm_setzero_ps();
	const __m128 xmm_inv_area = _mm_set1_ps(1.0f / (float)iarea);

	// Barycentric coordinate deltas for the X direction
	const __m128i imm_x_duvw_ = _mm_set_epi32(0, dy01_dy20, -dy20, -dy01);
	const __m128 xmm_x_duvw_1 = _mm_mul_ps(_mm_cvtepi32_ps(imm_x_duvw_), xmm_inv_area);
	const __m128 xmm_x_duvw_2 = _mm_add_ps(xmm_x_duvw_1, xmm_x_duvw_1);
	const __m128 xmm_x_duvw_3 = _mm_add_ps(xmm_x_duvw_1, xmm_x_duvw_2);
	const __m128 xmm_x_duvw_4 = _mm_add_ps(xmm_x_duvw_2, xmm_x_duvw_2);

	// UV deltas for the 1st and 2nd pixel
	const __m128 xmm_x_duv0_duv1 = _mm_shuffle_ps(xmm_zero, xmm_x_duvw_1, _MM_SHUFFLE(1, 0, 1, 0));

	// UV deltas for the 3rd and 4th pixel
	const __m128 xmm_x_duv2_duv3 = _mm_shuffle_ps(xmm_x_duvw_2, xmm_x_duvw_3, _MM_SHUFFLE(1, 0, 1, 0));

	const __m128 xmm_x_du = _mm_shuffle_ps(xmm_x_duvw_1, xmm_x_duvw_1, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_x_dv = _mm_shuffle_ps(xmm_x_duvw_1, xmm_x_duvw_1, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_x_du4 = _mm_shuffle_ps(xmm_x_duvw_4, xmm_x_duvw_4, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_x_dv4 = _mm_shuffle_ps(xmm_x_duvw_4, xmm_x_duvw_4, _MM_SHUFFLE(1, 1, 1, 1));

	// Barycentric coordinate deltas for the Y direction
	const __m128i imm_y_duvw_ = _mm_set_epi32(0, -(dx01 + dx20), dx20, dx01);

	// Calculate unnormalized barycentric coordinates of the bounding box min.
	const int32_t bboxMin_u = (x0 - bboxMinX) * dy01 - (y0 - bboxMinY) * dx01;
	const int32_t bboxMin_v = (x0 - bboxMinX) * dy20 - (y0 - bboxMinY) * dx20;
	const int32_t bboxMin_w = iarea - bboxMin_u - bboxMin_v;
	__m128i imm_row_uvw_ = _mm_set_epi32(0, bboxMin_w, bboxMin_v, bboxMin_u);

	// 
	const __m128 xmm_row_uvw_scale = _mm_set_ps(0.0f, 1.0f / (float)dy01_dy20, 1.0f / (float)dy20, 1.0f / (float)dy01);

	uint32_t* framebufferRow = &ctx->m_FrameBuffer[bboxMinX + bboxMinY * ctx->m_Width];
	for (int32_t iy = 0; iy <= bboxHeight; ++iy) {
		int32_t ixmin = 0;
		int32_t ixmax = (uint32_t)bboxWidth;

		// Calculate ixmin and ixmax
		{
			int32_t row_uvw_[4];
			_mm_storeu_si128((__m128i*) & row_uvw_[0], imm_row_uvw_);

			const __m128 xmm_row_uvw_ = _mm_mul_ps(_mm_cvtepi32_ps(imm_row_uvw_), xmm_row_uvw_scale);
			const __m128i imm_row_uvw_floor = _mm_cvtps_epi32(_mm_floor_ps_SSE2(xmm_row_uvw_));
			const __m128i imm_row_uvw_ceil = _mm_cvtps_epi32(_mm_ceil_ps_SSE2(xmm_row_uvw_));

			int32_t row_uvw_floor[4];
			_mm_storeu_si128((__m128i*) & row_uvw_floor[0], imm_row_uvw_floor);

			int32_t row_uvw_ceil[4];
			_mm_storeu_si128((__m128i*) & row_uvw_ceil[0], imm_row_uvw_ceil);

			if (dy01 > 0) {
				ixmax = swr_mini(ixmax, row_uvw_floor[0]);
			} else if (row_uvw_[0] != 0) {
				ixmin = swr_maxi(ixmin, row_uvw_ceil[0]);
			}

			if (dy20 > 0) {
				ixmax = swr_mini(ixmax, row_uvw_floor[1]);
			} else if (row_uvw_[1] != 0) {
				ixmin = swr_maxi(ixmin, row_uvw_ceil[1]);
			}

			if (dy01_dy20 < 0 && row_uvw_[2] >= 0) {
				ixmax = swr_mini(ixmax, -row_uvw_ceil[2]);
			} else if (dy01_dy20 > 0 && row_uvw_[2] < 0) {
				ixmin = swr_maxi(ixmin, -row_uvw_floor[2]);
			}
		}

		if (ixmin <= ixmax) {
			// Calculate normalized barycentric coordinates at ixmin of the current row of pixels.
			const __m128i imm_p0uvw_ = _mm_add_epi32(imm_row_uvw_, _mm_mullo_epi32_SSE2(_mm_set1_epi32(ixmin), imm_x_duvw_));
			const __m128 xmm_p0uvw_ = _mm_mul_ps(_mm_cvtepi32_ps(imm_p0uvw_), xmm_inv_area);
			const __m128 xmm_p0uvuv = _mm_shuffle_ps(xmm_p0uvw_, xmm_p0uvw_, _MM_SHUFFLE(1, 0, 1, 0));

			// Calculate barycentric coordinates for the 4 pixels.
			const __m128 xmm_p0uv_p1uv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv0_duv1); // Barycentric coordinates of 1st and 2nd pixels
			const __m128 xmm_p2uv_p3uv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv2_duv3); // Barycentric coordinates of 3rd and 4th pixels

			// Extract barycentric coordinates for each pixel
			__m128 xmm_u0123 = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p2uv_p3uv, _MM_SHUFFLE(2, 0, 2, 0));
			__m128 xmm_v0123 = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p2uv_p3uv, _MM_SHUFFLE(3, 1, 3, 1));

			uint32_t* frameBuffer = &framebufferRow[ixmin];
			const uint32_t numPixels = (uint32_t)((ixmax - ixmin) + 1);
			const uint32_t numIter = numPixels >> 2; // 4 pixels per iteration
			for (uint32_t iIter = 0; iIter < numIter; ++iIter) {
				// Calculate the color of each pixel
				const __m128 xmm_r_p0123 = _mm_add_ps(xmm_r0, _mm_add_ps(_mm_mul_ps(xmm_dr20, xmm_u0123), _mm_mul_ps(xmm_dr10, xmm_v0123)));
				const __m128 xmm_g_p0123 = _mm_add_ps(xmm_g0, _mm_add_ps(_mm_mul_ps(xmm_dg20, xmm_u0123), _mm_mul_ps(xmm_dg10, xmm_v0123)));
				const __m128 xmm_b_p0123 = _mm_add_ps(xmm_b0, _mm_add_ps(_mm_mul_ps(xmm_db20, xmm_u0123), _mm_mul_ps(xmm_db10, xmm_v0123)));
				const __m128 xmm_a_p0123 = _mm_add_ps(xmm_a0, _mm_add_ps(_mm_mul_ps(xmm_da20, xmm_u0123), _mm_mul_ps(xmm_da10, xmm_v0123)));

				// Pack into uint8_t
				// (uint8_t){ r0, r1, r2, r3, g0, g1, g2, g3, b0, b1, b2, b3, a0, a1, a2, a3 }
				const __m128i imm_r0123_g0123_b0123_a0123_u8 = _mm_packus_epi16(
					_mm_packs_epi32(_mm_cvtps_epi32(xmm_r_p0123), _mm_cvtps_epi32(xmm_g_p0123)), 
					_mm_packs_epi32(_mm_cvtps_epi32(xmm_b_p0123), _mm_cvtps_epi32(xmm_a_p0123))
				);

				// https://stackoverflow.com/questions/24595003/permuting-bytes-inside-sse-m128i-register
				// _mm_shuffle_epi8() with SSE2
				__m128i mask = _mm_set_epi8(0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF);

				// (uint8_t){ r0, r2, g0, g2, b0, b2, a0, a2, r1, r3, g1, g3, b1, b3, a1, a3 }
				const __m128i imm_r02_g02_b02_a02_r13_g13_b13_a13_u8 = 
					_mm_packus_epi16(
						_mm_and_si128(imm_r0123_g0123_b0123_a0123_u8, mask), 
						_mm_srli_epi16(imm_r0123_g0123_b0123_a0123_u8, 8)
					);

				// (uint8_t){ r0, g0, b0, a0, r1, g1, b1, a1, r2, g2, b3, a2, r3, g3, b3, a3 }
				const __m128i imm_rgba_p0123_u8 = 
					_mm_packus_epi16(
						_mm_and_si128(imm_r02_g02_b02_a02_r13_g13_b13_a13_u8, mask), 
						_mm_srli_epi16(imm_r02_g02_b02_a02_r13_g13_b13_a13_u8, 8)
					);

				// Store
				_mm_storeu_si128((__m128i*)frameBuffer, imm_rgba_p0123_u8);

				// Move on to the next set of pixels
				xmm_u0123 = _mm_add_ps(xmm_u0123, xmm_x_du4);
				xmm_v0123 = _mm_add_ps(xmm_v0123, xmm_x_dv4);
				frameBuffer += 4;
			}

			// Handle the remainder of pixels for this row
			const uint32_t rem = numPixels & 3;
			switch (rem) {
			case 3: {
				const __m128 xmm_r_p0 = _mm_add_ss(xmm_r0, _mm_add_ss(_mm_mul_ss(xmm_dr20, xmm_u0123), _mm_mul_ss(xmm_dr10, xmm_v0123)));
				const __m128 xmm_g_p0 = _mm_add_ss(xmm_g0, _mm_add_ss(_mm_mul_ss(xmm_dg20, xmm_u0123), _mm_mul_ss(xmm_dg10, xmm_v0123)));
				const __m128 xmm_b_p0 = _mm_add_ss(xmm_b0, _mm_add_ss(_mm_mul_ss(xmm_db20, xmm_u0123), _mm_mul_ss(xmm_db10, xmm_v0123)));
				const __m128 xmm_a_p0 = _mm_add_ss(xmm_a0, _mm_add_ss(_mm_mul_ss(xmm_da20, xmm_u0123), _mm_mul_ss(xmm_da10, xmm_v0123)));
				const __m128 xmm_rrgg_p0 = _mm_shuffle_ps(xmm_r_p0, xmm_g_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_bbaa_p0 = _mm_shuffle_ps(xmm_b_p0, xmm_a_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_rgba_p0 = _mm_shuffle_ps(xmm_rrgg_p0, xmm_bbaa_p0, _MM_SHUFFLE(2, 0, 2, 0));
				const __m128i imm_rgba_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_rgba_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_rgba_p0);
				xmm_u0123 = _mm_add_ss(xmm_u0123, xmm_x_du);
				xmm_v0123 = _mm_add_ss(xmm_v0123, xmm_x_dv);
				frameBuffer++;
			} // fallthrough
			case 2: {
				const __m128 xmm_r_p0 = _mm_add_ss(xmm_r0, _mm_add_ss(_mm_mul_ss(xmm_dr20, xmm_u0123), _mm_mul_ss(xmm_dr10, xmm_v0123)));
				const __m128 xmm_g_p0 = _mm_add_ss(xmm_g0, _mm_add_ss(_mm_mul_ss(xmm_dg20, xmm_u0123), _mm_mul_ss(xmm_dg10, xmm_v0123)));
				const __m128 xmm_b_p0 = _mm_add_ss(xmm_b0, _mm_add_ss(_mm_mul_ss(xmm_db20, xmm_u0123), _mm_mul_ss(xmm_db10, xmm_v0123)));
				const __m128 xmm_a_p0 = _mm_add_ss(xmm_a0, _mm_add_ss(_mm_mul_ss(xmm_da20, xmm_u0123), _mm_mul_ss(xmm_da10, xmm_v0123)));
				const __m128 xmm_rrgg_p0 = _mm_shuffle_ps(xmm_r_p0, xmm_g_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_bbaa_p0 = _mm_shuffle_ps(xmm_b_p0, xmm_a_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_rgba_p0 = _mm_shuffle_ps(xmm_rrgg_p0, xmm_bbaa_p0, _MM_SHUFFLE(2, 0, 2, 0));
				const __m128i imm_rgba_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_rgba_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_rgba_p0);
				xmm_u0123 = _mm_add_ss(xmm_u0123, xmm_x_du);
				xmm_v0123 = _mm_add_ss(xmm_v0123, xmm_x_dv);
				frameBuffer++;
			} // fallthrough
			case 1: {
				const __m128 xmm_r_p0 = _mm_add_ss(xmm_r0, _mm_add_ss(_mm_mul_ss(xmm_dr20, xmm_u0123), _mm_mul_ss(xmm_dr10, xmm_v0123)));
				const __m128 xmm_g_p0 = _mm_add_ss(xmm_g0, _mm_add_ss(_mm_mul_ss(xmm_dg20, xmm_u0123), _mm_mul_ss(xmm_dg10, xmm_v0123)));
				const __m128 xmm_b_p0 = _mm_add_ss(xmm_b0, _mm_add_ss(_mm_mul_ss(xmm_db20, xmm_u0123), _mm_mul_ss(xmm_db10, xmm_v0123)));
				const __m128 xmm_a_p0 = _mm_add_ss(xmm_a0, _mm_add_ss(_mm_mul_ss(xmm_da20, xmm_u0123), _mm_mul_ss(xmm_da10, xmm_v0123)));
				const __m128 xmm_rrgg_p0 = _mm_shuffle_ps(xmm_r_p0, xmm_g_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_bbaa_p0 = _mm_shuffle_ps(xmm_b_p0, xmm_a_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_rgba_p0 = _mm_shuffle_ps(xmm_rrgg_p0, xmm_bbaa_p0, _MM_SHUFFLE(2, 0, 2, 0));
				const __m128i imm_rgba_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_rgba_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_rgba_p0);
				xmm_u0123 = _mm_add_ss(xmm_u0123, xmm_x_du);
				xmm_v0123 = _mm_add_ss(xmm_v0123, xmm_x_dv);
			} // fallthrough
			case 0:
			default:
				break;
			}
		}

		// Move on to the next row of pixels.
		imm_row_uvw_ = _mm_add_epi32(imm_row_uvw_, imm_y_duvw_);
		framebufferRow += ctx->m_Width;
	}
}
```

Handling of remainder pixels gets a bit more complicated but it's all scalar code, so probably not that bad. I still have to check if there is a better way to do it in this case.

`swrDrawTriangle()` stats:
- Average: 0.94ms
- Min: 0.90ms
- Max: 1.01ms

### Ver. 11: SSSE3

If SSSE3 is available, RGBA32F to RGBA8 packing from the last version can be done using `_mm_shuffle_epi8()`.

```c
static void swrDrawTriangle_SSSE3(swr_context* ctx, int32_t x0, int32_t y0, int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color0, uint32_t color1, uint32_t color2)
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

	const int32_t bboxMinX = swr_maxi(swr_min3i(x0, x1, x2), 0);
	const int32_t bboxMinY = swr_maxi(swr_min3i(y0, y1, y2), 0);
	const int32_t bboxMaxX = swr_mini(swr_max3i(x0, x1, x2), (int32_t)ctx->m_Width - 1);
	const int32_t bboxMaxY = swr_mini(swr_max3i(y0, y1, y2), (int32_t)ctx->m_Height - 1);
	const int32_t bboxWidth = bboxMaxX - bboxMinX;
	const int32_t bboxHeight = bboxMaxY - bboxMinY;

	const __m128i imm_zero = _mm_setzero_si128();
	const __m128 xmm_rgba0 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color0), imm_zero), imm_zero));
	const __m128 xmm_rgba1 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color1), imm_zero), imm_zero));
	const __m128 xmm_rgba2 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color2), imm_zero), imm_zero));
	const __m128 xmm_drgba20 = _mm_sub_ps(xmm_rgba2, xmm_rgba0);
	const __m128 xmm_drgba10 = _mm_sub_ps(xmm_rgba1, xmm_rgba0);

	const __m128 xmm_r0 = _mm_shuffle_ps(xmm_rgba0, xmm_rgba0, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_g0 = _mm_shuffle_ps(xmm_rgba0, xmm_rgba0, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_b0 = _mm_shuffle_ps(xmm_rgba0, xmm_rgba0, _MM_SHUFFLE(2, 2, 2, 2));
	const __m128 xmm_a0 = _mm_shuffle_ps(xmm_rgba0, xmm_rgba0, _MM_SHUFFLE(3, 3, 3, 3));
	const __m128 xmm_dr20 = _mm_shuffle_ps(xmm_drgba20, xmm_drgba20, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_dg20 = _mm_shuffle_ps(xmm_drgba20, xmm_drgba20, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_db20 = _mm_shuffle_ps(xmm_drgba20, xmm_drgba20, _MM_SHUFFLE(2, 2, 2, 2));
	const __m128 xmm_da20 = _mm_shuffle_ps(xmm_drgba20, xmm_drgba20, _MM_SHUFFLE(3, 3, 3, 3));
	const __m128 xmm_dr10 = _mm_shuffle_ps(xmm_drgba10, xmm_drgba10, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_dg10 = _mm_shuffle_ps(xmm_drgba10, xmm_drgba10, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_db10 = _mm_shuffle_ps(xmm_drgba10, xmm_drgba10, _MM_SHUFFLE(2, 2, 2, 2));
	const __m128 xmm_da10 = _mm_shuffle_ps(xmm_drgba10, xmm_drgba10, _MM_SHUFFLE(3, 3, 3, 3));

	const int32_t dy01 = y0 - y1;
	const int32_t dx01 = x0 - x1;
	const int32_t dx20 = x2 - x0;
	const int32_t dy20 = y2 - y0;
	const int32_t dy01_dy20 = dy01 + dy20;

	const __m128 xmm_zero = _mm_setzero_ps();
	const __m128 xmm_inv_area = _mm_set1_ps(1.0f / (float)iarea);

	// Barycentric coordinate deltas for the X direction
	const __m128i imm_x_duvw_ = _mm_set_epi32(0, dy01_dy20, -dy20, -dy01);
	const __m128 xmm_x_duvw_1 = _mm_mul_ps(_mm_cvtepi32_ps(imm_x_duvw_), xmm_inv_area);
	const __m128 xmm_x_duvw_2 = _mm_add_ps(xmm_x_duvw_1, xmm_x_duvw_1);
	const __m128 xmm_x_duvw_3 = _mm_add_ps(xmm_x_duvw_1, xmm_x_duvw_2);
	const __m128 xmm_x_duvw_4 = _mm_add_ps(xmm_x_duvw_2, xmm_x_duvw_2);

	// UV deltas for the 1st and 2nd pixel
	const __m128 xmm_x_duv0_duv1 = _mm_shuffle_ps(xmm_zero, xmm_x_duvw_1, _MM_SHUFFLE(1, 0, 1, 0));

	// UV deltas for the 3rd and 4th pixel
	const __m128 xmm_x_duv2_duv3 = _mm_shuffle_ps(xmm_x_duvw_2, xmm_x_duvw_3, _MM_SHUFFLE(1, 0, 1, 0));

	const __m128 xmm_x_du = _mm_shuffle_ps(xmm_x_duvw_1, xmm_x_duvw_1, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_x_dv = _mm_shuffle_ps(xmm_x_duvw_1, xmm_x_duvw_1, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_x_du4 = _mm_shuffle_ps(xmm_x_duvw_4, xmm_x_duvw_4, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_x_dv4 = _mm_shuffle_ps(xmm_x_duvw_4, xmm_x_duvw_4, _MM_SHUFFLE(1, 1, 1, 1));

	// Barycentric coordinate deltas for the Y direction
	const __m128i imm_y_duvw_ = _mm_set_epi32(0, -(dx01 + dx20), dx20, dx01);

	// Calculate unnormalized barycentric coordinates of the bounding box min.
	const int32_t bboxMin_u = (x0 - bboxMinX) * dy01 - (y0 - bboxMinY) * dx01;
	const int32_t bboxMin_v = (x0 - bboxMinX) * dy20 - (y0 - bboxMinY) * dx20;
	const int32_t bboxMin_w = iarea - bboxMin_u - bboxMin_v;
	__m128i imm_row_uvw_ = _mm_set_epi32(0, bboxMin_w, bboxMin_v, bboxMin_u);

	// 
	const __m128 xmm_row_uvw_scale = _mm_set_ps(0.0f, 1.0f / (float)dy01_dy20, 1.0f / (float)dy20, 1.0f / (float)dy01);

	uint32_t* framebufferRow = &ctx->m_FrameBuffer[bboxMinX + bboxMinY * ctx->m_Width];
	for (int32_t iy = 0; iy <= bboxHeight; ++iy) {
		int32_t ixmin = 0;
		int32_t ixmax = (uint32_t)bboxWidth;

		// Calculate ixmin and ixmax
		{
			int32_t row_uvw_[4];
			_mm_storeu_si128((__m128i*) & row_uvw_[0], imm_row_uvw_);

			const __m128 xmm_row_uvw_ = _mm_mul_ps(_mm_cvtepi32_ps(imm_row_uvw_), xmm_row_uvw_scale);
			const __m128i imm_row_uvw_floor = _mm_cvtps_epi32(_mm_floor_ps_SSE2(xmm_row_uvw_));
			const __m128i imm_row_uvw_ceil = _mm_cvtps_epi32(_mm_ceil_ps_SSE2(xmm_row_uvw_));

			int32_t row_uvw_floor[4];
			_mm_storeu_si128((__m128i*) & row_uvw_floor[0], imm_row_uvw_floor);

			int32_t row_uvw_ceil[4];
			_mm_storeu_si128((__m128i*) & row_uvw_ceil[0], imm_row_uvw_ceil);

			if (dy01 > 0) {
				ixmax = swr_mini(ixmax, row_uvw_floor[0]);
			} else if (row_uvw_[0] != 0) {
				ixmin = swr_maxi(ixmin, row_uvw_ceil[0]);
			}

			if (dy20 > 0) {
				ixmax = swr_mini(ixmax, row_uvw_floor[1]);
			} else if (row_uvw_[1] != 0) {
				ixmin = swr_maxi(ixmin, row_uvw_ceil[1]);
			}

			if (dy01_dy20 < 0 && row_uvw_[2] >= 0) {
				ixmax = swr_mini(ixmax, -row_uvw_ceil[2]);
			} else if (dy01_dy20 > 0 && row_uvw_[2] < 0) {
				ixmin = swr_maxi(ixmin, -row_uvw_floor[2]);
			}
		}

		if (ixmin <= ixmax) {
			// Calculate normalized barycentric coordinates at ixmin of the current row of pixels.
			const __m128i imm_p0uvw_ = _mm_add_epi32(imm_row_uvw_, _mm_mullo_epi32_SSE2(_mm_set1_epi32(ixmin), imm_x_duvw_));
			const __m128 xmm_p0uvw_ = _mm_mul_ps(_mm_cvtepi32_ps(imm_p0uvw_), xmm_inv_area);
			const __m128 xmm_p0uvuv = _mm_shuffle_ps(xmm_p0uvw_, xmm_p0uvw_, _MM_SHUFFLE(1, 0, 1, 0));

			// Calculate barycentric coordinates for the 4 pixels.
			const __m128 xmm_p0uv_p1uv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv0_duv1); // Barycentric coordinates of 1st and 2nd pixels
			const __m128 xmm_p2uv_p3uv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv2_duv3); // Barycentric coordinates of 3rd and 4th pixels

			// Extract barycentric coordinates for each pixel
			__m128 xmm_u0123 = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p2uv_p3uv, _MM_SHUFFLE(2, 0, 2, 0));
			__m128 xmm_v0123 = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p2uv_p3uv, _MM_SHUFFLE(3, 1, 3, 1));

			uint32_t* frameBuffer = &framebufferRow[ixmin];
			const uint32_t numPixels = (uint32_t)((ixmax - ixmin) + 1);
			const uint32_t numIter = numPixels >> 2; // 4 pixels per iteration
			for (uint32_t iIter = 0; iIter < numIter; ++iIter) {
				// Calculate the color of each pixel
				const __m128 xmm_r_p0123 = _mm_add_ps(xmm_r0, _mm_add_ps(_mm_mul_ps(xmm_dr20, xmm_u0123), _mm_mul_ps(xmm_dr10, xmm_v0123)));
				const __m128 xmm_g_p0123 = _mm_add_ps(xmm_g0, _mm_add_ps(_mm_mul_ps(xmm_dg20, xmm_u0123), _mm_mul_ps(xmm_dg10, xmm_v0123)));
				const __m128 xmm_b_p0123 = _mm_add_ps(xmm_b0, _mm_add_ps(_mm_mul_ps(xmm_db20, xmm_u0123), _mm_mul_ps(xmm_db10, xmm_v0123)));
				const __m128 xmm_a_p0123 = _mm_add_ps(xmm_a0, _mm_add_ps(_mm_mul_ps(xmm_da20, xmm_u0123), _mm_mul_ps(xmm_da10, xmm_v0123)));

				// Pack into uint8_t
				// (uint8_t){ r0, r1, r2, r3, g0, g1, g2, g3, b0, b1, b2, b3, a0, a1, a2, a3 }
				const __m128i imm_r0123_g0123_b0123_a0123_u8 = _mm_packus_epi16(
					_mm_packs_epi32(_mm_cvtps_epi32(xmm_r_p0123), _mm_cvtps_epi32(xmm_g_p0123)),
					_mm_packs_epi32(_mm_cvtps_epi32(xmm_b_p0123), _mm_cvtps_epi32(xmm_a_p0123))
				);

				// Shuffle into RGBA uint32_t
				const __m128i mask = _mm_set_epi8(15, 11, 7, 3, 14, 10, 6, 2, 13, 9, 5, 1, 12, 8, 4, 0);
				const __m128i imm_rgba_p0123_u8 = _mm_shuffle_epi8(imm_r0123_g0123_b0123_a0123_u8, mask);

				// Store
				_mm_storeu_si128((__m128i*)frameBuffer, imm_rgba_p0123_u8);

				// Move on to the next set of pixels
				xmm_u0123 = _mm_add_ps(xmm_u0123, xmm_x_du4);
				xmm_v0123 = _mm_add_ps(xmm_v0123, xmm_x_dv4);
				frameBuffer += 4;
			}

			// Handle the remainder of pixels for this row
			const uint32_t rem = numPixels & 3;
			switch (rem) {
			case 3: {
				const __m128 xmm_r_p0 = _mm_add_ss(xmm_r0, _mm_add_ss(_mm_mul_ss(xmm_dr20, xmm_u0123), _mm_mul_ss(xmm_dr10, xmm_v0123)));
				const __m128 xmm_g_p0 = _mm_add_ss(xmm_g0, _mm_add_ss(_mm_mul_ss(xmm_dg20, xmm_u0123), _mm_mul_ss(xmm_dg10, xmm_v0123)));
				const __m128 xmm_b_p0 = _mm_add_ss(xmm_b0, _mm_add_ss(_mm_mul_ss(xmm_db20, xmm_u0123), _mm_mul_ss(xmm_db10, xmm_v0123)));
				const __m128 xmm_a_p0 = _mm_add_ss(xmm_a0, _mm_add_ss(_mm_mul_ss(xmm_da20, xmm_u0123), _mm_mul_ss(xmm_da10, xmm_v0123)));
				const __m128 xmm_rrgg_p0 = _mm_shuffle_ps(xmm_r_p0, xmm_g_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_bbaa_p0 = _mm_shuffle_ps(xmm_b_p0, xmm_a_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_rgba_p0 = _mm_shuffle_ps(xmm_rrgg_p0, xmm_bbaa_p0, _MM_SHUFFLE(2, 0, 2, 0));
				const __m128i imm_rgba_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_rgba_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_rgba_p0);
				xmm_u0123 = _mm_add_ss(xmm_u0123, xmm_x_du);
				xmm_v0123 = _mm_add_ss(xmm_v0123, xmm_x_dv);
				frameBuffer++;
			} // fallthrough
			case 2: {
				const __m128 xmm_r_p0 = _mm_add_ss(xmm_r0, _mm_add_ss(_mm_mul_ss(xmm_dr20, xmm_u0123), _mm_mul_ss(xmm_dr10, xmm_v0123)));
				const __m128 xmm_g_p0 = _mm_add_ss(xmm_g0, _mm_add_ss(_mm_mul_ss(xmm_dg20, xmm_u0123), _mm_mul_ss(xmm_dg10, xmm_v0123)));
				const __m128 xmm_b_p0 = _mm_add_ss(xmm_b0, _mm_add_ss(_mm_mul_ss(xmm_db20, xmm_u0123), _mm_mul_ss(xmm_db10, xmm_v0123)));
				const __m128 xmm_a_p0 = _mm_add_ss(xmm_a0, _mm_add_ss(_mm_mul_ss(xmm_da20, xmm_u0123), _mm_mul_ss(xmm_da10, xmm_v0123)));
				const __m128 xmm_rrgg_p0 = _mm_shuffle_ps(xmm_r_p0, xmm_g_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_bbaa_p0 = _mm_shuffle_ps(xmm_b_p0, xmm_a_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_rgba_p0 = _mm_shuffle_ps(xmm_rrgg_p0, xmm_bbaa_p0, _MM_SHUFFLE(2, 0, 2, 0));
				const __m128i imm_rgba_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_rgba_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_rgba_p0);
				xmm_u0123 = _mm_add_ss(xmm_u0123, xmm_x_du);
				xmm_v0123 = _mm_add_ss(xmm_v0123, xmm_x_dv);
				frameBuffer++;
			} // fallthrough
			case 1: {
				const __m128 xmm_r_p0 = _mm_add_ss(xmm_r0, _mm_add_ss(_mm_mul_ss(xmm_dr20, xmm_u0123), _mm_mul_ss(xmm_dr10, xmm_v0123)));
				const __m128 xmm_g_p0 = _mm_add_ss(xmm_g0, _mm_add_ss(_mm_mul_ss(xmm_dg20, xmm_u0123), _mm_mul_ss(xmm_dg10, xmm_v0123)));
				const __m128 xmm_b_p0 = _mm_add_ss(xmm_b0, _mm_add_ss(_mm_mul_ss(xmm_db20, xmm_u0123), _mm_mul_ss(xmm_db10, xmm_v0123)));
				const __m128 xmm_a_p0 = _mm_add_ss(xmm_a0, _mm_add_ss(_mm_mul_ss(xmm_da20, xmm_u0123), _mm_mul_ss(xmm_da10, xmm_v0123)));
				const __m128 xmm_rrgg_p0 = _mm_shuffle_ps(xmm_r_p0, xmm_g_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_bbaa_p0 = _mm_shuffle_ps(xmm_b_p0, xmm_a_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_rgba_p0 = _mm_shuffle_ps(xmm_rrgg_p0, xmm_bbaa_p0, _MM_SHUFFLE(2, 0, 2, 0));
				const __m128i imm_rgba_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_rgba_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_rgba_p0);
				xmm_u0123 = _mm_add_ss(xmm_u0123, xmm_x_du);
				xmm_v0123 = _mm_add_ss(xmm_v0123, xmm_x_dv);
			} // fallthrough
			case 0:
			default:
				break;
			}
		}

		// Move on to the next row of pixels.
		imm_row_uvw_ = _mm_add_epi32(imm_row_uvw_, imm_y_duvw_);
		framebufferRow += ctx->m_Width;
	}
}
```

`swrDrawTriangle_SSSE3()` stats:
- Average: 0.85ms
- Min: 0.80ms
- Max: 0.91ms

### Ver. 12: SSE 4.1

If SSE 4.1 is available floors and ceils in the outer loop can be calculated using `_mm_round_ps()`. Also `_mm_mullo_epi32_SSE2()` can be replaced with `_mm_mullo_epi32()` SSE4.1 intrinsic.

```c
static void swrDrawTriangle_SSE41(swr_context* ctx, int32_t x0, int32_t y0, int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color0, uint32_t color1, uint32_t color2)
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

	const int32_t bboxMinX = swr_maxi(swr_min3i(x0, x1, x2), 0);
	const int32_t bboxMinY = swr_maxi(swr_min3i(y0, y1, y2), 0);
	const int32_t bboxMaxX = swr_mini(swr_max3i(x0, x1, x2), (int32_t)ctx->m_Width - 1);
	const int32_t bboxMaxY = swr_mini(swr_max3i(y0, y1, y2), (int32_t)ctx->m_Height - 1);
	const int32_t bboxWidth = bboxMaxX - bboxMinX;
	const int32_t bboxHeight = bboxMaxY - bboxMinY;

	const __m128i imm_zero = _mm_setzero_si128();
	const __m128 xmm_rgba0 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color0), imm_zero), imm_zero));
	const __m128 xmm_rgba1 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color1), imm_zero), imm_zero));
	const __m128 xmm_rgba2 = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_loadu_si32(&color2), imm_zero), imm_zero));
	const __m128 xmm_drgba20 = _mm_sub_ps(xmm_rgba2, xmm_rgba0);
	const __m128 xmm_drgba10 = _mm_sub_ps(xmm_rgba1, xmm_rgba0);

	const __m128 xmm_r0 = _mm_shuffle_ps(xmm_rgba0, xmm_rgba0, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_g0 = _mm_shuffle_ps(xmm_rgba0, xmm_rgba0, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_b0 = _mm_shuffle_ps(xmm_rgba0, xmm_rgba0, _MM_SHUFFLE(2, 2, 2, 2));
	const __m128 xmm_a0 = _mm_shuffle_ps(xmm_rgba0, xmm_rgba0, _MM_SHUFFLE(3, 3, 3, 3));
	const __m128 xmm_dr20 = _mm_shuffle_ps(xmm_drgba20, xmm_drgba20, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_dg20 = _mm_shuffle_ps(xmm_drgba20, xmm_drgba20, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_db20 = _mm_shuffle_ps(xmm_drgba20, xmm_drgba20, _MM_SHUFFLE(2, 2, 2, 2));
	const __m128 xmm_da20 = _mm_shuffle_ps(xmm_drgba20, xmm_drgba20, _MM_SHUFFLE(3, 3, 3, 3));
	const __m128 xmm_dr10 = _mm_shuffle_ps(xmm_drgba10, xmm_drgba10, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_dg10 = _mm_shuffle_ps(xmm_drgba10, xmm_drgba10, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_db10 = _mm_shuffle_ps(xmm_drgba10, xmm_drgba10, _MM_SHUFFLE(2, 2, 2, 2));
	const __m128 xmm_da10 = _mm_shuffle_ps(xmm_drgba10, xmm_drgba10, _MM_SHUFFLE(3, 3, 3, 3));

	const int32_t dy01 = y0 - y1;
	const int32_t dx01 = x0 - x1;
	const int32_t dx20 = x2 - x0;
	const int32_t dy20 = y2 - y0;
	const int32_t dy01_dy20 = dy01 + dy20;

	const __m128 xmm_zero = _mm_setzero_ps();
	const __m128 xmm_inv_area = _mm_set1_ps(1.0f / (float)iarea);

	// Barycentric coordinate deltas for the X direction
	const __m128i imm_x_duvw_ = _mm_set_epi32(0, dy01_dy20, -dy20, -dy01);
	const __m128 xmm_x_duvw_1 = _mm_mul_ps(_mm_cvtepi32_ps(imm_x_duvw_), xmm_inv_area);
	const __m128 xmm_x_duvw_2 = _mm_add_ps(xmm_x_duvw_1, xmm_x_duvw_1);
	const __m128 xmm_x_duvw_3 = _mm_add_ps(xmm_x_duvw_1, xmm_x_duvw_2);
	const __m128 xmm_x_duvw_4 = _mm_add_ps(xmm_x_duvw_2, xmm_x_duvw_2);

	// UV deltas for the 1st and 2nd pixel
	const __m128 xmm_x_duv0_duv1 = _mm_shuffle_ps(xmm_zero, xmm_x_duvw_1, _MM_SHUFFLE(1, 0, 1, 0));

	// UV deltas for the 3rd and 4th pixel
	const __m128 xmm_x_duv2_duv3 = _mm_shuffle_ps(xmm_x_duvw_2, xmm_x_duvw_3, _MM_SHUFFLE(1, 0, 1, 0));

	const __m128 xmm_x_du = _mm_shuffle_ps(xmm_x_duvw_1, xmm_x_duvw_1, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_x_dv = _mm_shuffle_ps(xmm_x_duvw_1, xmm_x_duvw_1, _MM_SHUFFLE(1, 1, 1, 1));
	const __m128 xmm_x_du4 = _mm_shuffle_ps(xmm_x_duvw_4, xmm_x_duvw_4, _MM_SHUFFLE(0, 0, 0, 0));
	const __m128 xmm_x_dv4 = _mm_shuffle_ps(xmm_x_duvw_4, xmm_x_duvw_4, _MM_SHUFFLE(1, 1, 1, 1));

	// Barycentric coordinate deltas for the Y direction
	const __m128i imm_y_duvw_ = _mm_set_epi32(0, -(dx01 + dx20), dx20, dx01);

	// Calculate unnormalized barycentric coordinates of the bounding box min.
	const int32_t bboxMin_u = (x0 - bboxMinX) * dy01 - (y0 - bboxMinY) * dx01;
	const int32_t bboxMin_v = (x0 - bboxMinX) * dy20 - (y0 - bboxMinY) * dx20;
	const int32_t bboxMin_w = iarea - bboxMin_u - bboxMin_v;
	__m128i imm_row_uvw_ = _mm_set_epi32(0, bboxMin_w, bboxMin_v, bboxMin_u);

	// 
	const __m128 xmm_row_uvw_scale = _mm_set_ps(0.0f, 1.0f / (float)dy01_dy20, 1.0f / (float)dy20, 1.0f / (float)dy01);

	uint32_t* framebufferRow = &ctx->m_FrameBuffer[bboxMinX + bboxMinY * ctx->m_Width];
	for (int32_t iy = 0; iy <= bboxHeight; ++iy) {
		int32_t ixmin = 0;
		int32_t ixmax = (uint32_t)bboxWidth;

		// Calculate ixmin and ixmax
		{
			int32_t row_uvw_[4];
			_mm_storeu_si128((__m128i*) & row_uvw_[0], imm_row_uvw_);

			const __m128 xmm_row_uvw_ = _mm_mul_ps(_mm_cvtepi32_ps(imm_row_uvw_), xmm_row_uvw_scale);
			const __m128i imm_row_uvw_floor = _mm_cvtps_epi32(_mm_floor_ps(xmm_row_uvw_));
			const __m128i imm_row_uvw_ceil = _mm_cvtps_epi32(_mm_ceil_ps(xmm_row_uvw_));

			int32_t row_uvw_floor[4];
			_mm_storeu_si128((__m128i*) & row_uvw_floor[0], imm_row_uvw_floor);

			int32_t row_uvw_ceil[4];
			_mm_storeu_si128((__m128i*) & row_uvw_ceil[0], imm_row_uvw_ceil);

			if (dy01 > 0) {
				ixmax = swr_mini(ixmax, row_uvw_floor[0]);
			} else if (row_uvw_[0] != 0) {
				ixmin = swr_maxi(ixmin, row_uvw_ceil[0]);
			}

			if (dy20 > 0) {
				ixmax = swr_mini(ixmax, row_uvw_floor[1]);
			} else if (row_uvw_[1] != 0) {
				ixmin = swr_maxi(ixmin, row_uvw_ceil[1]);
			}

			if (dy01_dy20 < 0 && row_uvw_[2] >= 0) {
				ixmax = swr_mini(ixmax, -row_uvw_ceil[2]);
			} else if (dy01_dy20 > 0 && row_uvw_[2] < 0) {
				ixmin = swr_maxi(ixmin, -row_uvw_floor[2]);
			}
		}

		if (ixmin <= ixmax) {
			// Calculate normalized barycentric coordinates at ixmin of the current row of pixels.
			const __m128i imm_p0uvw_ = _mm_add_epi32(imm_row_uvw_, _mm_mullo_epi32(_mm_set1_epi32(ixmin), imm_x_duvw_));
			const __m128 xmm_p0uvw_ = _mm_mul_ps(_mm_cvtepi32_ps(imm_p0uvw_), xmm_inv_area);
			const __m128 xmm_p0uvuv = _mm_shuffle_ps(xmm_p0uvw_, xmm_p0uvw_, _MM_SHUFFLE(1, 0, 1, 0));

			// Calculate barycentric coordinates for the 4 pixels.
			const __m128 xmm_p0uv_p1uv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv0_duv1); // Barycentric coordinates of 1st and 2nd pixels
			const __m128 xmm_p2uv_p3uv = _mm_add_ps(xmm_p0uvuv, xmm_x_duv2_duv3); // Barycentric coordinates of 3rd and 4th pixels

			// Extract barycentric coordinates for each pixel
			__m128 xmm_u0123 = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p2uv_p3uv, _MM_SHUFFLE(2, 0, 2, 0));
			__m128 xmm_v0123 = _mm_shuffle_ps(xmm_p0uv_p1uv, xmm_p2uv_p3uv, _MM_SHUFFLE(3, 1, 3, 1));

			uint32_t* frameBuffer = &framebufferRow[ixmin];
			const uint32_t numPixels = (uint32_t)((ixmax - ixmin) + 1);
			const uint32_t numIter = numPixels >> 2; // 4 pixels per iteration
			for (uint32_t iIter = 0; iIter < numIter; ++iIter) {
				// Calculate the color of each pixel
				const __m128 xmm_r_p0123 = _mm_add_ps(xmm_r0, _mm_add_ps(_mm_mul_ps(xmm_dr20, xmm_u0123), _mm_mul_ps(xmm_dr10, xmm_v0123)));
				const __m128 xmm_g_p0123 = _mm_add_ps(xmm_g0, _mm_add_ps(_mm_mul_ps(xmm_dg20, xmm_u0123), _mm_mul_ps(xmm_dg10, xmm_v0123)));
				const __m128 xmm_b_p0123 = _mm_add_ps(xmm_b0, _mm_add_ps(_mm_mul_ps(xmm_db20, xmm_u0123), _mm_mul_ps(xmm_db10, xmm_v0123)));
				const __m128 xmm_a_p0123 = _mm_add_ps(xmm_a0, _mm_add_ps(_mm_mul_ps(xmm_da20, xmm_u0123), _mm_mul_ps(xmm_da10, xmm_v0123)));

				// Pack into uint8_t
				// (uint8_t){ r0, r1, r2, r3, g0, g1, g2, g3, b0, b1, b2, b3, a0, a1, a2, a3 }
				const __m128i imm_r0123_g0123_b0123_a0123_u8 = _mm_packus_epi16(
					_mm_packs_epi32(_mm_cvtps_epi32(xmm_r_p0123), _mm_cvtps_epi32(xmm_g_p0123)),
					_mm_packs_epi32(_mm_cvtps_epi32(xmm_b_p0123), _mm_cvtps_epi32(xmm_a_p0123))
				);

				// Shuffle into RGBA uint32_t
				const __m128i mask = _mm_set_epi8(15, 11, 7, 3, 14, 10, 6, 2, 13, 9, 5, 1, 12, 8, 4, 0);
				const __m128i imm_rgba_p0123_u8 = _mm_shuffle_epi8(imm_r0123_g0123_b0123_a0123_u8, mask);

				// Store
				_mm_storeu_si128((__m128i*)frameBuffer, imm_rgba_p0123_u8);

				// Move on to the next set of pixels
				xmm_u0123 = _mm_add_ps(xmm_u0123, xmm_x_du4);
				xmm_v0123 = _mm_add_ps(xmm_v0123, xmm_x_dv4);
				frameBuffer += 4;
			}

			// Handle the remainder of pixels for this row
			const uint32_t rem = numPixels & 3;
			switch (rem) {
			case 3: {
				const __m128 xmm_r_p0 = _mm_add_ss(xmm_r0, _mm_add_ss(_mm_mul_ss(xmm_dr20, xmm_u0123), _mm_mul_ss(xmm_dr10, xmm_v0123)));
				const __m128 xmm_g_p0 = _mm_add_ss(xmm_g0, _mm_add_ss(_mm_mul_ss(xmm_dg20, xmm_u0123), _mm_mul_ss(xmm_dg10, xmm_v0123)));
				const __m128 xmm_b_p0 = _mm_add_ss(xmm_b0, _mm_add_ss(_mm_mul_ss(xmm_db20, xmm_u0123), _mm_mul_ss(xmm_db10, xmm_v0123)));
				const __m128 xmm_a_p0 = _mm_add_ss(xmm_a0, _mm_add_ss(_mm_mul_ss(xmm_da20, xmm_u0123), _mm_mul_ss(xmm_da10, xmm_v0123)));
				const __m128 xmm_rrgg_p0 = _mm_shuffle_ps(xmm_r_p0, xmm_g_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_bbaa_p0 = _mm_shuffle_ps(xmm_b_p0, xmm_a_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_rgba_p0 = _mm_shuffle_ps(xmm_rrgg_p0, xmm_bbaa_p0, _MM_SHUFFLE(2, 0, 2, 0));
				const __m128i imm_rgba_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_rgba_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_rgba_p0);
				xmm_u0123 = _mm_add_ss(xmm_u0123, xmm_x_du);
				xmm_v0123 = _mm_add_ss(xmm_v0123, xmm_x_dv);
				frameBuffer++;
			} // fallthrough
			case 2: {
				const __m128 xmm_r_p0 = _mm_add_ss(xmm_r0, _mm_add_ss(_mm_mul_ss(xmm_dr20, xmm_u0123), _mm_mul_ss(xmm_dr10, xmm_v0123)));
				const __m128 xmm_g_p0 = _mm_add_ss(xmm_g0, _mm_add_ss(_mm_mul_ss(xmm_dg20, xmm_u0123), _mm_mul_ss(xmm_dg10, xmm_v0123)));
				const __m128 xmm_b_p0 = _mm_add_ss(xmm_b0, _mm_add_ss(_mm_mul_ss(xmm_db20, xmm_u0123), _mm_mul_ss(xmm_db10, xmm_v0123)));
				const __m128 xmm_a_p0 = _mm_add_ss(xmm_a0, _mm_add_ss(_mm_mul_ss(xmm_da20, xmm_u0123), _mm_mul_ss(xmm_da10, xmm_v0123)));
				const __m128 xmm_rrgg_p0 = _mm_shuffle_ps(xmm_r_p0, xmm_g_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_bbaa_p0 = _mm_shuffle_ps(xmm_b_p0, xmm_a_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_rgba_p0 = _mm_shuffle_ps(xmm_rrgg_p0, xmm_bbaa_p0, _MM_SHUFFLE(2, 0, 2, 0));
				const __m128i imm_rgba_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_rgba_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_rgba_p0);
				xmm_u0123 = _mm_add_ss(xmm_u0123, xmm_x_du);
				xmm_v0123 = _mm_add_ss(xmm_v0123, xmm_x_dv);
				frameBuffer++;
			} // fallthrough
			case 1: {
				const __m128 xmm_r_p0 = _mm_add_ss(xmm_r0, _mm_add_ss(_mm_mul_ss(xmm_dr20, xmm_u0123), _mm_mul_ss(xmm_dr10, xmm_v0123)));
				const __m128 xmm_g_p0 = _mm_add_ss(xmm_g0, _mm_add_ss(_mm_mul_ss(xmm_dg20, xmm_u0123), _mm_mul_ss(xmm_dg10, xmm_v0123)));
				const __m128 xmm_b_p0 = _mm_add_ss(xmm_b0, _mm_add_ss(_mm_mul_ss(xmm_db20, xmm_u0123), _mm_mul_ss(xmm_db10, xmm_v0123)));
				const __m128 xmm_a_p0 = _mm_add_ss(xmm_a0, _mm_add_ss(_mm_mul_ss(xmm_da20, xmm_u0123), _mm_mul_ss(xmm_da10, xmm_v0123)));
				const __m128 xmm_rrgg_p0 = _mm_shuffle_ps(xmm_r_p0, xmm_g_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_bbaa_p0 = _mm_shuffle_ps(xmm_b_p0, xmm_a_p0, _MM_SHUFFLE(0, 0, 0, 0));
				const __m128 xmm_rgba_p0 = _mm_shuffle_ps(xmm_rrgg_p0, xmm_bbaa_p0, _MM_SHUFFLE(2, 0, 2, 0));
				const __m128i imm_rgba_p0 = _mm_packus_epi16(_mm_packs_epi32(_mm_cvtps_epi32(xmm_rgba_p0), imm_zero), imm_zero);
				_mm_storeu_si32(frameBuffer, imm_rgba_p0);
				xmm_u0123 = _mm_add_ss(xmm_u0123, xmm_x_du);
				xmm_v0123 = _mm_add_ss(xmm_v0123, xmm_x_dv);
			} // fallthrough
			case 0:
			default:
				break;
			}
		}

		// Move on to the next row of pixels.
		imm_row_uvw_ = _mm_add_epi32(imm_row_uvw_, imm_y_duvw_);
		framebufferRow += ctx->m_Width;
	}
}
```

`swrDrawTriangle_SSE41()` stats:
- Average: 0.82ms
- Min: 0.78ms
- Max: 0.94ms

### Final remarks

The last thing still bothering me about the code above is the asymetry of the conditionals which calculate `xmin` and `xmax`. The reason for this asymmetry is that the equations I used for the first two (independent) barycentric coordinates were in the form:

```
bc = (A - xB) / C
```

and the equation for the third (dependent) barycentric coordinate was in the form:

```
bc = (A + xB) / C
```

Unfortunately, I have to redo the truth table using the same equation for all 3 of them, but I'm too lazy to do it at this moment. I don't expect any difference in performance, but it would be nice if everything were symmetric.

Anyway, here is the table with all the results. Thanks for reading!

Function | Average | Min | Max | Description
--|:-:|:-:|:-: |---
`swrDrawTriangle_1()` | 8.50ms | 8.43ms | 9.49ms | Original C implementation
`swrDrawTriangle_2()` | 4.76ms | 4.71ms | 4.98ms | SIMD color calculations
`swrDrawTriangle_3()` | 2.71ms | 2.65ms | 3.01ms | Moved constant calcs out of inner loops
`swrDrawTriangle_4()` | 1.59ms | 1.53ms | 1.81ms | Calculate `xmin` and `xmax` for which all coordinates will be greater than or equal to 0.
`swrDrawTriangle_5()` | 1.34ms | 1.28ms | 1.46ms | SIMD counters
`swrDrawTriangle_6()` | 1.17ms | 1.10ms | 1.23ms | 4 pixels per inner loop iteration
`swrDrawTriangle_7()` | 1.11ms | 1.06ms | 1.21ms | Rearranging the inner loop
`swrDrawTriangle_8()` | 1.03ms | 0.98ms | 1.13ms | Outer loop: `_mm_floor_ps_SSE2()` and `_mm_ceil_ps_SSE2()`
`swrDrawTriangle_9()` | 0.90ms | 0.86ms | 1.00ms | Faster RGBA32F -> RGBA8 packing
`swrDrawTriangle()` | 0.94ms | 0.90ms | 1.01ms | Treat colors as 4 independent varying attributes
`swrDrawTriangle_SSSE3()` | 0.85ms | 0.80ms | 0.91ms | SSSE3: `_mm_shuffle_epi8()`
`swrDrawTriangle_SSE41()` | 0.82ms | 0.78ms | 0.94ms | SSE4.1: `_mm_floor_ps()` and `_mm_ceil_ps()`

Bonus: If you made it this far, [here is a gist](https://gist.github.com/jdryg/725e408b948ed1d0c02f675e825c7e97#file-software_rasterizer-c) with all the code. It uses [MiniFB](https://github.com/emoon/minifb) for drawing the framebuffer into a window.
