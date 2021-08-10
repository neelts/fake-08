//
//  Lol Engine
//
//  Copyright © 2010—2020 Sam Hocevar <sam@hocevar.net>
//            © 2013—2015 Benjamin “Touky” Huet <huet.benjamin@gmail.com>
//            © 2013—2015 Guillaume Bittoun <guillaume.bittoun@gmail.com>
//
//  Lol Engine is free software. It comes without any warranty, to
//  the extent permitted by applicable law. You can redistribute it
//  and/or modify it under the terms of the Do What the Fuck You Want
//  to Public License, Version 2, as published by the WTFPL Task Force.
//  See http://www.wtfpl.net/ for more details.
//

//Quick and dirty port removing vec_t and assuming only a single value passed
//jtothebell August 2021
//source file: https://github.com/lolengine/lol/blob/2c61a5e80521c91e52c2c1ad8a1e621e43ae1a1a/include/lol/private/math/perlin.h
#pragma once

//#include "functions.h" // for clamp()
#include "gradient.h"

#include <cmath> // std::sqrt

float perlin_clamp (float val, float lo, float hi) {
	val = (val > hi) ? hi : val;
    return (val < lo) ? lo : val;
}

namespace lol
{

template<int N>
class perlin_noise : public gradient_provider<N>
{
public:
    perlin_noise()
      : gradient_provider<N>()
    {
    }

    perlin_noise(int seed)
      : gradient_provider<N>(seed)
    {
    }

    /* Evaluate noise at a given point */
    inline float eval(float position) const
    {
        using std::sqrt;

        /* Compute the containing hypercube origin */
        int origin = (int)position - (position < 0);

        float delta = position - (float)origin;

        /* Apply a smooth step to delta and store it in “t”. */
        float t = delta;
        t = ((6.f * t - 15.f)
                  * t + 10.f) * (t * t * t);
        /* DEBUG: original Perlin noise polynomial */
        //t = (3.f - 2.f * t) * t * t;

        /* Premultiply and predivide (1-t)/t and t/(1-t) into “u” and “v”. */
        float u;
        float v;
        float multiplier = 1.f;
        float f = perlin_clamp(t, 0.001f, 0.999f);

        multiplier *= (1.f - f);
        u = (1.f - f) / f;
        v = f / (1.f - f);

        float ret = 0.f;

        /* Compute all gradient contributions, for each of the 2^N corners
         * of the hypercube. */
        for (int i = 0; ; ++i)
        {
            /* Accumulate Perlin noise */
            ret += multiplier * delta * this->get_gradient(origin);

            /* Avoid buffer overflow below in origin[bit] */
            if (i + 1 == (1 << 1))
                break;

            /* Don’t use the binary pattern for “i” but use its Gray code
             * “j” instead, so we know we only have one component to alter
             * in “origin” and in “delta”. We know which bit was flipped by
             * looking at “k”, the Gray code for the next value of “i”. */
            int j = i ^ (i >> 1);
            int k = (i + 1) ^ ((i + 1) >> 1);

            int bit = 0;
            while ((j ^ k) > (1 << bit))
                ++bit;

            origin += j > k ? -1 : 1;
            delta += j > k ? 1.f : -1.f;
            multiplier *= (j > k ? u : v);
        }

        return sqrt(2.f) * ret;
    }
};

}

