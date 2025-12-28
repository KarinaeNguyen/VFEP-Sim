// Math3D.h
// Minimal shared math types/utilities for the VFEP visualizer.

#pragma once

struct Vec3f {
    float x, y, z;
};

static inline Vec3f v3(float x, float y, float z) { return {x, y, z}; }
