#pragma once

#include <cmath>
#include <algorithm>

namespace vfep {

// Minimal 3D vector (double precision) for deterministic, dependency-free aero telemetry.
struct Vec3d {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct AeroInputs {
    double mdot_kgps = 0.0;     // agent mass flow
    double vfep_rpm  = 0.0;     // actual VFEP RPM (actuator truth)
    Vec3d draft_vel_mps;        // cross-draft velocity vector in world frame
    Vec3d nozzle_dir_unit;      // nominal nozzle axis (unit vector)
};

// Constants intentionally small and stable; tune later without touching call sites.
struct AeroConstants {
    // Effective jet speed model: v_exit = v_exit_min + rpm * v_exit_per_rpm
    double v_exit_min_mps   = 5.0;
    double v_exit_per_rpm   = 0.005; // 3000 rpm -> +15 m/s

    // Cross-draft drag proxy: drag_N = drag_N_per_mps2 * |draft|^2
    double drag_N_per_mps2  = 2.0;

    // Direction deflection gain: higher => more deflection for same draft
    double deflect_gain     = 1.0;

    // Efficiency mapping bias (N): prevents division by zero and softens low-flow behavior
    double eff_bias_N       = 1.0;
};

struct AeroOutputs {
    double hit_efficiency_0_1 = 0.0; // deterministic, monotonic in momentum vs drag
    Vec3d  spray_dir_unit;          // effective spray direction (unit vector)

    // Optional debug telemetry for HUD
    double jet_momentum_N = 0.0;
    double draft_drag_N   = 0.0;
};

// Deterministic, monotonic "momentum-to-drag" mapping consistent with Phase 2A narrative.
// No randomness; dependency-free; stable under zero/degenerate inputs.
AeroOutputs computeAerodynamics(const AeroInputs& in, const AeroConstants& c = AeroConstants());

} // namespace vfep
