#include "Aerodynamics.h"

namespace vfep {

namespace {
constexpr double kTiny = 1e-12;

static inline bool isFinite(double x) {
    return std::isfinite(x);
}

static inline double clamp01(double x) {
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

static inline Vec3d add(Vec3d a, Vec3d b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
static inline Vec3d sub(Vec3d a, Vec3d b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static inline Vec3d mul(Vec3d a, double s) { return {a.x * s, a.y * s, a.z * s}; }

static inline double dot(Vec3d a, Vec3d b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline double len(Vec3d a) { return std::sqrt(std::max(0.0, dot(a,a))); }

static inline Vec3d normalize(Vec3d v) {
    const double L = len(v);
    if (!(L > kTiny) || !isFinite(L)) return {0.0, 0.0, 0.0};
    return {v.x / L, v.y / L, v.z / L};
}

static inline Vec3d finiteOrZero(Vec3d v) {
    if (!isFinite(v.x)) v.x = 0.0;
    if (!isFinite(v.y)) v.y = 0.0;
    if (!isFinite(v.z)) v.z = 0.0;
    return v;
}
} // namespace

AeroOutputs computeAerodynamics(const AeroInputs& in, const AeroConstants& c) {
    AeroOutputs out{};

    // Sanitize inputs (visualizer expects NaN-safe telemetry)
    const double mdot = (isFinite(in.mdot_kgps) && in.mdot_kgps > 0.0) ? in.mdot_kgps : 0.0;
    const double rpm  = (isFinite(in.vfep_rpm)  && in.vfep_rpm  > 0.0) ? in.vfep_rpm  : 0.0;

    const Vec3d draft = finiteOrZero(in.draft_vel_mps);
    const Vec3d nozzle = normalize(finiteOrZero(in.nozzle_dir_unit));

    // Jet momentum proxy (N): mdot * v_exit (v_exit in m/s)
    const double v_exit = std::max(0.0, c.v_exit_min_mps + c.v_exit_per_rpm * rpm);
    out.jet_momentum_N = mdot * v_exit;

    // Draft "drag" proxy (N): monotonic in |draft|^2
    const double v_d = len(draft);
    out.draft_drag_N = std::max(0.0, c.drag_N_per_mps2 * v_d * v_d);

    // Effective spray direction: nozzle axis plus draft-induced deflection term.
    // Deflection magnitude scales with draft speed and decreases as jet speed increases.
    const double denom = std::max(kTiny, v_exit);
    Vec3d deflect = mul(draft, -c.deflect_gain / denom);
    // If nozzle axis is degenerate, fall back to opposing draft direction (still deterministic).
    Vec3d base = (len(nozzle) > kTiny) ? nozzle : normalize(mul(draft, -1.0));
    out.spray_dir_unit = normalize(add(base, deflect));

    // Efficiency: momentum/(momentum + drag + bias) with a gentle mdot gate to keep 0 flow -> 0 eff.
    const double denomEff = out.jet_momentum_N + out.draft_drag_N + std::max(kTiny, c.eff_bias_N);
    double eff = (denomEff > kTiny) ? (out.jet_momentum_N / denomEff) : 0.0;
    // Gate by mdot so eff->0 as mdot->0 (without being a proxy in the visualizer).
    const double mdot_gate = clamp01(mdot / 0.05); // 0.05 kg/s is a reasonable "on" threshold
    out.hit_efficiency_0_1 = clamp01(eff * mdot_gate);

    return out;
}

} // namespace vfep
