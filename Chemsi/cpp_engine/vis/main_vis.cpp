
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <deque>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>

#include "Simulation.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "Math3D.h"
#include "OrbitCamera.h"

// Platform GL headers: on Windows, <GL/gl.h> requires Windows types/macros (APIENTRY/WINGDIAPI).
// Include <windows.h> first to avoid syntax errors in the Windows SDK gl.h.
#ifdef _WIN32
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  include <windows.h>
#endif

#include <GLFW/glfw3.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

static void glfw_error_callback(int error, const char* description) {
    std::fprintf(stderr, "GLFW Error %d: %s\n", error, description ? description : "(null)");
}

static int fail(const char* msg) {
    std::fprintf(stderr, "FATAL: %s\n", msg ? msg : "(null)");
    std::fprintf(stderr, "\n");
    return EXIT_FAILURE;
}

// ============================================================
// Minimal Phase-1 3D Twin (fixed-pipeline, deterministic, no assets)
// Adds Phase 1D: hit quality indicator (proxy from mdot + draft deflection + centering)
// ============================================================

// Vec3f + v3() are shared with OrbitCamera via Math3D.h
static Vec3f add(Vec3f a, Vec3f b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
static Vec3f sub(Vec3f a, Vec3f b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
static Vec3f mul(Vec3f a, float s)  { return {a.x*s, a.y*s, a.z*s}; }

static float clampf(float x, float lo, float hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}

static float dot(Vec3f a, Vec3f b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static Vec3f cross(Vec3f a, Vec3f b) { return { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x }; }
static float len(Vec3f a) { return std::sqrt(dot(a,a)); }
static Vec3f norm(Vec3f a) {
    float l = len(a);
    return (l > 1e-6f) ? mul(a, 1.0f/l) : v3(0,0,0);
}

static void set_perspective(float fovy_deg, float aspect, float znear, float zfar) {
    // OpenGL fixed pipeline expects column-major matrix.
    const float fovy_rad = fovy_deg * 3.1415926535f / 180.0f;
    const float f = 1.0f / std::tan(0.5f * fovy_rad);

    float m[16] = {};
    m[0]  = f / aspect;
    m[5]  = f;
    m[10] = (zfar + znear) / (znear - zfar);
    m[11] = -1.0f;
    m[14] = (2.0f * zfar * znear) / (znear - zfar);

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(m);
}

static void look_at(Vec3f eye, Vec3f center, Vec3f up) {
    // Minimal lookAt for fixed pipeline.
    Vec3f fwd = sub(center, eye);
    float fl = len(fwd);
    if (fl > 1e-6f) fwd = mul(fwd, 1.0f / fl);

    float ul = len(up);
    if (ul > 1e-6f) up = mul(up, 1.0f / ul);

    Vec3f s = cross(fwd, up);
    float sl = len(s);
    if (sl > 1e-6f) s = mul(s, 1.0f / sl);

    Vec3f u = cross(s, fwd);

    float m[16] = {
        s.x,  u.x,  -fwd.x, 0.0f,
        s.y,  u.y,  -fwd.y, 0.0f,
        s.z,  u.z,  -fwd.z, 0.0f,
        0.0f, 0.0f, 0.0f,   1.0f
    };

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(m);
    glTranslatef(-eye.x, -eye.y, -eye.z);
}

static void draw_wire_box(Vec3f c, Vec3f half) {
    const float x0 = c.x - half.x, x1 = c.x + half.x;
    const float y0 = c.y - half.y, y1 = c.y + half.y;
    const float z0 = c.z - half.z, z1 = c.z + half.z;

    glBegin(GL_LINES);
    // bottom
    glVertex3f(x0,y0,z0); glVertex3f(x1,y0,z0);
    glVertex3f(x1,y0,z0); glVertex3f(x1,y0,z1);
    glVertex3f(x1,y0,z1); glVertex3f(x0,y0,z1);
    glVertex3f(x0,y0,z1); glVertex3f(x0,y0,z0);
    // top
    glVertex3f(x0,y1,z0); glVertex3f(x1,y1,z0);
    glVertex3f(x1,y1,z0); glVertex3f(x1,y1,z1);
    glVertex3f(x1,y1,z1); glVertex3f(x0,y1,z1);
    glVertex3f(x0,y1,z1); glVertex3f(x0,y1,z0);
    // verticals
    glVertex3f(x0,y0,z0); glVertex3f(x0,y1,z0);
    glVertex3f(x1,y0,z0); glVertex3f(x1,y1,z0);
    glVertex3f(x1,y0,z1); glVertex3f(x1,y1,z1);
    glVertex3f(x0,y0,z1); glVertex3f(x0,y1,z1);
    glEnd();
}

static void draw_solid_box(Vec3f c, Vec3f half) {
    const float x0 = c.x - half.x, x1 = c.x + half.x;
    const float y0 = c.y - half.y, y1 = c.y + half.y;
    const float z0 = c.z - half.z, z1 = c.z + half.z;

    glBegin(GL_QUADS);
    // +Z
    glVertex3f(x0,y0,z1); glVertex3f(x1,y0,z1); glVertex3f(x1,y1,z1); glVertex3f(x0,y1,z1);
    // -Z
    glVertex3f(x1,y0,z0); glVertex3f(x0,y0,z0); glVertex3f(x0,y1,z0); glVertex3f(x1,y1,z0);
    // +X
    glVertex3f(x1,y0,z1); glVertex3f(x1,y0,z0); glVertex3f(x1,y1,z0); glVertex3f(x1,y1,z1);
    // -X
    glVertex3f(x0,y0,z0); glVertex3f(x0,y0,z1); glVertex3f(x0,y1,z1); glVertex3f(x0,y1,z0);
    // +Y
    glVertex3f(x0,y1,z1); glVertex3f(x1,y1,z1); glVertex3f(x1,y1,z0); glVertex3f(x0,y1,z0);
    // -Y
    glVertex3f(x0,y0,z0); glVertex3f(x1,y0,z0); glVertex3f(x1,y0,z1); glVertex3f(x0,y0,z1);
    glEnd();
}

static void draw_line(Vec3f a, Vec3f b) {
    glBegin(GL_LINES);
    glVertex3f(a.x,a.y,a.z);
    glVertex3f(b.x,b.y,b.z);
    glEnd();
}

static void temp_to_color(float tempC, float& r, float& g, float& b) {
    // Investor-safe, conservative ramp. 24C neutral -> hotter -> red.
    float t = clampf((tempC - 24.0f) / (120.0f - 24.0f), 0.0f, 1.0f);
    // gray -> yellow -> orange -> red
    r = 0.25f + 0.75f * t;
    g = 0.25f + 0.50f * (1.0f - std::abs(2.0f*t - 1.0f)); // peak mid
    b = 0.25f * (1.0f - t);
}

static const char* suppression_regime_text(int r) {
    switch (r) {
        case 0: return "None";
        case 1: return "Ineffective";
        case 2: return "Marginal";
        case 3: return "Effective";
        case 4: return "Overkill";
        default: return "Unknown";
    }
}

static float fire_scale_from_HRR_W(double hrrW) {
    // cube-root scaling for stability (use kW ref)
    const float hrr_kW = (float)(hrrW * 0.001);
    const float ref_kW = 1000.0f; // 1 MW reference
    float s = std::pow(std::max(hrr_kW, 0.0f) / ref_kW, 1.0f / 3.0f);
    return clampf(s, 0.10f, 2.00f);
}

// Draw a simple cone aligned to dir_unit in world space (apex at nozzle).
static void draw_cone_world(Vec3f apex, Vec3f dir_unit, float length_m, float radius_m, int slices = 16) {
    Vec3f d = norm(dir_unit);
    if (len(d) < 1e-6f || length_m <= 1e-4f || radius_m <= 1e-4f) return;

    // Orthonormal basis with z = d
    Vec3f up = (std::abs(d.y) < 0.9f) ? v3(0,1,0) : v3(1,0,0);
    Vec3f x = norm(cross(up, d));
    Vec3f y = cross(d, x);

    Vec3f base_center = add(apex, mul(d, length_m));

    glBegin(GL_TRIANGLES);
    for (int i = 0; i < slices; ++i) {
        const float a0 = (2.0f * 3.1415926535f * (float)i) / (float)slices;
        const float a1 = (2.0f * 3.1415926535f * (float)(i+1)) / (float)slices;

        Vec3f p0 = add(base_center, add(mul(x, radius_m * std::cos(a0)), mul(y, radius_m * std::sin(a0))));
        Vec3f p1 = add(base_center, add(mul(x, radius_m * std::cos(a1)), mul(y, radius_m * std::sin(a1))));

        glVertex3f(apex.x, apex.y, apex.z);
        glVertex3f(p0.x,   p0.y,   p0.z);
        glVertex3f(p1.x,   p1.y,   p1.z);
    }
    glEnd();
}

// Very simple arrow: shaft + line head.
static void draw_arrow(Vec3f origin, Vec3f dir_unit, float length_m) {
    Vec3f d = norm(dir_unit);
    if (len(d) < 1e-6f || length_m <= 1e-4f) return;

    Vec3f tip = add(origin, mul(d, length_m));

    Vec3f up = (std::abs(d.y) < 0.9f) ? v3(0,1,0) : v3(1,0,0);
    Vec3f x = norm(cross(up, d));
    Vec3f y = cross(d, x);

    const float head_len = length_m * 0.18f;
    const float head_w   = length_m * 0.06f;

    Vec3f h0 = add(tip, add(mul(d, -head_len), mul(x,  head_w)));
    Vec3f h1 = add(tip, add(mul(d, -head_len), mul(x, -head_w)));
    Vec3f h2 = add(tip, add(mul(d, -head_len), mul(y,  head_w)));
    Vec3f h3 = add(tip, add(mul(d, -head_len), mul(y, -head_w)));

    draw_line(origin, tip);
    draw_line(tip, h0);
    draw_line(tip, h1);
    draw_line(tip, h2);
    draw_line(tip, h3);
}

// Ray vs axis-aligned box intersection (slab method).
// Returns true if intersects; t_hit is distance along ray to first hit (>=0).
static bool ray_aabb_intersect(Vec3f ro, Vec3f rd_unit, Vec3f box_center, Vec3f box_half, float& t_hit) {
    Vec3f rd = rd_unit;
    // Avoid divide-by-zero; treat near-zero components as huge inv.
    auto inv = [&](float v) -> float { return (std::abs(v) > 1e-8f) ? (1.0f / v) : 1e30f; };

    float tmin = -1e30f;
    float tmax =  1e30f;

    const float bmin_x = box_center.x - box_half.x;
    const float bmax_x = box_center.x + box_half.x;
    const float bmin_y = box_center.y - box_half.y;
    const float bmax_y = box_center.y + box_half.y;
    const float bmin_z = box_center.z - box_half.z;
    const float bmax_z = box_center.z + box_half.z;

    float tx1 = (bmin_x - ro.x) * inv(rd.x);
    float tx2 = (bmax_x - ro.x) * inv(rd.x);
    tmin = std::max(tmin, std::min(tx1, tx2));
    tmax = std::min(tmax, std::max(tx1, tx2));

    float ty1 = (bmin_y - ro.y) * inv(rd.y);
    float ty2 = (bmax_y - ro.y) * inv(rd.y);
    tmin = std::max(tmin, std::min(ty1, ty2));
    tmax = std::min(tmax, std::max(ty1, ty2));

    float tz1 = (bmin_z - ro.z) * inv(rd.z);
    float tz2 = (bmax_z - ro.z) * inv(rd.z);
    tmin = std::max(tmin, std::min(tz1, tz2));
    tmax = std::min(tmax, std::max(tz1, tz2));

    if (tmax < 0.0f) return false;       // box behind ray
    if (tmin > tmax) return false;

    t_hit = (tmin >= 0.0f) ? tmin : tmax; // if inside box, take exiting hit
    return t_hit >= 0.0f;
}

int main() {
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit()) return fail("glfwInit failed");

    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

    GLFWwindow* window = glfwCreateWindow(1280, 720, "VFEP Visualizer", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return fail("glfwCreateWindow failed");
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // vsync

    // Validate OpenGL context exists.
    const GLubyte* gl_version = glGetString(GL_VERSION);
    if (!gl_version) {
        glfwDestroyWindow(window);
        glfwTerminate();
        return fail("OpenGL context validation failed (glGetString(GL_VERSION) returned null)");
    }
    std::fprintf(stderr, "OpenGL Vendor:   %s\n", glGetString(GL_VENDOR));
    std::fprintf(stderr, "OpenGL Renderer: %s\n", glGetString(GL_RENDERER));
    std::fprintf(stderr, "OpenGL Version:  %s\n", gl_version);

    bool imgui_ctx = false;
    bool imgui_glfw = false;
    bool imgui_gl3 = false;

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    imgui_ctx = true;

    ImGui::StyleColorsDark();

    if (!ImGui_ImplGlfw_InitForOpenGL(window, true)) {
        if (imgui_ctx) ImGui::DestroyContext();
        glfwDestroyWindow(window);
        glfwTerminate();
        return fail("ImGui_ImplGlfw_InitForOpenGL failed");
    }
    imgui_glfw = true;

    if (!ImGui_ImplOpenGL3_Init(glsl_version)) {
        if (imgui_glfw) ImGui_ImplGlfw_Shutdown();
        if (imgui_ctx) ImGui::DestroyContext();
        glfwDestroyWindow(window);
        glfwTerminate();
        return fail("ImGui_ImplOpenGL3_Init failed");
    }
    imgui_gl3 = true;

    vfep::Simulation sim;
    bool running = false;

    // --- Phase-1 3D Twin camera (deterministic, ImGui-controlled) ---
    OrbitCamera cam;
    cam.yaw_deg   = 35.0f;
    cam.pitch_deg = 20.0f;
    cam.dist      = 8.0f;
    cam.target    = v3(0.0f, 1.2f, 0.0f);
// --- Rack color (smooth HRR-driven red ramp) ---
float rack_red_t = 0.0f; // 0=dark rack, 1=full red

    // OrbitCamera stores all orbit tuning + mouse state internally.
    // --- Phase-1 scene layout (meters, simple boxes) ---
    Vec3f warehouse_half = v3(6.0f, 3.0f, 6.0f);          // 12m x 6m x 12m volume
    Vec3f rack_center    = v3(0.0f, 1.0f, 0.0f);          // centered
    Vec3f rack_half      = v3(0.6f, 1.0f, 0.4f);          // ~1.2m x 2.0m x 0.8m
    Vec3f fire_center    = v3(0.0f, 0.6f, 0.7f);          // in front of rack

    // --- Spray/nozzle parameters (Phase 1A + 1B) ---
    Vec3f nozzle_pos     = v3(-2.0f, 1.5f, -2.0f);
    Vec3f nozzle_dir     = v3(0.7f, -0.15f, 0.7f); // will be normalized
    float mdot_ref       = 0.15f;  // kg/s at "full scale" visualization
    float spray_L0       = 0.6f;   // base cone length
    float spray_L1       = 3.2f;   // added length at eff=1
    float spray_R0       = 0.10f;  // base cone radius
    float spray_R1       = 0.28f;  // added radius at eff=1
    float spray_max_len  = 8.0f;

    // --- Hit marker parameters (Phase 1B) ---
    float hit_marker_base = 0.06f;  // meters (half extent scale base)
    float hit_marker_gain = 0.20f;  // meters added at eff=1

    // --- Cross-draft (Phase 1A) ---
    Vec3f draft_vel_mps  = v3(0.0f, 0.0f, 0.0f);
    float draft_arrow_scale = 0.7f; // meters per (m/s), clamped later

    // Phase 1C/1D: simple deterministic draft deflection
    float draft_deflect_gain = 0.35f; // dimensionless

    // Fixed simulation timestep
    double dt = 0.05;       // seconds
    double simTime = 0.0;   // simulated time (seconds)

    // Wall-time stepping (decouple sim stepping from render FPS)
    double wall_prev = glfwGetTime();
    double accum_s = 0.0;

    // Operator diagnostics
    int last_substeps = 0;
    bool dropped_accum = false;


// ---------------- VFEP-Updated-Console-0.2 UI state ----------------
struct UIEvent { double t_s; std::string text; };

struct UIState {
    // Core controls
    bool running = false;
    double dt_s = 0.05;

    bool step_requested = false;
    bool ignite_requested = false;
    bool suppress_requested = false;
    bool reset_accumulator = false;

    // Fonts (optional; nullptr => default font)
    ImFont* font_base = nullptr;
    ImFont* font_big  = nullptr;
    ImFont* font_mono = nullptr;

    // Trend history (last sample)
    bool has_prev = false;
    double prev_HRR_kW = 0.0;
    double prev_T_C = 0.0;
    double prev_agent_mdot = 0.0;

    // Event log (ring buffer)
    std::deque<UIEvent> events;
    size_t events_cap = 200;

    // Transition tracking for event generation
    bool prev_concluded = false;
    bool prev_dropped_accum = false;
};

UIState ui;
ui.running = running;    // bind to existing sim loop var
ui.dt_s    = dt;

auto push_event = [&](double t_s, const char* fmt, ...) {
    char buf[512];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);

    ui.events.push_back(UIEvent{ t_s, std::string(buf) });
    while (ui.events.size() > ui.events_cap) ui.events.pop_front();
};

enum class Trend { Down, Flat, Up };

auto trend_of = [&](double curr, double prev, double rel_eps, double abs_eps) -> Trend {
    const double eps = std::max(std::abs(curr) * rel_eps, abs_eps);
    if (curr > prev + eps) return Trend::Up;
    if (curr < prev - eps) return Trend::Down;
    return Trend::Flat;
};

auto trend_glyph = [&](Trend t) -> const char* {
    switch (t) {
        case Trend::Up:   return "↑";
        case Trend::Down: return "↓";
        default:          return "→";
    }
};

// Font + style (demo-ready). If custom font files are unavailable, ImGui will fall back safely.
{
    ImGuiIO& io = ImGui::GetIO();

    // ---- Demo readability defaults (font-based scaling) ----
    // Tuned so the console is readable from a distance on a typical 1080p/1440p projector.
    constexpr float kFontBasePx = 30.0f;   // main UI text
    constexpr float kFontBigPx  = 46.0f;   // headers / key metrics
    constexpr float kFontMonoPx = 28.0f;   // console / numbers

    // Try local assets first (recommended). If missing, fall back to Segoe UI (Windows).
    ui.font_base = io.Fonts->AddFontFromFileTTF("assets/fonts/Inter-Regular.ttf",    kFontBasePx);
    ui.font_big  = io.Fonts->AddFontFromFileTTF("assets/fonts/Inter-SemiBold.ttf",  kFontBigPx);
    ui.font_mono = io.Fonts->AddFontFromFileTTF("assets/fonts/JetBrainsMono-Regular.ttf", kFontMonoPx);

    // Windows fallback (Segoe UI)
    if (!ui.font_base) ui.font_base = io.Fonts->AddFontFromFileTTF("C:\Windows\Fonts\segoeui.ttf",  kFontBasePx);
    if (!ui.font_big)  ui.font_big  = io.Fonts->AddFontFromFileTTF("C:\Windows\Fonts\segoeuib.ttf", kFontBigPx);
    if (!ui.font_mono) ui.font_mono = io.Fonts->AddFontFromFileTTF("C:\Windows\Fonts\consola.ttf",  kFontMonoPx);

    // If everything failed, ImGui will still use its default font; keep running regardless.

    ImGuiStyle& style = ImGui::GetStyle();

    // Give the layout more breathing room at larger fonts.
    style.WindowPadding     = ImVec2(18, 18);
    style.FramePadding      = ImVec2(14, 10);
    style.ItemSpacing       = ImVec2(14, 12);
    style.ItemInnerSpacing  = ImVec2(10, 8);
    style.IndentSpacing     = 18.0f;
    style.ScrollbarSize     = 22.0f;
    style.GrabMinSize       = 18.0f;

    // Important: avoid global scaling (we want font-based scaling to preserve crispness).
    // io.FontGlobalScale = 1.0f;
}

// -------------------------------------------------------------------
    vfep::Observation last_obs = sim.observe();
    // VFEP-Updated-Console-0.2: no plot history (trends computed from last sample only)

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        // --- advance sim (wall-time accumulator) ---
        const double wall_now = glfwGetTime();
        double wall_dt = wall_now - wall_prev;
        wall_prev = wall_now;

        wall_dt = std::clamp(wall_dt, 0.0, 0.1);
        dt = std::clamp(dt, 0.001, 1.0);

        bool advanced_this_frame = false;

        if (running && !sim.isConcluded()) {
            accum_s += wall_dt;

            constexpr int kMaxSubstepsPerFrame = 20;
            int substeps = 0;
            dropped_accum = false;

            while (accum_s >= dt && substeps < kMaxSubstepsPerFrame && !sim.isConcluded()) {
                sim.step(dt);
                simTime += dt;
                last_obs = sim.observe();
                // VFEP-Updated-Console-0.2: no plot history (trends computed from last sample only)
                accum_s -= dt;
                ++substeps;
                advanced_this_frame = true;
            }

            last_substeps = substeps;

            if (substeps == kMaxSubstepsPerFrame) {
                accum_s = 0.0;
                dropped_accum = true;
            }
        } else {
            last_substeps = 0;
            dropped_accum = false;
        }

        if (!advanced_this_frame) {
            last_obs = sim.observe();
        }

        // Update truth telemetry vectors used by the visualizer (Phase 2A)
        draft_vel_mps = v3((float)last_obs.draft_vel_mps_x,
                           (float)last_obs.draft_vel_mps_y,
                           (float)last_obs.draft_vel_mps_z);

        // --- ImGui frame ---
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();


    // --- Orbit camera controls (encapsulated), gated by ImGui capture ---
    {
        ImGuiIO& io = ImGui::GetIO();
        cam.updateFromInput(window, io, wall_dt, rack_center, /*lock_target=*/true);
    }

}


        // Controls

// ---------------- VFEP-Updated-Console-0.2: Single console dashboard ----------------
auto hrr_kW = [&](const vfep::Observation& o) -> double { return 1e-3 * (double)o.effective_HRR_W; };
auto temp_C = [&](const vfep::Observation& o) -> double { return (double)o.T_K - 273.15; };
auto agent_mdot = [&](const vfep::Observation& o) -> double { return (double)o.agent_mdot_kgps; };

ImGui::SetNextWindowSize(ImVec2(560, 920), ImGuiCond_FirstUseEver);
ImGui::Begin("VFEP Console");

// Header strip (time-centric "now")
{
    if (ui.font_big) ImGui::PushFont(ui.font_big);
    ImGui::Text("t = %.2f s", simTime);
    if (ui.font_big) ImGui::PopFont();

    ImGui::SameLine();
    ImGui::Dummy(ImVec2(18, 0));
    ImGui::SameLine();

    const bool concluded = sim.isConcluded();
    const char* state_txt = concluded ? "CONCLUDED" : (running ? "RUNNING" : "PAUSED");
    ImVec4 state_col = concluded ? ImVec4(0.70f, 0.70f, 0.70f, 1.0f)
                                 : (running ? ImVec4(0.20f, 0.85f, 0.35f, 1.0f)
                                            : ImVec4(0.95f, 0.80f, 0.20f, 1.0f));
    ImGui::TextColored(state_col, "%s", state_txt);

    ImGui::SameLine();
    ImGui::Dummy(ImVec2(18, 0));
    ImGui::SameLine();
    ImGui::Text("dt=%.3f s", dt);

    ImGui::SameLine();
    ImGui::Dummy(ImVec2(12, 0));
    ImGui::SameLine();
    ImGui::Text("substeps=%d", last_substeps);

    if (dropped_accum) {
        ImGui::SameLine();
        ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "Realtime: DROPPED");
    }
}

ImGui::Separator();

if (ImGui::BeginTabBar("vfep_console_tabs"))
{
    // ---------------- Control tab ----------------
    if (ImGui::BeginTabItem("Control"))
    {
        // Controls: primary actions (no scrolling)
        ImGui::TextDisabled("Controls");
        {
            const float btn_w = 140.0f;
            if (ImGui::Button(running ? "Pause" : "Run", ImVec2(btn_w, 0))) {
                running = !running;
                ui.running = running;
                push_event(simTime, running ? "Run" : "Pause");
            }

            ImGui::SameLine();
            if (ImGui::Button("Step", ImVec2(btn_w, 0))) {
                if (!sim.isConcluded()) {
                    sim.step(dt);
                    simTime += dt;
                    last_obs = sim.observe();
                    accum_s = 0.0;
                    last_substeps = 1;
                    dropped_accum = false;
                    push_event(simTime, "Step");
                }
            }

            ImGui::SameLine();
            if (ImGui::Button("Ignite Fire", ImVec2(btn_w, 0))) {
                if (!sim.isConcluded()) {
                    sim.commandIgniteOrIncreasePyrolysis();
                    push_event(simTime, "Ignite command issued");
                }
            }

            if (ImGui::Button("Start Suppression", ImVec2(btn_w * 1.5f, 0))) {
                if (!sim.isConcluded()) {
                    sim.commandStartSuppression();
                    push_event(simTime, "Suppression started");
                }
            }
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        // dt control (core only)
        ImGui::TextDisabled("Timestep");
        {
            float dt_f = (float)dt;
            if (ImGui::SliderFloat("dt (s)", &dt_f, 0.005f, 0.200f, "%.3f")) {
                dt = (double)dt_f;
                ui.dt_s = dt;
                accum_s = 0.0;
                dropped_accum = false;
                push_event(simTime, "dt set to %.3f s", dt);
            }

            ImGui::SameLine();
            if (ImGui::Button("-")) {
                dt = std::max(0.001, dt - 0.005);
                ui.dt_s = dt;
                accum_s = 0.0;
                dropped_accum = false;
                push_event(simTime, "dt set to %.3f s", dt);
            }
            ImGui::SameLine();
            if (ImGui::Button("+")) {
                dt = std::min(1.0, dt + 0.005);
                ui.dt_s = dt;
                accum_s = 0.0;
                dropped_accum = false;
                push_event(simTime, "dt set to %.3f s", dt);
            }
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        // Live telemetry (no plots): value + trend arrows
        ImGui::TextDisabled("Live Telemetry");
        {
            const double curr_hrr = hrr_kW(last_obs);
            const double curr_Tc  = temp_C(last_obs);
            const double curr_mdot = agent_mdot(last_obs);

            Trend th = Trend::Flat, tT = Trend::Flat, tm = Trend::Flat;
            if (ui.has_prev) {
                th = trend_of(curr_hrr, ui.prev_HRR_kW, 0.01, 0.05);   // 1% or 0.05 kW
                tT = trend_of(curr_Tc,  ui.prev_T_C,   0.005, 0.02);  // 0.5% or 0.02 C
                tm = trend_of(curr_mdot,ui.prev_agent_mdot, 0.02, 1e-4);
            }

            if (ImGui::BeginTable("telemetry", 3, ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_RowBg))
            {
                ImGui::TableSetupColumn("Signal");
                ImGui::TableSetupColumn("Value");
                ImGui::TableSetupColumn("Trend");
                ImGui::TableHeadersRow();

                auto row = [&](const char* label, const char* fmt, double v, const char* trend) {
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(label);
                    ImGui::TableSetColumnIndex(1); ImGui::Text(fmt, v);
                    ImGui::TableSetColumnIndex(2); ImGui::TextUnformatted(trend);
                };

                row("HRR",      "%.1f kW", curr_hrr, ui.has_prev ? trend_glyph(th) : "—");
                row("Temp",     "%.2f C",  curr_Tc,  ui.has_prev ? trend_glyph(tT) : "—");
                row("Agent Flow","%.4f kg/s", curr_mdot, ui.has_prev ? trend_glyph(tm) : "—");

                ImGui::EndTable();
            }

            // Update trend memory AFTER rendering using the current observation
            ui.prev_HRR_kW = curr_hrr;
            ui.prev_T_C = curr_Tc;
            ui.prev_agent_mdot = curr_mdot;
            ui.has_prev = true;
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        // Status / events: fixed-height log, minimal scrolling
        ImGui::TextDisabled("Status & Events");
        {
            ImGui::BeginChild("events_child", ImVec2(0, 190), true, ImGuiWindowFlags_None);
            if (ui.events.empty()) {
                ImGui::TextDisabled("No events yet.");
            } else {
                // newest last (more natural reading)
                for (const auto& e : ui.events) {
                    ImGui::Text("[%.2fs] %s", e.t_s, e.text.c_str());
                }
                if (ImGui::GetScrollY() >= ImGui::GetScrollMaxY() - 5.0f) {
                    ImGui::SetScrollHereY(1.0f);
                }
            }
            ImGui::EndChild();
        }

        ImGui::EndTabItem();
    }

    // ---------------- Physics & Chemistry tab ----------------
    if (ImGui::BeginTabItem("Physics & Chemistry"))
    {
        ImGui::TextDisabled("Environment");
        if (ImGui::BeginTable("env_tbl", 2, ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_RowBg)) {
            ImGui::TableNextRow(); ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted("Temperature");
            ImGui::TableSetColumnIndex(1); ImGui::Text("%.2f K  (%.2f C)", (double)last_obs.T_K, temp_C(last_obs));

            ImGui::TableNextRow(); ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted("O2");
            ImGui::TableSetColumnIndex(1); ImGui::Text("%.3f vol%%", (double)last_obs.O2_volpct);

            ImGui::TableNextRow(); ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted("CO2");
            ImGui::TableSetColumnIndex(1); ImGui::Text("%.3f vol%%", (double)last_obs.CO2_volpct);

            ImGui::TableNextRow(); ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted("H2O");
            ImGui::TableSetColumnIndex(1); ImGui::Text("%.3f vol%%", (double)last_obs.H2O_volpct);
            ImGui::EndTable();
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        ImGui::TextDisabled("Fire / Suppression");
        if (ImGui::BeginTable("fire_tbl", 2, ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_RowBg)) {
            ImGui::TableNextRow(); ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted("HRR (effective)");
            ImGui::TableSetColumnIndex(1); ImGui::Text("%.1f kW", hrr_kW(last_obs));

            ImGui::TableNextRow(); ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted("Knockdown");
            ImGui::TableSetColumnIndex(1); ImGui::Text("%.3f", (double)last_obs.knockdown_0_1);

            ImGui::TableNextRow(); ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted("Effective Exposure");
            ImGui::TableSetColumnIndex(1); ImGui::Text("%.4f kg", (double)last_obs.effective_exposure_kg);

            ImGui::TableNextRow(); ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted("Regime");
            ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(suppression_regime_text(last_obs.suppression_regime));
            ImGui::EndTable();
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        ImGui::TextDisabled("Actuation / Delivery");
        if (ImGui::BeginTable("act_tbl", 2, ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_RowBg)) {
            ImGui::TableNextRow(); ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted("Agent mdot");
            ImGui::TableSetColumnIndex(1); ImGui::Text("%.6f kg/s", (double)last_obs.agent_mdot_kgps);

            ImGui::TableNextRow(); ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted("Net delivered mdot");
            ImGui::TableSetColumnIndex(1); ImGui::Text("%.6f kg/s", (double)last_obs.net_delivered_mdot_kgps);

            ImGui::TableNextRow(); ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted("Draft vel (z)");
            ImGui::TableSetColumnIndex(1); ImGui::Text("%.3f m/s", (double)last_obs.draft_vel_mps_z);
            ImGui::EndTable();
        }

        ImGui::EndTabItem();
    }

    // ---------------- Data tab ----------------
    if (ImGui::BeginTabItem("Data"))
    {
        if (ui.font_mono) ImGui::PushFont(ui.font_mono);

        ImGui::TextDisabled("Current Observation (selected fields)");

        if (ImGui::BeginTable("obs_tbl", 3, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_SizingFixedFit))
        {
            ImGui::TableSetupColumn("Field");
            ImGui::TableSetupColumn("Value");
            ImGui::TableSetupColumn("Units");
            ImGui::TableHeadersRow();

            auto row = [&](const char* k, double v, const char* units, const char* fmt) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(k);
                ImGui::TableSetColumnIndex(1); ImGui::Text(fmt, v);
                ImGui::TableSetColumnIndex(2); ImGui::TextUnformatted(units);
            };

            row("t", simTime, "s", "%.3f");
            row("T_K", (double)last_obs.T_K, "K", "%.3f");
            row("HRR_W (effective)", (double)last_obs.effective_HRR_W, "W", "%.1f");
            row("O2", (double)last_obs.O2_volpct, "vol%", "%.4f");
            row("CO2", (double)last_obs.CO2_volpct, "vol%", "%.4f");
            row("H2O", (double)last_obs.H2O_volpct, "vol%", "%.4f");
            row("agent_mdot", (double)last_obs.agent_mdot_kgps, "kg/s", "%.6f");
            row("net_delivered_mdot", (double)last_obs.net_delivered_mdot_kgps, "kg/s", "%.6f");
            row("effective_exposure", (double)last_obs.effective_exposure_kg, "kg", "%.6f");
            row("knockdown", (double)last_obs.knockdown_0_1, "0..1", "%.6f");

            ImGui::EndTable();
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        // Optional: sector arrays (moved here to avoid density in the console)
        if (ImGui::CollapsingHeader("Suppression sectors (arrays)", ImGuiTreeNodeFlags_DefaultOpen))
        {
            constexpr int N = vfep::Observation::kNumSuppressionSectors;
            if (ImGui::BeginTable("sector_tbl", 4, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_SizingFixedFit))
            {
                ImGui::TableSetupColumn("i");
                ImGui::TableSetupColumn("occ");
                ImGui::TableSetupColumn("loa");
                ImGui::TableSetupColumn("KD_target");
                ImGui::TableHeadersRow();

                for (int i = 0; i < N; ++i) {
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("%d", i);
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%.3f", (double)last_obs.sector_occlusion_0_1[i]);
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%.3f", (double)last_obs.sector_line_attack_0_1[i]);
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%.3f", (double)last_obs.sector_knockdown_target_0_1[i]);
                }
                ImGui::EndTable();
            }
        }

        if (ui.font_mono) ImGui::PopFont();
        ImGui::EndTabItem();
    }

    ImGui::EndTabBar();
}

// VFEP-Updated-Console-0.2 event generation (automatic)
{
    const bool concluded = sim.isConcluded();
    if (!ui.prev_concluded && concluded) {
        push_event(simTime, "Simulation concluded");
    }
    ui.prev_concluded = concluded;

    if (!ui.prev_dropped_accum && dropped_accum) {
        push_event(simTime, "Realtime accumulator dropped (frame cap)");
    }
    ui.prev_dropped_accum = dropped_accum;

    
static bool prev_warn_fb = false;
static bool prev_warn_gh = false;

if (last_obs.warn_fully_blocked) {
    if (!prev_warn_fb) push_event(simTime, "Warning: fully blocked sectors present");
    prev_warn_fb = true;
} else {
    prev_warn_fb = false;
}

if (last_obs.warn_glancing_hold) {
    if (!prev_warn_gh) push_event(simTime, "Warning: glancing hold sectors present");
    prev_warn_gh = true;
} else {
    prev_warn_gh = false;
}
}

ImGui::End();
// -------------------------------------------------------------------
        ImGui::Render();

        int fb_w = 0, fb_h = 0;
        glfwGetFramebufferSize(window, &fb_w, &fb_h);

        if (fb_w > 0 && fb_h > 0) {
            glViewport(0, 0, fb_w, fb_h);

            glEnable(GL_DEPTH_TEST);
            glDepthFunc(GL_LESS);
            glDisable(GL_CULL_FACE);

            glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            const float aspect = (fb_h > 0) ? (float)fb_w / (float)fb_h : 1.0f;
            set_perspective(55.0f, aspect, 0.05f, 100.0f);

            const Vec3f eye = cam.eye();
            look_at(eye, cam.target, v3(0.0f, 1.0f, 0.0f));

            // Warehouse wireframe
            glColor3f(0.35f, 0.35f, 0.35f);
            draw_wire_box(v3(0.0f, warehouse_half.y, 0.0f), warehouse_half);

            // Rack (dark by default; smoothly turns red with HRR during fire)
{
    // Base rack color (near black, but not pure black for projector contrast)
    const float base_r = 0.10f;
    const float base_g = 0.10f;
    const float base_b = 0.10f;

    // Target "fire red" color (demo-friendly, not neon)
    const float fire_r = 0.88f;
    const float fire_g = 0.08f;
    const float fire_b = 0.08f;

    // Normalize HRR into [0,1] intensity. Tune hrr_full_W to your typical demo range.
    constexpr float hrr_full_W = 300000.0f; // 300 kW => full red (adjust if needed)
    const float hrr_W = (float)last_obs.effective_HRR_W;
    const float target_t = std::clamp(hrr_W / hrr_full_W, 0.0f, 1.0f);

    // Smooth with an exponential low-pass to avoid flicker
    // tau = response time in seconds (smaller = snappier, larger = smoother)
    constexpr float tau_s = 0.45f;
    const float dt_s = (float)wall_dt; // wall_dt already clamped earlier
    const float alpha = 1.0f - std::exp(-dt_s / tau_s);

    rack_red_t = rack_red_t + (target_t - rack_red_t) * alpha;

    const float rr = base_r + (fire_r - base_r) * rack_red_t;
    const float rg = base_g + (fire_g - base_g) * rack_red_t;
    const float rb = base_b + (fire_b - base_b) * rack_red_t;

    glColor3f(rr, rg, rb);
    draw_solid_box(rack_center, rack_half);
}
glColor3f(0.05f, 0.05f, 0.05f);
            draw_wire_box(rack_center, rack_half);

            // Fire proxy (Phase 2B/2C): visualize effective HRR and spatial knockdown
            const double hrr_vis_W = (std::isfinite(last_obs.effective_HRR_W) && last_obs.effective_HRR_W > 0.0)
                ? last_obs.effective_HRR_W
                : last_obs.HRR_W;

            const float fire_s = fire_scale_from_HRR_W(hrr_vis_W);
            Vec3f fire_half = mul(v3(0.35f, 0.45f, 0.35f), fire_s);

            if (hrr_vis_W > 1.0) {
                // Render as 2x2 sector cubes (no new geometry types) to make spatial suppression visible.
                Vec3f sub_half = v3(fire_half.x * 0.48f, fire_half.y, fire_half.z * 0.48f);

                const float sx[4] = {-1.0f, +1.0f, -1.0f, +1.0f};
                const float sz[4] = {-1.0f, -1.0f, +1.0f, +1.0f};

                for (int i = 0; i < 4; ++i) {
                    const float kd = clampf((float)last_obs.sector_knockdown_0_1[i], 0.0f, 1.0f);
                    const float intensity = clampf(0.20f + 0.80f * (1.0f - kd), 0.0f, 1.0f);

                    glColor3f(0.85f * intensity, 0.25f * intensity, 0.05f * intensity);

                    Vec3f c = fire_center;
                    c.x += sx[i] * sub_half.x;
                    c.z += sz[i] * sub_half.z;
                    draw_solid_box(c, sub_half);
                }

                // Outline whole fire volume for readability.
                glColor3f(0.15f, 0.05f, 0.02f);
                draw_wire_box(fire_center, fire_half);
            }

// Cross-draft arrow (centered near rack)
            {
                const float mag = len(draft_vel_mps);
                const float L = clampf(draft_arrow_scale * mag, 0.2f, 4.0f);
                glColor3f(0.10f, 0.75f, 0.75f);
                draw_arrow(add(rack_center, v3(0.0f, rack_half.y + 0.3f, 0.0f)), draft_vel_mps, L);
            }

            // Nozzle marker
            glColor3f(0.55f, 0.55f, 0.60f);
            draw_wire_box(nozzle_pos, v3(0.06f, 0.06f, 0.06f));

            
// Spray cone (truth telemetry)
const float eff_draw = clampf((float)last_obs.hit_efficiency_0_1, 0.0f, 1.0f);
const float cone_len = clampf(spray_L0 + spray_L1 * eff_draw, 0.0f, spray_max_len);
const float cone_rad = clampf(spray_R0 + spray_R1 * eff_draw, 0.0f, 3.0f);

// Phase 2A: effective spray direction comes from sim telemetry.
const Vec3f eff_dir = norm(v3((float)last_obs.spray_dir_unit_x,
                              (float)last_obs.spray_dir_unit_y,
                              (float)last_obs.spray_dir_unit_z));

const Vec3f nozzle_dir_n =
(len(nozzle_dir) > 1e-6f) ? norm(nozzle_dir) : eff_dir;

if (last_obs.agent_mdot_kgps > 1e-6) {
    // Purple-ish tint for Purple K (no blending required)
    glColor3f(0.55f, 0.25f, 0.70f);
    draw_cone_world(nozzle_pos, eff_dir, cone_len, cone_rad, 18);

    // Centerline for readability (effective direction)
    glColor3f(0.35f, 0.18f, 0.45f);
    draw_line(nozzle_pos, add(nozzle_pos, mul(eff_dir, cone_len)));

    // Nominal nozzle direction (faint) to show draft deflection
    glColor3f(0.22f, 0.22f, 0.25f);
    draw_line(nozzle_pos, add(nozzle_pos, mul(nozzle_dir_n, cone_len)));
}

// Phase 1B/1D: Hit marker on rack where effective spray ray intersects the rack box.
// Phase 1D adds a "hit quality" proxy: eff * alignment * centering.
{
    float t_hit = 0.0f;
    if (len(eff_dir) > 1e-6f && ray_aabb_intersect(nozzle_pos, eff_dir, rack_center, rack_half, t_hit)) {
        Vec3f hit = add(nozzle_pos, mul(eff_dir, t_hit));

        const float marker_half = clampf(hit_marker_base + hit_marker_gain * eff_draw, 0.01f, 0.5f);
        Vec3f mh = v3(marker_half, marker_half, marker_half);

        // Hit quality components
        const float hit_align = clampf(dot(nozzle_dir_n, eff_dir), 0.0f, 1.0f);

        // Centering score on the impacted face
        const float nx = (rack_half.x > 1e-6f) ? (hit.x - rack_center.x) / rack_half.x : 0.0f;
        const float ny = (rack_half.y > 1e-6f) ? (hit.y - rack_center.y) / rack_half.y : 0.0f;
        const float nz = (rack_half.z > 1e-6f) ? (hit.z - rack_center.z) / rack_half.z : 0.0f;

        const float ax = std::abs(nx), ay = std::abs(ny), az = std::abs(nz);

        float u = 0.0f, v = 0.0f;
        if (ax >= ay && ax >= az) { u = ny; v = nz; }
        else if (ay >= ax && ay >= az) { u = nx; v = nz; }
        else { u = nx; v = ny; }

        const float dist = std::sqrt(u*u + v*v);
        const float hit_center = clampf(1.0f - dist / 1.41421356f, 0.0f, 1.0f);

        const float hit_quality = eff_draw; // truth-driven (already includes aero + draft effects)

        // Color encodes hit quality (low -> warm/purple, high -> greenish)
        const float r = 0.80f * (1.0f - hit_quality) + 0.25f * hit_quality;
        const float g = 0.20f * (1.0f - hit_quality) + 0.90f * hit_quality;
        const float b = 0.75f * (1.0f - hit_quality) + 0.30f * hit_quality;

        glColor3f(r, g, b);
        draw_solid_box(hit, mh);

        glColor3f(0.08f, 0.03f, 0.10f);
        draw_wire_box(hit, mh);
    }
}

            // UI overlay
            glDisable(GL_DEPTH_TEST);
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        }

        glfwSwapBuffers(window);
    }
    if (imgui_gl3) ImGui_ImplOpenGL3_Shutdown();
    if (imgui_glfw) ImGui_ImplGlfw_Shutdown();
    if (imgui_ctx) ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
