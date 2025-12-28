// OrbitCamera.h
// Encapsulated orbit camera used by main_vis.cpp.

#pragma once

#include <algorithm>
#include <cmath>

#include <GLFW/glfw3.h>
#include <imgui.h>

#include "Math3D.h"

// A simple orbit camera (yaw/pitch/distance around a target).
//
// Controls (when ImGui is not capturing the relevant input):
//   - RMB drag: orbit
//   - Mouse wheel: exponential zoom
//   - Arrow keys: orbit
//   - W/S: zoom
//   - Shift: fast modifier
//
// Intended usage:
//   OrbitCamera cam;
//   cam.setTarget(...);
//   cam.updateFromInput(window, ImGui::GetIO(), wall_dt, fixedTarget);
//   Vec3f eye = cam.eye();

struct OrbitCamera {
    // Pose
    float yaw_deg   = 35.0f;
    float pitch_deg = 20.0f;
    float dist      = 8.0f;
    Vec3f target    = v3(0.0f, 1.2f, 0.0f);

    // Mouse tuning
    float orbit_sens_deg_per_px = 0.18f;
    float zoom_step             = 0.90f; // multiplicative per wheel "tick"

    // Distance clamp
    float min_dist = 2.0f;
    float max_dist = 20.0f;

    // Keyboard tuning
    float yaw_speed_deg_per_s   = 70.0f;
    float pitch_speed_deg_per_s = 55.0f;
    float zoom_speed_per_s      = 6.0f;
    float fast_multiplier       = 2.25f;

    // Internal mouse state
    bool   rmb_was_down = false;
    double last_x = 0.0;
    double last_y = 0.0;

    void setTarget(Vec3f t) { target = t; }

    // Updates yaw/pitch/dist based on current input. If fixed_target is provided,
    // it becomes the current target (useful when you want a fixed orbit pivot).
    void updateFromInput(GLFWwindow* window,
                         const ImGuiIO& io,
                         double wall_dt,
                         Vec3f fixed_target,
                         bool lock_target = true)
    {
        if (lock_target) {
            target = fixed_target;
        }

        // Mouse orbit/zoom only when ImGui does not want the mouse
        const bool imgui_captures_mouse = io.WantCaptureMouse;
        const bool rmb_down = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS);

        if (!imgui_captures_mouse && rmb_down) {
            double mx = 0.0, my = 0.0;
            glfwGetCursorPos(window, &mx, &my);

            if (!rmb_was_down) {
                // Start tracking on first frame of RMB down to avoid jump
                last_x = mx;
                last_y = my;
                rmb_was_down = true;
            } else {
                const double dx = mx - last_x;
                const double dy = my - last_y;

                yaw_deg   += (float)(dx * orbit_sens_deg_per_px);
                pitch_deg -= (float)(dy * orbit_sens_deg_per_px);

                last_x = mx;
                last_y = my;
            }
        } else {
            rmb_was_down = false;
        }

        if (!imgui_captures_mouse && io.MouseWheel != 0.0f) {
            // Exponential zoom is stable across wheel hardware
            const float k = std::pow(zoom_step, io.MouseWheel);
            dist *= k;
        }

        // Angle sanitation
        pitch_deg = std::clamp(pitch_deg, -85.0f, 85.0f);
        if (yaw_deg > 360.0f || yaw_deg < -360.0f) {
            yaw_deg = std::fmod(yaw_deg, 360.0f);
        }

        // Keyboard orbit/zoom only when ImGui does not want the keyboard
        const bool imgui_captures_kb = io.WantCaptureKeyboard || io.WantTextInput;
        if (!imgui_captures_kb) {
            const float dt = (float)wall_dt;
            const bool fast = ImGui::IsKeyDown(ImGuiKey_LeftShift) || ImGui::IsKeyDown(ImGuiKey_RightShift);
            const float mul = fast ? fast_multiplier : 1.0f;

            if (ImGui::IsKeyDown(ImGuiKey_LeftArrow))  yaw_deg   -= yaw_speed_deg_per_s   * dt * mul;
            if (ImGui::IsKeyDown(ImGuiKey_RightArrow)) yaw_deg   += yaw_speed_deg_per_s   * dt * mul;
            if (ImGui::IsKeyDown(ImGuiKey_UpArrow))    pitch_deg += pitch_speed_deg_per_s * dt * mul;
            if (ImGui::IsKeyDown(ImGuiKey_DownArrow))  pitch_deg -= pitch_speed_deg_per_s * dt * mul;

            pitch_deg = std::clamp(pitch_deg, -85.0f, 85.0f);

            if (ImGui::IsKeyDown(ImGuiKey_W)) dist -= zoom_speed_per_s * dt * mul;
            if (ImGui::IsKeyDown(ImGuiKey_S)) dist += zoom_speed_per_s * dt * mul;
        }

        dist = std::clamp(dist, min_dist, max_dist);
    }

    // Returns current eye position in world space.
    Vec3f eye() const {
        const float yaw   = yaw_deg   * 3.1415926535f / 180.0f;
        const float pitch = pitch_deg * 3.1415926535f / 180.0f;

        return v3(
            target.x + dist * std::cos(pitch) * std::sin(yaw),
            target.y + dist * std::sin(pitch),
            target.z + dist * std::cos(pitch) * std::cos(yaw)
        );
    }
};
