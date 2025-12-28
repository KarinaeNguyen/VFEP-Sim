// vis/Visualizer.cpp
#include "Visualizer.h"

#include <cmath>
#include <algorithm>

namespace vfepvis {

namespace {
static bool isFinitePositive(double x) {
    return std::isfinite(x) && x > 0.0;
}
} // namespace

Visualizer::Visualizer() {
    hist_.reserve(5000);
    reset();
}

void Visualizer::setDt(double dt_s) {
    // Clamp dt to a reasonable range for stability + UI responsiveness
    if (!std::isfinite(dt_s)) return;
    dt_s_ = std::clamp(dt_s, 0.001, 0.5);
}

void Visualizer::reset() {
    sim_.resetToDataCenterRackScenario();
    running_ = false;
    simTime_s_ = 0.0;
    hist_.clear();
    sample();
}

void Visualizer::igniteOrIncreasePyrolysis() {
    sim_.commandIgniteOrIncreasePyrolysis();
    sample();
}

void Visualizer::startSuppression() {
    sim_.commandStartSuppression();
    sample();
}

void Visualizer::sample() {
    hist_.push(simTime_s_, sim_.observe(), maxHistory_);
}

void Visualizer::stepOnce() {
    if (!isFinitePositive(dt_s_)) return;
    if (sim_.isConcluded()) return;

    sim_.step(dt_s_);
    simTime_s_ += dt_s_;
    sample();
}

void Visualizer::tick() {
    if (!running_) return;
    stepOnce();
}

} // namespace vfepvis
