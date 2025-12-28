// vis/Visualizer.h
#pragma once

#include <vector>
#include <cstddef>
#include <algorithm>

#include "Simulation.h"

namespace vfepvis {

// Simple history buffer for plotting
struct HistoryBuffer {
    std::vector<double> t;
    std::vector<double> T_K;
    std::vector<double> HRR_W;
    std::vector<double> O2_pct;
    std::vector<double> CO2_pct;
    std::vector<double> agent_kgps;

    void clear() {
        t.clear(); T_K.clear(); HRR_W.clear(); O2_pct.clear(); CO2_pct.clear(); agent_kgps.clear();
    }

    void reserve(std::size_t n) {
        t.reserve(n); T_K.reserve(n); HRR_W.reserve(n); O2_pct.reserve(n); CO2_pct.reserve(n); agent_kgps.reserve(n);
    }

    void push(double time_s, const vfep::Observation& o, std::size_t maxN) {
        t.push_back(time_s);
        T_K.push_back(o.T_K);
        HRR_W.push_back(o.HRR_W);
        O2_pct.push_back(o.O2_volpct);
        CO2_pct.push_back(o.CO2_volpct);
        agent_kgps.push_back(o.agent_mdot_kgps);

        // Trim (simple FIFO) to keep memory bounded
        if (maxN > 0 && t.size() > maxN) {
            const std::size_t extra = t.size() - maxN;
            t.erase(t.begin(), t.begin() + (std::ptrdiff_t)extra);
            T_K.erase(T_K.begin(), T_K.begin() + (std::ptrdiff_t)extra);
            HRR_W.erase(HRR_W.begin(), HRR_W.begin() + (std::ptrdiff_t)extra);
            O2_pct.erase(O2_pct.begin(), O2_pct.begin() + (std::ptrdiff_t)extra);
            CO2_pct.erase(CO2_pct.begin(), CO2_pct.begin() + (std::ptrdiff_t)extra);
            agent_kgps.erase(agent_kgps.begin(), agent_kgps.begin() + (std::ptrdiff_t)extra);
        }
    }
};

// Visualizer state and stepping policy (no graphics code here)
class Visualizer {
public:
    Visualizer();

    void reset();
    void toggleRun() { running_ = !running_; }
    void setRunning(bool r) { running_ = r; }
    bool isRunning() const { return running_; }

    void igniteOrIncreasePyrolysis();
    void startSuppression();

    // Advances simulation by one fixed dt step (if possible)
    void stepOnce();

    // Advances simulation only if running and not concluded
    void tick();

    // Settings
    void setDt(double dt_s);
    double dt() const { return dt_s_; }

    double timeS() const { return simTime_s_; }

    const vfep::Simulation& sim() const { return sim_; }
    vfep::Simulation& sim() { return sim_; }

    const HistoryBuffer& history() const { return hist_; }
    HistoryBuffer& history() { return hist_; }

    std::size_t maxHistory() const { return maxHistory_; }
    void setMaxHistory(std::size_t n) { maxHistory_ = n; }

private:
    void sample();

    vfep::Simulation sim_;
    bool running_ = false;

    double dt_s_ = 0.05;
    double simTime_s_ = 0.0;

    std::size_t maxHistory_ = 5000;
    HistoryBuffer hist_;
};

} // namespace vfepvis
