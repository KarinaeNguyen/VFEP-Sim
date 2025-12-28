import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Configuration for the Innovation Pitch
DATA_PATH = '../data/outputs/aerodynamic_ml_set.csv'

def launch_dashboard():
    if not os.path.exists(DATA_PATH):
        print(f"CRITICAL ERROR: {DATA_PATH} not found.")
        print("Please run the C++ VFEP_Sim.exe first to generate the Physics Data.")
        return

    # Load the high-fidelity dataset
    df = pd.read_csv(DATA_PATH)

    # Initialize a professional 3-panel dashboard
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 1, hspace=0.4)
    
    # --- PANEL 1: THERMAL PERFORMANCE ---
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(df['Time'], df['Temp'], color='#e74c3c', linewidth=2, label='Core Temperature (°C)')
    ax1.axhline(y=70, color='gray', linestyle='--', label='AI Trigger Threshold')
    ax1.set_ylabel('Temperature (°C)')
    ax1.set_title('THERMAL ANALYSIS: Fire Growth vs. Suppression', fontsize=12, fontweight='bold')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)

    # --- PANEL 2: AERODYNAMICS & EFFICIENCY ---
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.fill_between(df['Time'], df['Hit_Efficiency'] * 100, color='#2ecc71', alpha=0.2, label='Target Hit Rate')
    ax2.plot(df['Time'], df['Wind_Speed'] * 5, color='#3498db', linestyle=':', label='Wind Turbulence (Scaled)')
    ax2.set_ylabel('Efficiency % / Wind')
    ax2.set_title('AERODYNAMIC INTEGRITY: Agent Delivery Efficiency under Wind Stress', fontsize=12, fontweight='bold')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)

    # --- PANEL 3: AI RPM RESPONSE ---
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.step(df['Time'], df['RPM'], color='#8e44ad', linewidth=2, label='Variable RPM (AI Action)')
    ax3.set_ylabel('Motor RPM')
    ax3.set_xlabel('Simulation Time (Seconds)')
    ax3.set_title('AI DECISION LOG: Dynamic Pressure Compensation', fontsize=12, fontweight='bold')
    ax3.legend(loc='upper right')
    ax3.grid(True, alpha=0.3)

    plt.suptitle("VFEP GUARDIAN: INDUSTRIAL SUPPRESSION ANALYTICS", fontsize=16, fontweight='bold')
    print("--- SUCCESS: DASHBOARD GENERATED FOR PITCH ---")
    plt.show()

if __name__ == "__main__":
    launch_dashboard()
