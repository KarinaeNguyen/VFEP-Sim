#ifndef DATALOGGER_H
#define DATALOGGER_H

#include <iostream>
#include <fstream>
#include <string>
#include "Projectile.h"

class DataLogger {
private:
    std::ofstream file; 

public:
    void startLog(std::string filename) {
        file.open(filename);
        file << "Time,Type,Mass,Temp,X,Y,Z\n"; 
        // CLEAN LOG MESSAGE (No Emojis)
        std::cout << "--- [LOG] Data Recording Started: " << filename << " ---" << std::endl;
    }

    void logFrame(double time, Projectile& p, double temperature) {
        if (file.is_open()) {
            file << time << "," 
                 << p.payload.name << "," 
                 << p.mass << "," 
                 << temperature << ","
                 << p.x << "," << p.y << "," << p.z << "\n";
        }
    }

    void endLog() {
        file.close();
        std::cout << "--- [LOG] Data Saved Successfully ---" << std::endl;
    }
};

#endif