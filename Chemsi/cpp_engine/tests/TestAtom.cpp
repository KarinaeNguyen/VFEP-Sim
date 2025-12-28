#include <iostream>
#include <cassert>   // This library gives us the 'assert' command
#include "../include/Atom.h" // We need to go up one folder (..) to find include

void testAtomCreation() {
    // 1. Setup: Create a Nitrogen atom
    Atom a(1, "N", 14.007, 0.0, 0.0, 0.0);

    // 2. Verification (The "Test")
    // If these are true, nothing happens.
    // If these are false, the program CRASHES immediately (which is good!)
    assert(a.id == 1);
    assert(a.type == "N");
    assert(a.mass == 14.007);
    
    // Check initial position
    assert(a.x == 0.0);

    std::cout << "✅ Atom Creation Test Passed!" << std::endl;
}

void testAtomMovement() {
    Atom a(1, "H", 1.008, 0.0, 0.0, 0.0);
    
    // Manually change position (simulating movement)
    a.x = 5.0;
    
    assert(a.x == 5.0);
    assert(a.y == 0.0); // Should still be 0

    std::cout << "✅ Atom Movement Test Passed!" << std::endl;
}

int main() {
    std::cout << "--- Running Unit Tests ---" << std::endl;
    
    testAtomCreation();
    testAtomMovement();

    std::cout << "--- All Tests Passed Successfully ---" << std::endl;
    return 0;
}
