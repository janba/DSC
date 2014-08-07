//
//  Deformabel Simplicial Complex (DSC) method
//  Copyright (C) 2013  Technical University of Denmark
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  See licence.txt for a copy of the GNU General Public License.

#include "user_interface.h"
#include <iostream>
using namespace std;

void printMenu(){
    cout << "--------- DSC Demo ---------\n";
    cout << " *** SELECT PARAMETERS ***  \n"<<
            ",:         Decreases discretization by 0.5 to a minimum of 1.\n"<<
            ".:         Increases discretization by 0.5 to a maximum of 100.\n"<<
            "-:         Decreases velocity by 1 to a minimum of 1.\n"<<
            "+:         Increases velocity by 1 to a maximum of 100.\n"<<
            ">:         Decreases accuracy by 1 to a minimum of 1.\n"<<
            "<:         Increases accuracy by 1 to a maximum of 100.\n"<<
            "\n"<<
            "*** SELECT MOTION ***\n"<<
            "1:         Selects rotation.\n"<<
            "2:         Selects smoothing.\n"<<
            "3:         Selects expansion.\n"<<
            "\n"<<
            "*** START/STOP MOTION ***\n"<<
            "SPACE:     Starts/pauses the current motion.\n"<<
            "0:         Stops the current motion.\n"<<
            "ESCAPE:    Stops the current motion and exits the application\n"<<
            "m:         Moves the interface vertices one time step according to the current velocity function.\n"<<
            "\n"<<
            "*** MISCELLANEOUS ***\n"<<
            "r:         Reloads the model.\n"<<
            "t:         Performs a test on the current velocity function.\n"<<
            "s:         Takes a screen shot.\n"<<
            "e:         Export the simplicial complex to a .dsc file.\n"<<
            "i:         Export the surface mesh to a .obj file.\n"<<
            "w:         Switch wireframe rendering on and off.\n"<<
            "TAB:       Switches the display type (Surface, wireframe, edges, etc.).\n"<<endl;
}

int main(int argc, char** argv)
{
    UI ui(argc, argv);
    printMenu();
    glutMainLoop();
    return 0;
}
