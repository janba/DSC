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
    cout << "Esc     Quit\n";
    cout << "0       Set velocity function\n";
    cout << "1       Set rotate function\n";
    cout << "2       Set average function\n";
    cout << "Space   Toggle pause\n";
    cout << "m       Move\n";
    cout << "r       Reload model\n";
    cout << "t       Test velocity function\n";
    cout << "tab     Change display type (INTERFACE, WIRE_FRAME, BOUNDARY, EDGES, LOW_QUALITY, UNMOVED)\n";
    cout << "s       Save screenshot\n";
    cout << "e       Export tetrahedron mesh\n";
    cout << "i       Export triangle mesh\n";
    cout << "+       Increase velocity\n";
    cout << "-       Descrease velocity\n";
    cout << ".       Increase avg edge length\n";
    cout << ",       Descrease avg edge length\n";
    cout << "<       Increase accuracy\n";
    cout << ">       Descrease accuracy\n";



}

int main(int argc, char** argv)
{
    UI ui(argc, argv);
    printMenu();
    glutMainLoop();
    return 0;
}
