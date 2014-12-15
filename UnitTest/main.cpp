//
//  main.cpp
//  UnitTest
//
//  Created by Morten Nobel-Jørgensen on 12/12/14.
//  Copyright (c) 2014 Asger Nyman Christiansen. All rights reserved.
//

#include <iostream>
#include "tinytest/tinytest.h"
#include "DSC_Suite.h"

using namespace std;

TINYTEST_START_SUITE(DSCSuite);
    TINYTEST_ADD_TEST(tetGenTest);
TINYTEST_END_SUITE();

TINYTEST_MAIN_SINGLE_SUITE(DSCSuite);