//
// Created by Morten Nobel-Jørgensen on 12/12/14.
// Copyright (c) 2014 Asger Nyman Christiansen. All rights reserved.
//

#include "DSC_Suite.h"

int t1(void) {
    int* p = NULL;
    TINYTEST_ASSERT(!p);
    TINYTEST_ASSERT(!printf(""));
    return 1; // Always return a value different than 0 at test end.
}