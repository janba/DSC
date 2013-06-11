#ifndef OPENTISSUE_CONFIGURATIOM_H
#define OPENTISSUE_CONFIGURATIOM_H
//
// OpenTissue, A toolbox for physical based simulation and animation.
// Copyright (C) 2007 Department of Computer Science, University of Copenhagen
//
#if (_MSC_VER >= 1200)
# pragma once
# pragma warning(default: 56 61 62 191 263 264 265 287 289 296 347 529 686)
# pragma warning(disable: 503)
#endif

#ifdef WIN32
#  define WIN32_LEAN_AND_MEAN
#  define _USE_MATH_DEFINES
#  define NOMINMAX
#  include <windows.h>
#  undef WIN32_LEAN_AND_MEAN
#  undef NOMINMAX
#endif


/**
 * OpenTissue Version string
 */
#define OPENTISSUE_VERSION  "OpenTissue version 0.993"

#include <string>

/**
 * OpenTissue Path.
 * This is the path where OpenTissue was copied onto ones system. It can be used to locate shader programs or data resources.
 */
std::string const opentissue_path = "C:/DTU/OpenTissue/";

/**
 * OpenTissue Version String.
 * This string value can be used by end users for compatibility testing.
 */
std::string const opentissue_version = OPENTISSUE_VERSION;

//OPENTISSUE_CONFIGURATIOM_H
#endif
