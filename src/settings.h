// OpenCSG - library for image-based CSG rendering for OpenGL
// Copyright (C) 2002-2016, Florian Kirsch,
// Hasso-Plattner-Institute at the University of Potsdam, Germany
//
// This library is free software; you can redistribute it and/or 
// modify it under the terms of the GNU General Public License, 
// Version 2, as published by the Free Software Foundation.
// As a special exception, you have permission to link this library
// with the CGAL library and distribute executables.
//
// This library is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License 
// along with this program; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

//
// settings.h 
//
// global OpenCSG settings
//

#ifndef __OpenCSG__settings_h__
#define __OpenCSG__settings_h__

#include "opencsgConfig.h"
#include <opencsg.h>

namespace OpenCSG {

    /// Initializes the set of options. No need to call this from outside,
    /// because invoked from setOptioni() and getOptioni() when they are first called.
    void initOptions();
    /// redeclared from opencsg.h
    void setOption(OptionType option, int newSetting);
    /// redeclared from opencsg.h
    int  getOption(OptionType option);

} // namespace OpenCSG

#endif // __OpenCSG__settings_h__

