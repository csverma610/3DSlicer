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

#include "opencsgConfig.h"
#include "area.h"
#include "openglHelper.h"

namespace OpenCSG {

NDCVolume::operator PCArea() const {
    using OpenGL::canvasPos;

    float widthPixels  = 0.5f * (canvasPos[2] - canvasPos[0]);
    float heightPixels = 0.5f * (canvasPos[3] - canvasPos[1]);

    int iminx = static_cast<int>( (minx + 1.0f) * widthPixels );
    int imaxx = static_cast<int>( (maxx + 1.0f) * widthPixels );
    int iminy = static_cast<int>( (miny + 1.0f) * heightPixels );
    int imaxy = static_cast<int>( (maxy + 1.0f) * heightPixels );

    return PCArea(iminx, iminy, imaxx, imaxy);
}

} // namespace OpenCSG

