// OpenCSG - library for image-based CSG rendering for OpenGL
// Copyright (C) 2006-2016, Florian Kirsch
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
// pBufferTexture.h 
//
// RenderTexture wrapper class implementing the offscreen buffer interface
//

#ifndef __OpenCSG__pbuffer_texture_h__
#define __OpenCSG__pbuffer_texture_h__

#include "opencsgConfig.h"
#include "offscreenBuffer.h"

#ifdef OPENCSG_HAVE_PBUFFER

class RenderTexture;

namespace OpenCSG {

    namespace OpenGL {

        class PBufferTexture : public OffscreenBuffer {
        public:
            // ctor / dtor
            PBufferTexture();
            virtual ~PBufferTexture();

            /// Nothing to do here.
            virtual bool ReadCurrent();

            /// Call this once before use.  Set bShare to true to share lists, textures, 
            /// and program objects between the render texture context and the 
            /// current active GL context.
            virtual bool Initialize(int width, int height, bool shareObjects=true, bool copyContext=false);

            /// checks whether Initialize has been called before or not
            virtual bool IsInitialized() const;

            /// Change the render texture format.
            virtual bool Reset();
            /// Change the size of the render texture.
            virtual bool Resize(int width, int height);

            /// Begin drawing to the texture. (i.e. use as "output" texture)
            virtual bool BeginCapture();
            /// End drawing to the texture.
            virtual bool EndCapture();

            /// Bind the texture to the active texture unit for use as an "input" texture
            virtual void Bind() const;

            /// Enables the texture target appropriate for this render texture.
            virtual void EnableTextureTarget() const;
            /// Disables the texture target appropriate for this render texture.
            virtual void DisableTextureTarget() const;

            /// Returns the texture target this texture is bound to.
            virtual unsigned int GetTextureTarget() const;
            /// Returns the width of the offscreen buffer.
            virtual int GetWidth() const;
            /// Returns the width of the offscreen buffer.
            virtual int GetHeight() const;

            /// PBuffers have a separate rendering context (which is, in the case
            /// of OpenCSG, shared with the main context, i.e., we can use display
            /// lists uva. in both contexts
            virtual bool haveSeparateContext() const { return true; }

        protected:
            RenderTexture* r;
            const char*    s;

        private:
            PBufferTexture(const PBufferTexture&);
            PBufferTexture& operator=(const PBufferTexture&);

        };

    } // namespace OpenGL

} // namespace OpenCSG

#endif // OPENCSG_HAVE_PBUFFER

#endif // __OpenCSG__pbuffer_texture_h__
