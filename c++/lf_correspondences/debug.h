/** \file debug.h

    File imported from "common" lib, use if this library is not available.

    Implements some system dependent and debugging code.
    Ancient, but should work reasonably well.

    Copyright (C) 2001 Bastian Goldluecke,
    <first name>AT<last name>.net

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __DEBUG_H_STDEXT
#define __DEBUG_H_STDEXT

#include <assert.h>
#include <stdio.h>

#include <string>
#include <vector>
#include <fstream>

// #include "../modules.h"
// #include "../defs.h"

using namespace std;

namespace miplf
{





  /****************************************************
                    MISC TOOLS
  *****************************************************/

  /// Upper case a string
  void toUpper( string &str );
  /// Make RGB value
  inline unsigned int make_rgb(int r, int g, int b)
  {
    return (0xffu << 24) | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
  }


};


/// Write output to debug output stream. The stream is automatically flushed afterwards.
/** \ingroup sysdebug */
#define TRACE(s) TRACE0(s)
#define TRACE0(s) cout << s;

/// Write output to error output stream. The stream is automatically flushed afterwards.
/** \ingroup sysdebug */
#define ERROR(s) cout <<"ERROR "<< s << endl;


#endif //__DEBUG_H_STDEXT
