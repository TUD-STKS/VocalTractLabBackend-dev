// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2020, Peter Birkholz, Dresden, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************

#ifndef __XML_HELPER_H__
#define __XML_HELPER_H__

#include <string>

#include "XmlNode.h"

// Disable compiler warnings in Visual Studio, because the compiler
// does not yet support the throw() statement behind function declarations.
#ifdef WIN32
#pragma warning( disable : 4290 )
#endif

// ****************************************************************************
/// Static xml support class for the xml parser (XmlNode).
// ****************************************************************************

class XmlHelper
{
public:
  static XmlNode *getChildNode(XmlNode *node, const char *childName, int index = 0) throw (std::string);

  static void readAttribute(XmlNode *node, const char *attrName, double &attrValue) throw (std::string);
  static void readAttribute(XmlNode *node, const char *attrName, int &attrValue) throw (std::string);
  static void readAttribute(XmlNode *node, const char *attrName, std::string &attrValue) throw (std::string);
};

// ****************************************************************************

#endif
