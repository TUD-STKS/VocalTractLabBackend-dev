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

#include "VocalTractLabBackend/XmlHelper.h"
#include <cstdio>

#include <regex>
#include <sstream>

// ****************************************************************************
// ****************************************************************************

XmlNode *XmlHelper::getChildNode(XmlNode *node, const char *childName, int index)
{
  if ((node == NULL) || (childName == NULL))
  { 
    throw std::string("Invalid parameters for getChildNode(...).");
  }

  XmlNode *childNode = node->getChildElement(childName, index);

  if (childNode == NULL)
  {
    char st[512];
    sprintf(st, "The child element <%s> of the node <%s> at position %d"
      " does not exist!", childName, node->name.c_str(), index);
    throw std::string(st);
  }

  return childNode;
}


// ****************************************************************************
// ****************************************************************************

void XmlHelper::readAttribute(XmlNode *node, const char *attrName, double &attrValue)
{
  if ((node == NULL) || (attrName == NULL))
  {
    throw std::string("Invalid parameters for readAttribute(...).");
  }

  if (node->hasAttribute(attrName) == false)
  {
    throw std::string("The attribute '") + attrName + "' for the element <" + 
      node->name + "> does not exist!"; 
  }

  attrValue = node->getAttributeDouble(attrName);
}

// ****************************************************************************
// ****************************************************************************

void XmlHelper::readAttribute(XmlNode *node, const char *attrName, int &attrValue)
{
  if ((node == NULL) || (attrName == NULL))
  {
    throw std::string("Invalid parameters for readAttribute(...).");
  }

  if (node->hasAttribute(attrName) == false)
  {
    throw std::string("The attribute '") + attrName + "' for the element <" + 
      node->name + "> does not exist!"; 
  }

  attrValue = node->getAttributeInt(attrName);
}

// ****************************************************************************
// ****************************************************************************

void XmlHelper::readAttribute(XmlNode *node, const char *attrName, std::string &attrValue)
{
  if ((node == NULL) || (attrName == NULL))
  {
    throw std::string("Invalid parameters for readAttribute(...).");
  }

  if (node->hasAttribute(attrName) == false)
  {
    throw std::string("The attribute '") + attrName + "' for the element <" + 
      node->name + "> does not exist!"; 
  }

  attrValue = node->getAttributeString(attrName);
}

// ****************************************************************************
// A very basic XML formatter. Currently only does indentation and requires
// closing tags to be on their own line.
// ****************************************************************************

std::string XmlHelper::formatXmlString(const std::string& unformattedString)
{
    // Get rid of current indentation
    std::stringstream ss(std::regex_replace(unformattedString, std::regex(R"(\n\s+)"), "\n"));

    // Split into lines
	std::string line;
	std::vector<std::string> lines;
    while (std::getline(ss, line))
    {
        lines.push_back(line);
    }

    // Indent properly
    unsigned indent{ 0 };
    std::ostringstream os;
    for (auto& line : lines)
    {
        // If the line is a closing tag, reduce indent by one step
        if (line.substr(0, 2) == "</")
        {
            indent -= 2;
        }
        // Indent
        line = std::string(indent, ' ') + line;

    	// Add line to string
        os << line << std::endl;

        // If the current line is not a self-closing or closing tag, increase indentation by one step
	    if (line.substr(line.size() - 2, 2) != "/>" && line.substr(indent, 2) != "</")
	    {
            indent += 2;
	    }
    }

    return os.str();
}

// ****************************************************************************

