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

#ifndef __XML_EXCEPTIONS_H__
#define __XML_EXCEPTIONS_H__

#include <string>

using namespace std;


// ****************************************************************************
// Classes for custom XML rlated exceptions.
// ****************************************************************************

class XmlException : public std::exception
{
public:
    XmlException() = default;
    XmlException(const std::string& msg) : message(msg) {}
    const char* what() const throw () override;

private:
    std::string message{ "An XML exception occurred!" };
};

inline const char* XmlException::what() const throw ()
{
    return message.c_str();
}

// ****************************************************************************

class XmlAttributeNotFound : public XmlException
{
public:
    XmlAttributeNotFound(const std::string& attribute) : message(attribute) { }
    const char* what() const throw () override;

private:
    std::string message;
};

inline const char* XmlAttributeNotFound::what() const throw ()
{
    return (std::string("Tag '") + message + std::string("' could not be found in node!")).c_str();
}

// ****************************************************************************

class XmlBadAttribute : public XmlException
{
public:
    const char* what() const throw () override;
};

inline const char* XmlBadAttribute::what() const throw ()
{
    return "Attribute could not be converted to expected type!";
}

#endif
