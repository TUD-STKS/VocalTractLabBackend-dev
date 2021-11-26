// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2020, Peter Birkholz, Dresden, Germany
// www.vocaltractlab.de
// author: Simon Stone
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

#ifndef GLOTTISFACTORY_H
#define GLOTTISFACTORY_H

#include <map>
#include <string>

#include "Glottis.h"

/// @brief Types of available Glottis objects
enum GlottisModel
{
	GEOMETRIC_GLOTTIS,
	TWO_MASS_MODEL,
	TRIANGULAR_GLOTTIS,
	NUM_GLOTTIS_MODELS
};

/// @brief Contains factory functions to make Glottis-derived objects
///	These functions can be helpful in dynamic object construction, e.g. when
///	reading a series of sibling Glottis-derived objects from a file/stream
namespace GlottisFactory  // Factory design pattern
{
	/// @brief A map to convert a string representation of a glottis type
	///	to the corresponding enum value.
	const map<std::string, GlottisModel> glottis_name_to_enum
	{
	{"Geometric glottis", GEOMETRIC_GLOTTIS},
	{"Two-mass model", TWO_MASS_MODEL},
	{"Triangular glottis", TRIANGULAR_GLOTTIS}
	};

	/// @brief Create a Glottis-derived object based on an XML representation.
	/// @param xml An XmlNode object containing the <glottis_model> tag
	/// @return A Glottis-derived object initialized with the data
	/// from the XML representation. The type depends on the "type"
	/// attribute in the XML representation.
	Glottis* makeGlottis(XmlNode& xml);
	/// @brief Create a default Glottis-derived object
	/// @param type The type of Glottis to create (see GlottisModel for
	/// options). Uses the default constructor of the respective object.
	/// @return A default-initialized Glottis-derived object of type "type".
	Glottis* makeGlottis(GlottisModel type);
	/// @brief Create a Glottis-derived object of a specific type using
	///	and initialize it using data in an XML representation.
	/// @param type The type of Glottis to create (see GlottisModel for
	/// options).
	/// @param xml An XmlNode object containing the <glottis_model> tag
	/// @return A Glottis-derived object of type "type", initialized using
	///	the data from the XML representation.
	Glottis* makeGlottis(GlottisModel type, XmlNode& xml);
	/// @brief Create a default Glottis-derived object
	/// @param type The type of Glottis to create (see string_to_enum for
	/// options). Uses the default constructor of the respective object.
	/// @return A default-initialized Glottis-derived object of type "type".
	Glottis* makeGlottis(const std::string& type);
	/// @brief Create a Glottis-derived object of a specific type using
	///	and initialize it using data in an XML representation.
	/// @param type The type of Glottis to create (see string_to_enum for
	/// options). 
	/// @param xml An XmlNode object containing the <glottis_model> tag
	/// @return A Glottis-derived object of type "type", initialized using
	///	the data from the XML representation.
	Glottis* makeGlottis(const std::string& type, XmlNode& xml);
};

#endif // GLOTTISFACTORY_H
