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

#include "VocalTractLabBackend/GlottisFactory.h"

#include "VocalTractLabBackend/GeometricGlottis.h"
#include "VocalTractLabBackend/TriangularGlottis.h"
#include "VocalTractLabBackend/TwoMassModel.h"

#include <string>

Glottis* GlottisFactory::makeGlottis(XmlNode& xml)
{
	return makeGlottis(xml.getAttributeString("type"), xml);
}

Glottis* GlottisFactory::makeGlottis(GlottisModel type)
{
	switch(type)
	{
		case GEOMETRIC_GLOTTIS:
			return new GeometricGlottis();
		case TWO_MASS_MODEL: 
			return new TwoMassModel();
		case TRIANGULAR_GLOTTIS: 
			return new TriangularGlottis();
		default: 
			throw std::invalid_argument("[GlottisFactory::makeGlottis()] Invalid glottis type requested: " + std::to_string(type));
	}
}

Glottis* GlottisFactory::makeGlottis(GlottisModel type, XmlNode& xml)
{
	switch (type)
	{
	case GEOMETRIC_GLOTTIS:
		return new GeometricGlottis(xml);
	case TWO_MASS_MODEL:
		return new TwoMassModel(xml);
	case TRIANGULAR_GLOTTIS:
		return new TriangularGlottis(xml);
	default:
		throw std::invalid_argument("[GlottisFactory::makeGlottis()] Invalid glottis type requested: " + std::to_string(type));
	}
}

Glottis* GlottisFactory::makeGlottis(const std::string& type)
{
	const auto enumIt = glottis_name_to_enum.find(type);
	if (enumIt == glottis_name_to_enum.end()) throw std::invalid_argument("[GlottisFactory::makeGlottis()] Invalid glottis name: " + type);

	return makeGlottis(enumIt->second);
}

Glottis* GlottisFactory::makeGlottis(const std::string& type, XmlNode& xml)
{
	const auto enumIt = glottis_name_to_enum.find(type);
	if (enumIt == glottis_name_to_enum.end()) throw std::invalid_argument("[GlottisFactory::makeGlottis()] Invalid glottis name: " + type);

	return makeGlottis(enumIt->second, xml);
}
