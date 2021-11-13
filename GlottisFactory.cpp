#include "GlottisFactory.h"

#include "GeometricGlottis.h"
#include "TriangularGlottis.h"
#include "TwoMassModel.h"

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
			throw std::invalid_argument("[GlottisFactory::getGlottis()] Invalid glottis type requested: " + std::to_string(type));
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
		throw std::invalid_argument("[GlottisFactory::getGlottis()] Invalid glottis type requested: " + std::to_string(type));
	}
}

Glottis* GlottisFactory::makeGlottis(const std::string& type)
{
	const auto enumIt = string_to_enum.find(type);
	if (enumIt == string_to_enum.end()) throw std::invalid_argument("[GlottisFactory::getGlottis()] Invalid glottis name: " + type);

	return makeGlottis(enumIt->second);
}

Glottis* GlottisFactory::makeGlottis(const std::string& type, XmlNode& xml)
{
	const auto enumIt = string_to_enum.find(type);
	if (enumIt == string_to_enum.end()) throw std::invalid_argument("[GlottisFactory::getGlottis()] Invalid glottis name: " + type);

	return makeGlottis(enumIt->second, xml);
}
