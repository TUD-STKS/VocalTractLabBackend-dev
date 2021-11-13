#ifndef GLOTTISFACTORY_H
#define GLOTTISFACTORY_H

#include <map>
#include <string>

#include "Glottis.h"

enum GlottisModel
{
	GEOMETRIC_GLOTTIS,
	TWO_MASS_MODEL,
	TRIANGULAR_GLOTTIS,
	NUM_GLOTTIS_MODELS
};


namespace GlottisFactory  // Factory design pattern
{
	const map<std::string, GlottisModel> string_to_enum
	{
	{"Geometric glottis", GEOMETRIC_GLOTTIS},
	{"Two-mass model", TWO_MASS_MODEL},
	{"Triangular glottis", TRIANGULAR_GLOTTIS}
	};

	Glottis* makeGlottis(XmlNode& xml);
	Glottis* makeGlottis(GlottisModel type);
	Glottis* makeGlottis(GlottisModel type, XmlNode& xml);
	Glottis* makeGlottis(const std::string& type);
	Glottis* makeGlottis(const std::string& type, XmlNode& xml);
};

#endif // GLOTTISFACTORY_H
