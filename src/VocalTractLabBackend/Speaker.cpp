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

#include "VocalTractLabBackend/Speaker.h"

#include <fstream>
#include <sstream>

#include "VocalTractLabBackend/GlottisFactory.h"
#include "VocalTractLabBackend/XmlHelper.h"
#include "VocalTractLabBackend/XmlNode.h"


Speaker::Speaker(VocalTract* vocalTract, const std::vector<Glottis*>& glottisModels, size_t selectedGlottis) :
	glottisModels(glottisModels), selectedGlottis(selectedGlottis), vocalTract(vocalTract)
{

}

Speaker::Speaker(const std::string& path) : Speaker()
{
    this->read(path);
}

size_t Speaker::addGlottisModel(Glottis& newModel)
{
	glottisModels.push_back(&newModel);
	return glottisModels.size();
}

std::vector<Glottis*> Speaker::getGlottisModels() const
{
	return glottisModels;
}

size_t Speaker::getSelectedGlottis() const
{
	return selectedGlottis;
}

void Speaker::setGlottisModels(const std::vector<Glottis*>& newModels)
{
	glottisModels = newModels;
}

void Speaker::setSelectedGlottis(size_t idx)
{
	selectedGlottis = idx;
}

void Speaker::setVocalTract(VocalTract* newModel)
{
	vocalTract = newModel;
}

VocalTract* Speaker::getVocalTract() const
{
	return vocalTract;
}

std::ostream& Speaker::save(std::ostream& os) const
{
    std::ostringstream xml;
    xml << "<speaker>" << std::endl;

    xml << *vocalTract;

    xml << "<glottis_models" << std::endl;

    for (size_t i = 0; i < glottisModels.size(); ++i)
    {
        if (i != selectedGlottis)
        {
            xml << *glottisModels[i];
        }
        else
        {
            glottisModels[i]->writeToXml(xml, 0, true);
        }
			
    }

    xml << "</glottis_models>" << std::endl;
    xml << "</speaker>" << std::endl;

    return os << XmlHelper::formatXmlString(xml.str());
}

void Speaker::read(const std::string& path)
{
	// ****************************************************************
	// Load the XML data from the speaker file.
	// ****************************************************************

	std::vector<XmlError> xmlErrors;
	XmlNode* rootNode = xmlParseFile(path, "speaker", &xmlErrors);
	if (rootNode == nullptr)
	{
		xmlPrintErrors(xmlErrors);
		throw XmlException("[SpeakerFile::read()] Could not parse file " + path);
	}

	// ****************************************************************
	// Load the data for the glottis models.
	// ****************************************************************

    XmlNode* glottisModelsNode = rootNode->getChildElement("glottis_models");
    if (glottisModelsNode != nullptr)
    {
        glottisModels.clear();
	    for (const auto& node : glottisModelsNode->childElement)
        {
            // Create a new glottis from the XML item and add it to the list
            auto* newGlottis = GlottisFactory::makeGlottis(*node);
            if (node->getAttributeInt("selected") == 1)
            {
                this->setSelectedGlottis(glottisModels.size());
            }                
            glottisModels.push_back(newGlottis);
        }
    }
    else
    {
        printf("Warning: No glottis model data found in the speaker file %s!\n", path.c_str());
    }


    // ****************************************************************
	// Load the vocal tract anatomy and vocal tract shapes.
	// ****************************************************************

    try
    {
        XmlNode* vocalTractNode = rootNode->getChildElement("vocal_tract_model");
        vocalTract = new VocalTract(*vocalTractNode);
        vocalTract->calculateAll();
    }
    catch (std::string st)
    {
        printf("%s\n", st.c_str());
        printf("Error reading the anatomy data from %s.\n", path.c_str());
        return;
    }

    // Free the memory of the XML tree !
    delete rootNode;
}

void Speaker::save(const std::string& path) const
{
    std::ofstream os(path);
    if (!os)
    {
        throw std::runtime_error("[Speaker::save()] Could not open " + path + "for writing!");
    }

    os << *this;

    os.close();
}

std::ostream& operator<<(std::ostream& os, const Speaker& obj)
{
	return obj.save(os);
}
