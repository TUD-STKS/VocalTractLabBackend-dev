#include "Speaker.h"

#include "GlottisFactory.h"
#include "XmlNode.h"


Speaker::Speaker(const std::string& path) : Speaker()
{
    this->read(path);
}

size_t Speaker::addGlottisModel(Glottis& newModel)
{
	glottisModels.push_back(&newModel);
	return glottisModels.size();
}

std::pair<std::vector<Glottis*>, size_t> Speaker::getGlottisModels() const
{
	return make_pair(glottisModels, getSelectedGlottis());
}

size_t Speaker::getSelectedGlottis() const
{
	return selectedGlottis;
}

void Speaker::setGlottisModels(const std::vector<Glottis*>& newModels)
{
	glottisModels = newModels;
}

void Speaker::setSelectedGlottis(size_t index)
{
	selectedGlottis = index;
}

void Speaker::setVocalTract(VocalTract* newModel)
{
	vocalTract = newModel;
}

VocalTract* Speaker::getVocalTract() const
{
	return vocalTract;
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

	// This may be overwritten later.
	selectedGlottis = 0;
    XmlNode* glottisModelsNode = rootNode->getChildElement("glottis_models");
    if (glottisModelsNode != nullptr)
    {
        glottisModels.clear();
	    for (const auto& node : glottisModelsNode->childElement)
        {
            // Create a new glottis from the XML item and add it to the list
            auto* newGlottis = GlottisFactory::makeGlottis(*node);
            glottisModels.push_back(newGlottis);

            // Check if the current model is the default model
            if (node->getAttributeInt("selected") == 1)
            {
                selectedGlottis = glottisModels.size() - 1;
            }
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
