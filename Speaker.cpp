#include "Speaker.h"

#include <fstream>

#include "GlottisFactory.h"
#include "XmlNode.h"


Speaker::Speaker(VocalTract* vocalTract, const std::vector<Glottis*>& glottisModels) : glottisModels(glottisModels), vocalTract(vocalTract)
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

void Speaker::setGlottisModels(const std::vector<Glottis*>& newModels)
{
	glottisModels = newModels;
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
    os << "<speaker>" << std::endl;

    os << *vocalTract;

    os << "  <glottis_models" << std::endl;

    for (const auto glottisModel : glottisModels)
    {
        os << *glottisModel;
    }

    os << "  </glottis_models>" << std::endl;
    os << "</speaker>" << std::endl;

    return os;
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
