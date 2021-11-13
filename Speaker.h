#ifndef SPEAKERFILE_H
#define SPEAKERFILE_H
#include <string>

#include "Glottis.h"
#include "VocalTract.h"

class Speaker
{
public:
	Speaker() = default;
	Speaker(VocalTract* vocalTract, const std::vector<Glottis*>& glottisModels);
	Speaker(const std::string& path);

public:
	// Getter and Setter
	size_t addGlottisModel(Glottis& newModel);
	std::pair<std::vector<Glottis*>, size_t> getGlottisModels() const;
	size_t getSelectedGlottis() const;
	void setGlottisModels(const std::vector<Glottis*>& newModels);
	void setSelectedGlottis(size_t index);

	void setVocalTract(VocalTract* newModel);
	VocalTract* getVocalTract() const;

	// File I/O
	void read(const std::string& path);
	void save(const std::string& path);

private:
	std::vector<Glottis*> glottisModels;
	VocalTract* vocalTract;
	size_t selectedGlottis{ 0 };
};

#endif // SPEAKERFILE_H
