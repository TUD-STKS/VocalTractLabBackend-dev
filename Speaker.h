#ifndef SPEAKERFILE_H
#define SPEAKERFILE_H
#include <ostream>
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
	std::vector<Glottis*> getGlottisModels() const;
	void setGlottisModels(const std::vector<Glottis*>& newModels);

	void setVocalTract(VocalTract* newModel);
	VocalTract* getVocalTract() const;

	// Stream and file I/O
	std::ostream& save(std::ostream& os) const;
	void read(const std::string& path);
	void save(const std::string& path) const;

	friend std::ostream& operator<<(std::ostream& os, const Speaker& obj)
	{
		return obj.save(os);
	}

private:
	std::vector<Glottis*> glottisModels;
	VocalTract* vocalTract;
};

#endif // SPEAKERFILE_H
