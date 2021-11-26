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

#ifndef SPEAKER_H
#define SPEAKER_H
#include <ostream>
#include <string>

#include "Glottis.h"
#include "VocalTract.h"

class Speaker
{
public:
	// Construct an empty Speaker object
	Speaker() = default;
	// Construct a Speaker from existing information
	Speaker(VocalTract* vocalTract, const std::vector<Glottis*>& glottisModels, size_t selectedGlottis = 0);
	// Construct a Speaker from a file
	Speaker(const std::string& path);


public:
	/// @brief Add a glottis model to the speaker
	/// @param newModel The model to add. Must be some sort of Glottis-derived object.
	/// @return The number of the speaker's glottis models (after the addition)
	size_t addGlottisModel(Glottis& newModel);
	/// @brief Get the speaker's glottis models
	/// @return A vector containing pointers to Glottis-derived objects. 
	std::vector<Glottis*> getGlottisModels() const;
	/// @brief Get the currently selected glottis model.
	/// @return The index in the vector of glottis models (see getGlottisModels())
	/// of the currently selected glottis model
	size_t getSelectedGlottis() const;
	/// @brief Set the speaker's glottis models.
	/// @param newModels A vector containing pointers to Glottis-derived objects.
	/// Speaker does not assume ownership of the pointers and thus never deletes them!
	void setGlottisModels(const std::vector<Glottis*>& newModels);
	/// @brief Set the currently selected glottis.
	/// @param idx The index of the currently selected glottis w.r.t. the vector of
	///	glottis models (see getGlottisModels()).
	void setSelectedGlottis(size_t idx);

	/// @brief Set the speaker's vocal tract models.
	/// @param newModel A pointer to a VocalTract object. Speaker does not assume
	/// ownership of the pointer and thus never deletes it!
	void setVocalTract(VocalTract* newModel);
	/// @brief  Get the speaker's vocal tract.
	/// @return A pointer to the speaker's vocal tract.
	VocalTract* getVocalTract() const;

	/// @brief Save the Speaker object in XML format.
	/// @param os An opened output stream.
	/// @return The output stream.
	std::ostream& save(std::ostream& os) const;
	/// @brief Read the Speaker information from an XML file.
	/// @param path Path to the XML file.
	void read(const std::string& path);
	/// @brief Save the Speaker information to an XML file.
	/// @param path Path to the XML file.
	void save(const std::string& path) const;

	/// @brief Stream output operator for XML output.
	/// @param os An already opened output stream.
	/// @param obj The Speaker object to output.
	/// @return The output stream.
	friend std::ostream& operator<<(std::ostream& os, const Speaker& obj);

private:
	/// @brief Stores the glottis models. Does not assume ownership of the
	/// stored pointers and thus does not delete them!
	std::vector<Glottis*> glottisModels;
	/// @brief Index of the currently selected model in glottisModels.
	size_t selectedGlottis{ 0 };
	/// @brief Pointer to the speaker's VocalTract. Does not assume ownership
	/// of the stored pointer and thus does not delete it!
	VocalTract* vocalTract{ nullptr };
};

#endif // SPEAKERFILE_H
