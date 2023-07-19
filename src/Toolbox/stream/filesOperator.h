#pragma once
#include <string>
#include <vector>
#include <io.h>
void getFiles(std::string path, std::vector<std::string>& files);

#include <qstring.h>
void truncateFilePath(std::string& file);

void truncateFileExtension(std::string& file);

void truncateFileName(std::string& file);
