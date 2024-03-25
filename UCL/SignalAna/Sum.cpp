#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

double sumNumbersInFiles(const std::string& directory) {
    double totalSum = 0.0;
    // Add filenames to exclude
    std::vector<std::string> excludedFiles = {   
    "run000347-ana_vertexfit.root",
    "run000365-ana_vertexfit.root",
    "run000402-ana_vertexfit.root",
    "run000435-ana_vertexfit.root",
    "run000495-ana_vertexfit.root",
    "run000820-ana_vertexfit.root",
    "run000822-ana_vertexfit.root",
    "run000840-ana_vertexfit.root",
    "run000848-ana_vertexfit.root",
    "run000873-ana_vertexfit.root",
    "run000906-ana_vertexfit.root",
    "run000916-ana_vertexfit.root",
    "run000917-ana_vertexfit.root",
    "run000926-ana_vertexfit.root",
    "run000930-ana_vertexfit.root",
    "run000933-ana_vertexfit.root",
    "run000941-ana_vertexfit.root",
    "run000943-ana_vertexfit.root",
    "run000950-ana_vertexfit.root",
    "run000953-ana_vertexfit.root",
    "run000957-ana_vertexfit.root",
    "run000960-ana_vertexfit.root",
    "run000967-ana_vertexfit.root",
    "run000968-ana_vertexfit.root",
    "run000971-ana_vertexfit.root",
    "run000984-ana_vertexfit.root",
    "run000991-ana_vertexfit.root",
    "run000993-ana_vertexfit.root",
    "run000995-ana_vertexfit.root",
    "run000996-ana_vertexfit.root",
    }; // Add filenames to exclude here
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.path().extension() == ".txt") {
            bool excludeFile = false;
            for (const std::string& excludedFile : excludedFiles) {
                if (entry.path().filename() == excludedFile) {
                    excludeFile = true;
                    break;
                }
            }
            if (excludeFile)
                continue; // Skip this file

            std::ifstream file(entry.path());
            if (file.is_open()) {
                double number;
                if (file >> number) {
                    totalSum += number;
                } else {
                    std::cerr << "Error: Cannot read number from file " << entry.path() << std::endl;
                }
            } else {
                std::cerr << "Error: Cannot open file " << entry.path() << std::endl;
            }
        }
    }
    return totalSum;
}

int main() {
    std::string directory = "/unix/muons/mu3e/Samples/4.4.3/Signal";
    double totalSum = sumNumbersInFiles(directory);
    std::cout << "Total sum of numbers in all .txt files (excluding specified files): " << totalSum << std::endl;
    return 0;
}