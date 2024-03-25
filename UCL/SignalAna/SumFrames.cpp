#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <TLegend.h>
#include <TStyle.h>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TSystem.h> // For directory listing
#include <string>
#include <vector>
#include <TSystemDirectory.h>


double sumLeafValues(const std::string& directory, const char* treeName, const char* leafName) {
    // Get list of ROOT files in the directory
    std::vector<std::string> files;
    TSystemDirectory dir(directory.c_str(), directory.c_str());
    TList *fileList = dir.GetListOfFiles();
    if (fileList) {
        TSystemFile *file;
        TIter next(fileList);
        while ((file=(TSystemFile*)next())) {
            const char* fname = file->GetName();
            if (!file->IsDirectory() && TString(fname).EndsWith(".root"))
                files.push_back(directory + "/" + fname);
        }
        delete fileList;
    }

    // Variable to hold the total sum
    double totalSum = 0.0;

    // Loop over files in the directory
    for (const auto& file : files) {
        // Open the ROOT file
        TFile *rootFile = TFile::Open(file.c_str());
        if (!rootFile || rootFile->IsZombie()) {
            std::cerr << "Error: Could not open ROOT file or file is corrupted: " << file << std::endl;
            continue;
        }

        // Get the TTree
        TTree *tree = nullptr;
        rootFile->GetObject(treeName, tree);
        if (!tree) {
            std::cerr << "Error: Could not retrieve TTree from ROOT file: " << file << std::endl;
            rootFile->Close();
            continue;
        }

        // Access the TBranch corresponding to the leaf
        TBranch *branch = tree->GetBranch(leafName);
        if (!branch) {
            std::cerr << "Error: Could not retrieve TBranch from TTree: " << file << std::endl;
            rootFile->Close();
            continue;
        }

        // Create a buffer to hold the data
        double leafValue = 0.0;

        // Connect the buffer to the TBranch
        branch->SetAddress(&leafValue);

        // Loop over entries in the tree
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
            
            // Fill the branch buffer with data for this entry
            branch->GetEntry(iEntry);
            
            // Accumulate the leaf value
            totalSum += leafValue;
            
        }

        // Close the ROOT file
        rootFile->Close();
    }

    return totalSum;
}

int main() {
    const char* directory = "/unix/muons/mu3e/Samples/4.4.3/ICmichel";
    const char* treeName = "vertex";
    const char* leafName = "nEntries";

    double totalSum = sumLeafValues(directory, treeName, leafName);
    std::cout << "Total sum of values in leaf \"" << leafName << "\" across all files in directory \"" << directory << "\": " << totalSum << std::endl;

    return 0;
}