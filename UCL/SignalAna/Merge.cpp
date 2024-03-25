#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TSystemDirectory.h"

#include <iostream>
#include <TFile.h>
#include <TH1.h>

void MergeHistograms(const std::vector<std::string>& fileNames, const std::string& outputFileName) {
    // Create an output ROOT file
    TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");

    // Loop over input files
    for (const auto& fileName : fileNames) {
        // Open input ROOT file
        TFile* inputFile = TFile::Open(fileName.c_str());
        if (!inputFile) {
            std::cerr << "Error: Cannot open input file " << fileName << std::endl;
            continue;
        }

        // Loop over histograms in the input file
        TIter nextkey(inputFile->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)nextkey())) {
            TObject* obj = key->ReadObj();
            if (obj->IsA()->InheritsFrom("TH1")) {
                // Merge histogram into output file
                TH1* hist = (TH1*)obj;
                TH1* outputHist = dynamic_cast<TH1*>(outputFile->Get(hist->GetName()));
                if (!outputHist) {
                    outputHist = (TH1*)hist->Clone();
                    outputHist->SetDirectory(outputFile);
                } else {
                    outputHist->Add(hist);
                }
            }
            delete obj;
        }

        // Close input ROOT file
        inputFile->Close();
    }

    // Write and close output ROOT file
    outputFile->Write();
    outputFile->Close();
    delete outputFile;
}

int main() {
    std::vector<std::string> fileNames = {
        "signal_0_to_9.root",
        "signal_10_to_19.root",
        "signal_20_to_29.root",
        "signal_30_to_39.root",
        "signal_40_to_49.root",
        "signal_50_to_59.root",
        "signal_60_to_69.root",
        "signal_70_to_79.root",
        "signal_80_to_89.root",
        "signal_90_to_99.root",
        "signal_100_to_199.root",
        "signal_200_to_299.root",
        "signal_300_to_399.root",
        "signal_400_to_499.root",
        "signal_500_to_599.root",
        "signal_600_to_699.root",
        "signal_700_to_799.root",
        "signal_800_to_899.root",
        "signal_900_to_999.root",
        
        
        
        

    
        
        
        // Add more file names as needed
    };
    std::string outputFileName = "merged_signal_weightcut.root";
    MergeHistograms(fileNames, outputFileName);
    return 0;
}

