#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TSystemDirectory.h"

void MergeHistograms(const std::string& inputDir, const std::string& outputFileName) {
    std::cout << "Input directory: " << inputDir << std::endl;
    std::cout << "Output file name: " << outputFileName << std::endl;

    // Create an output ROOT file
    TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    if (!outputFile->IsOpen()) {
        std::cerr << "Error: Cannot open output file " << outputFileName << std::endl;
        return;
    }

    // Loop over files in the directory
    TSystemDirectory dir("dir", inputDir.c_str());
    TList *files = dir.GetListOfFiles();
    if (files) {
        TSystemFile *file;
        TString fileName;
        TIter next(files);
        while ((file=(TSystemFile*)next())) {
            fileName = file->GetName();
            if (!file->IsDirectory() && fileName.EndsWith(".root") && fileName != outputFileName) {
                std::cout << "Processing file: " << fileName << std::endl;

                // Open input ROOT file
                TFile* inputFile = TFile::Open((inputDir + "/" + fileName).Data());
                if (!inputFile || inputFile->IsZombie()) {
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
        }
    }

    // Write and close output ROOT file
    outputFile->Write();
    outputFile->Close();
    delete outputFile;
}

int main() {
    std::string inputDirectory = "/unix/muons/mu3e/Alex/Run/plot_files/parallelrun/signal_weightcut_nvert";
    std::string outputFileName = "/unix/muons/mu3e/Alex/Run/plot_files/parallelrun/merged_signal_histograms_nvert.root";
    MergeHistograms(inputDirectory, outputFileName);
    return 0;
}