using namespace RooFit;

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>
#include "ACCSEL.h"

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooProduct.h>
#include <RooFormulaVar.h>
#include <RooExponential.h>
#include <RooGenericPdf.h>  
#include <RooProdPdf.h>

#include <RooFitResult.h>
#include <RooFit.h>
#include <RooCmdArg.h>




// B0 Particle
void total_data_fit_Bd() {
    double min_signal = 5.1721;
    double max_signal = 5.3868;

    double xlow = 5.0;
    double xhigh = 6.0;
    int nbins = 150;
    double bin_width = (xhigh - xlow)/nbins;

    
    // Load real data histogram
    TFile* fData = TFile::Open("DATA_ppRef_Bmass_Bd_150bin_2.root");
    if (!fData || fData->IsZombie()) {
        std::cerr << "Error: Could not open real data file." << std::endl;
        return;
    }
    TH1F* hist = (TH1F*)fData->Get("DataBkg");
    if (!hist) {
        std::cerr << "Error: Histogram 'DataBkg' not found in file." << std::endl;
        return;
    }

    if (hist->GetEntries() == 0) {
        std::cerr << "Error: Histogram 'DataBkg' is empty." << std::endl;
        return;
    }

    // RooFit: variable and data
    RooRealVar Bmass("Bmass", "Bmass", xlow, xhigh);
    Bmass.setRange("gaussRange", min_signal, max_signal);
    RooDataHist data("data", "dataset", RooArgList(Bmass), hist);

    // Signal: Double Gaussian
    RooRealVar mean("mean", "Mean", 5.29, 5.275, 5.3);
    RooRealVar sigma("sigma", "Sigma1", 0.009, 0.001, 0.08);

    // Gaussians
    RooGaussian signal("signal", "Single Gaussian", Bmass, mean, sigma);

    RooRealVar Nsig("Nsig", "Signal Yield", 261, 0, 300000);
    RooExtendPdf signal_ext("signal_ext", "Extended Signal", signal, Nsig);

    // Background: Exponential model
    RooRealVar lambda("lambda", "Lambda", -2.72, -6.0, -0.1);
    RooExponential expo("expo", "Background", Bmass, lambda);

    //RooRealVar c0("c0", "c0", 1.0, -1e5, 1e5);
    //RooRealVar c1("c1", "c1", 0.0, -1e5, 1e5);
    //RooRealVar c2("c2", "c2", 0.0, -1e5, 1e5);
    //RooRealVar c3("c3", "c3", 0.0, -1e5, 1e5);
    //RooRealVar c4("c4", "c4", 0.0, -1e5, 1e5);
    //RooPolynomial poly("poly", "4th-degree polynomial background", Bmass, RooArgList(c0, c1, c2, c3, c4));

    RooRealVar Nbkg("Nbkg", "Background Yield", 1614, 0, 1700000);
    RooExtendPdf expo_ext("expo_ext", "Extended Exponential Background", expo, Nbkg);

    // Combined model
    RooAddPdf model("model", "Signal + Background", RooArgList(signal_ext, expo_ext));

    // Fit (Extended Maximum Likelihood Method)
    RooFitResult* result = model.fitTo(data, Save());

    // Compute background-only integrals in signal and sideband regions
    Bmass.setRange("signalRegion", min_signal, max_signal);
    Bmass.setRange("lowSideband", xlow, min_signal);
    Bmass.setRange("highSideband", max_signal, xhigh);

    double frac_bkg_signal = expo.createIntegral(Bmass, NormSet(Bmass), Range("signalRegion"))->getVal(); // Background in Signal Region
    double frac_bkg_low    = expo.createIntegral(Bmass, NormSet(Bmass), Range("lowSideband"))->getVal(); // Background in Left Noise Region
    double frac_bkg_high   = expo.createIntegral(Bmass, NormSet(Bmass), Range("highSideband"))->getVal(); // Background in Right Noise Region

    double total_bkg_yield = Nbkg.getVal(); // Total background

    double bkg_in_signal = total_bkg_yield * frac_bkg_signal; // Ammount of Noise in Signal Region 
    double bkg_out_signal = total_bkg_yield * (frac_bkg_low + frac_bkg_high); // Ammount of Noise in Sidebands
    double f_b = bkg_in_signal / bkg_out_signal; // Calculating F_b

    double frac_sig_in_signal = signal.createIntegral(Bmass, NormSet(Bmass), Range("signalRegion"))->getVal(); // Signal in Signal Region
    double sig_yield_in_region = Nsig.getVal() * frac_sig_in_signal; // Ammount of Signal in Signal Region



    // Opening and Checking MC File
    TFile *file_mc = TFile::Open("/lstore/cms/henrique/Bmeson/MC_DATA/MC_ppRef_Bmeson/Bd_phat5_Bfinder.root");
    if (!file_mc || file_mc->IsZombie()) {
        std::cerr << "Error: Could not open MC file." << std::endl;
        return;
    }

    TTree *treemc = nullptr;
    file_mc->GetObject("Bfinder/ntKstar", treemc);
    if (!treemc) {
        std::cerr << "Error: MC TTree not found!" << std::endl;
        return;
    }

    // Apply same cuts
    TString cut_mc = Form("Balpha<0.184842 && Bnorm_svpvDistance_2D>3.9916 && Bchi2cl>0.05 && (%s) && (%s) && (%s)",
                        ACCcuts_ppRef.Data(),
                        SELcuts_ppRef.Data(),
                        TRGmatching.Data());

    int nbins_mc = 150;

    TH1F *hist_mc = new TH1F("hist_mc", "MC Bmass in Signal Region;Bmass [GeV/c^{2}];Entries", nbins_mc, min_signal, max_signal);

    treemc->Draw("Bmass >> hist_mc", cut_mc + Form(" && Bmass > %.4f && Bmass < %.4f", min_signal, max_signal), "goff");

    double mc_yield_in_signal = hist_mc->Integral();

    delete hist_mc;
    file_mc->Close();
    // --------------------------------------------------------------------------------------------------
    
    double f_s = sig_yield_in_region / mc_yield_in_signal; // Calculating F_s


    // Plotting
    TCanvas *c = new TCanvas("c", "Bmass Fit", 800, 600);
    RooPlot* frame = Bmass.frame();
    data.plotOn(frame, MarkerStyle(20), MarkerSize(1.0));

    // Plotting block with correct styles and legend names
    data.plotOn(frame, MarkerStyle(20), MarkerSize(1.2), Name("data"));
    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global"));  // Total model
    model.plotOn(frame, Components(expo_ext), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("background"));
    model.plotOn(frame, Components(signal), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("signal"));

    
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("m_{J/#Psi K^{*}} [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle(Form("Events / ( %.4f )", bin_width));
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->Draw();

    // Vertical lines at signal region edges
    double binWidth = 0.01;
    double total_yield = Nsig.getVal() + Nbkg.getVal();

    Bmass.setVal(min_signal);
    double y_low = model.getVal(RooArgSet(Bmass)) * total_yield * binWidth;

    Bmass.setVal(max_signal);
    double y_high = model.getVal(RooArgSet(Bmass)) * total_yield * binWidth;

    // Create vertical dashed black lines
    TLine* line_low = new TLine(min_signal, 0, min_signal, y_low);
    TLine* line_high = new TLine(max_signal, 0, max_signal, y_high);


    line_low->SetLineColor(kBlack);
    line_low->SetLineStyle(2);
    line_low->SetLineWidth(2);

    line_high->SetLineColor(kBlack);
    line_high->SetLineStyle(2);
    line_high->SetLineWidth(2);

    // Draw lines
    line_low->Draw("same");
    line_high->Draw("same");


    // Chi2 and fit parameters
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare(nParams);
    
    // Custom Legend (top-right, styled as image)
    TLegend* legend = new TLegend(0.58, 0.66, 0.88, 0.88);  // Adjust position as needed
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetFillStyle(0);  // transparent background

    legend->AddEntry(frame->findObject("data"), "Data (B^{0} )", "lep");
    legend->AddEntry(frame->findObject("background"), "Background Fit (Exponential)", "l");
    legend->AddEntry(frame->findObject("signal"), "Signal Fit (Single Gaussian)", "l");
    legend->AddEntry(frame->findObject("global"), "Signal + Background Fit", "l");

    legend->Draw();
        

    // TPaveText (bottom)
    TPaveText* pave = new TPaveText(0.63, 0.38, 0.88, 0.66, "NDC");
    pave->SetTextAlign(12);  // left-aligned
    pave->SetTextFont(42);   // Helvetica
    pave->SetTextSize(0.025);
    pave->SetFillColor(0);   // transparent
    pave->SetBorderSize(1);  // thin box

    // Add formatted lines
    pave->AddText(Form("Mean = %.5f #pm %.5f", mean.getVal(), mean.getError()));
    pave->AddText(Form("#sigma = %.5f #pm %.5f", sigma.getVal(), sigma.getError()));
    pave->AddText(Form("N_{sig} = %.1f #pm %.1f", Nsig.getVal(), Nsig.getError()));
    pave->AddText(Form("#lambda = %.5f #pm %.5f", lambda.getVal(), lambda.getError()));
    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2));

    pave->Draw();    // draw on top


    // xxAdditional Legend-Like Box for f_s and f_b
    TPaveText* pave_fb_fs = new TPaveText(0.46, 0.77, 0.58, 0.88, "NDC");  // Adjust position if needed
    pave_fb_fs->SetTextAlign(12);  // Left-aligned
    pave_fb_fs->SetTextFont(42);
    pave_fb_fs->SetTextSize(0.025);
    pave_fb_fs->SetFillColor(0);
    pave_fb_fs->SetBorderSize(1);  // thin frame

    pave_fb_fs->AddText(Form("f_{b} = %.3f", f_b));
    pave_fb_fs->AddText(Form("f_{s} = %.3f", f_s));

    pave_fb_fs->Draw();  // Draw it on top of the plot

    c->SaveAs("Bd_Total_Fit_Binned.pdf");

    std::cout << std::fixed << std::setprecision(2);
    std::cout << " " << std::endl;
    std::cout << "R3 = " << bkg_in_signal << " events" << std::endl;
    std::cout << "R1 + R2 = " << bkg_out_signal << " events" << std::endl;
    std::cout << "f_{b} = " << f_b << std::endl;

    std::cout << " " << std::endl;

    std::cout << "S_{data} " << sig_yield_in_region << " events" << std::endl;
    std::cout << "S_{MC} = " << mc_yield_in_signal << " events" << std::endl;
    std::cout << "f_{s} = " << f_s << std::endl;
    std::cout << " " << std::endl; 

    std::cout << "Gaussian + Exponential fit complete. Output saved to 'Bd_Total_Fit_Binned.pdf'" << std::endl;

    delete c;
    delete hist;
    delete line_low;
    delete line_high;
}



void data_fit_Bzero() {
    total_data_fit_Bd();
}