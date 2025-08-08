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
#include <RooCurve.h>




// Aux function to get the y-value of the drawn model at a given mass (for vertical lines)
double getYatMass(RooPlot* frame, double mass) {
    for (int i = 0; i < frame->numItems(); ++i) {
        RooCurve* curve = dynamic_cast<RooCurve*>(frame->getObject(i));
        if (!curve) continue;
        int n = curve->GetN();
        double* x = curve->GetX();
        double* y = curve->GetY();
        for (int j = 0; j < n - 1; ++j) {
            if (x[j] <= mass && mass <= x[j+1]) {
                double slope = (y[j+1] - y[j]) / (x[j+1] - x[j]);
                return y[j] + slope * (mass - x[j]);
            }
        }
    }
    return 0.0;
}

// B+ Particle
void total_data_fit_Bu() {
    double min_signal = 5.1732;
    double max_signal = 5.3876;

    double xlow = 5.0;
    double xhigh = 6.0;
    double bin_width = 0.01;
    int nbins = int((xhigh - xlow) / bin_width);

    // Load ROOT file and TTree
    TFile *file = TFile::Open("data_unbinned_Bu_final.root");
    // dataCutted_Bu.root
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open real data file." << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cerr << "Error: Tree 'tree' not found in file." << std::endl;
        return;
    }

    // Define the mass variable and dataset
    RooRealVar B_mass("B_mass", "B_mass", xlow, xhigh);

    // Create unbinned RooDataSet from TTree
    B_mass.setRange("gaussRange", min_signal, max_signal);
    RooDataSet dataset("dataset", "Unbinned dataset from TTree", tree, RooArgSet(B_mass));

    // Signal model: Double Gaussian
    RooRealVar mean("mean", "Mean", 5.28, 5.277, 5.283);
    RooRealVar sigma1("sigma1", "Sigma1", 0.019, 0.005, 0.04);
    RooRealVar sigma2("sigma2", "Sigma2", 0.04, 0.02, 0.06);
    RooRealVar c1("c1", "Fraction of Gaussian1", 0.61, 0.55, 0.65);
    RooGaussian gauss1("gauss1", "Narrow Gaussian", B_mass, mean, sigma1);
    RooGaussian gauss2("gauss2", "Wide Gaussian", B_mass, mean, sigma2);
    RooAddPdf signal("signal", "Double Gaussian Model", RooArgList(gauss1, gauss2), RooArgList(c1));
    RooRealVar Nsig("Nsig", "Signal Yield", 1178, 1000, 1300);
    RooExtendPdf signal_ext("signal_ext", "Extended Signal", signal, Nsig);

    // Background model: Exponential + Gaussian (left-side component)
    RooRealVar lambda("lambda", "Lambda", -2.48, -3.00, -2.00);
    RooExponential expo("expo", "Exponential Background", B_mass, lambda);
    RooRealVar Nbkg("Nbkg", "Exponential Background Yield", 530, 300, 800);
    RooExtendPdf expo_ext("expo_ext", "Extended Exponential", expo, Nbkg);

    RooAddPdf background("background", "Total Background", RooArgList(expo_ext));

    // Full model = signal + background
    RooAddPdf model("model", "Signal + Background", RooArgList(signal_ext, background));

    // Fit the model to data (Extended Maximum Likelihood)
    RooFitResult* result = model.fitTo(dataset, Save(), Range(min_signal, xhigh));

    // Define regions for integrals
    B_mass.setRange("signalRegion", min_signal, max_signal);
    B_mass.setRange("lowSideband", xlow, min_signal);
    B_mass.setRange("highSideband", max_signal, xhigh);

    // Compute background-only integrals in signal and sideband regions
    double frac_bkg_signal = expo_ext.createIntegral(B_mass, NormSet(B_mass), Range("signalRegion"))->getVal();
    double frac_bkg_high   = expo_ext.createIntegral(B_mass, NormSet(B_mass), Range("highSideband"))->getVal();

    double total_bkg_yield = Nbkg.getVal();
    double bkg_in_signal   = total_bkg_yield * frac_bkg_signal;   // Background in signal region (R3)
    double bkg_out_signal  = total_bkg_yield * frac_bkg_high;     // Background in upper sideband (R2 only!)
    double f_b = bkg_in_signal / bkg_out_signal;                 // f_b calculation

    // Compute signal yield in signal region
    double frac_sig_in_signal = signal.createIntegral(B_mass, NormSet(B_mass), Range("signalRegion"))->getVal();
    double sig_yield_in_region = Nsig.getVal() * frac_sig_in_signal;  // Signal in signal region (S_data)

    // Open and process the MC file for signal region yield
    TFile *file_mc = TFile::Open("/lstore/cms/hlegoinha/Bmeson/MC_DATA/MC_ppRef_Bmeson/Bu_phat5_Bfinder.root");
    if (!file_mc || file_mc->IsZombie()) {
        std::cerr << "Error: Could not open MC file." << std::endl;
        return;
    }
    TTree *treemc = nullptr;
    file_mc->GetObject("Bfinder/ntKp", treemc);
    if (!treemc) {
        std::cerr << "Error: MC TTree not found!" << std::endl;
        file_mc->Close();
        return;
    }

    // Apply the same cuts as data
    TString cut_mc = Form("Bchi2cl>0.05 && Balpha<0.015 && Btrk1dR<1.66266 && (%s) && (%s) && (%s)",
                           ACCcuts_ppRef_Bu.Data(),
                           SELcuts_ppRef_Bu.Data(),
                           TRGmatching.Data());
    int nbins_mc = int((max_signal - min_signal) / 0.01);
    TH1F *hist_mc = new TH1F("hist_mc", "MC Bmass in Signal Region; Bmass [GeV/c^{2}]; Entries", nbins_mc, min_signal, max_signal);
    treemc->Draw("Bmass >> hist_mc", cut_mc + Form(" && Bmass > %.4f && Bmass < %.4f", min_signal, max_signal), "goff");
    double mc_yield_in_signal = hist_mc->Integral();  // S_MC
    delete hist_mc;
    file_mc->Close();

    double f_s = sig_yield_in_region / mc_yield_in_signal;  // f_s calculation


    TCanvas* c = new TCanvas("c", "Bmass Fit with Pulls", 800, 800);
    c->Divide(1, 2);

    // Top pad (fit)
    TPad* p1 = (TPad*)c->cd(1);
    p1->SetPad(0.0, 0.25, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->Draw();
    p1->cd();

    // Draw the fit on the top frame (match styles to X3872 code)
    RooPlot* frame = B_mass.frame();
    dataset.plotOn(frame, MarkerStyle(20), MarkerSize(1.2), Name("data"), DataError(RooAbsData::Poisson));
    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global"));
    model.plotOn(frame, Components(expo_ext), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("background"));
    model.plotOn(frame, Components(signal), LineColor(kGreen+2), LineStyle(kDashed), LineWidth(2), Name("signal"));

    frame->SetTitle("");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0); // hide x labels on top pad
    frame->Draw();

    // Vertical dashed lines at signal-region edges (line heights taken from the drawn curve)
    double y_low  = getYatMass(frame, min_signal);
    double y_high = getYatMass(frame, max_signal);

    TLine* line_low  = new TLine(min_signal, 0, min_signal, y_low);
    TLine* line_high = new TLine(max_signal, 0, max_signal, y_high);
    for (TLine* l : {line_low, line_high}) {
        l->SetLineColor(kBlack);
        l->SetLineStyle(2);
        l->SetLineWidth(2);
        l->Draw("same");
    }



    // Calculate chi2/ndf for the fit
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare(nParams);

    // Create a legend (top-right) for the plot
    p1->cd();
    TLegend* legend = new TLegend(0.56, 0.66, 0.88, 0.88);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetFillStyle(0);
    legend->AddEntry(frame->findObject("data"), "Data (B^{+}) Unbinned", "lep");
    legend->AddEntry(frame->findObject("background"), "Background Fit (Exponential)", "l");
    legend->AddEntry(frame->findObject("signal"), "Signal Fit (Double Gaussian)", "l");
    legend->AddEntry(frame->findObject("global"), "Total Fit (Signal + Background)", "l");
    legend->Draw();

    // TPaveText for fit parameters (bottom-right)
    p1->cd();
    TPaveText* pave = new TPaveText(0.64, 0.30, 0.88, 0.66, "NDC");
    pave->SetTextAlign(12);
    pave->SetTextFont(42);
    pave->SetTextSize(0.025);
    pave->SetFillColor(0);
    pave->SetBorderSize(1);

    // Signal: Double Gaussian
    pave->AddText(Form("Mean = %.5f #pm %.5f", mean.getVal(), mean.getError()));
    pave->AddText(Form("#sigma_{1} = %.5f #pm %.5f", sigma1.getVal(), sigma1.getError()));
    pave->AddText(Form("#sigma_{2} = %.5f #pm %.5f", sigma2.getVal(), sigma2.getError()));
    pave->AddText(Form("c_{1} = %.4f #pm %.4f", c1.getVal(), c1.getError()));
    pave->AddText(Form("N_{sig} = %.1f #pm %.1f", Nsig.getVal(), Nsig.getError()));
    // Exponential background
    pave->AddText(Form("#lambda = %.4f #pm %.4f", lambda.getVal(), lambda.getError()));
    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    // Chi2
    pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2));

    pave->Draw();

    // TPaveText for f_b and f_s (top-left or another suitable position)
    p1->cd();
    TPaveText* pave_fb_fs = new TPaveText(0.44, 0.77, 0.56, 0.88, "NDC");
    pave_fb_fs->SetTextAlign(12);
    pave_fb_fs->SetTextFont(42);
    pave_fb_fs->SetTextSize(0.025);
    pave_fb_fs->SetFillColor(0);
    pave_fb_fs->SetBorderSize(1);
    pave_fb_fs->AddText(Form("f_{b} = %.3f", f_b));
    pave_fb_fs->AddText(Form("f_{s} = %.3f", f_s));
    pave_fb_fs->Draw();


    // ---------- Bottom pad (pulls) ----------
    TPad* p2 = (TPad*)c->cd(2);
    p2->SetPad(0.0, 0.0, 1.0, 0.25);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.25);
    p2->Draw();
    p2->cd();

    RooPlot* pullFrame = B_mass.frame();
    RooHist* pullHist = frame->pullHist("data", "global"); // names must match plotOn Name(...)
    pullHist->SetMarkerSize(0.6);
    pullFrame->addPlotable(pullHist, "XP");

    pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetNdivisions(505);
    pullFrame->GetYaxis()->SetTitleSize(0.10);
    pullFrame->GetYaxis()->SetTitleOffset(0.40);
    pullFrame->GetYaxis()->SetLabelSize(0.08);
    pullFrame->GetXaxis()->SetTitle("m_{J/#Psi K^{+}} [GeV/c^{2}]");
    pullFrame->GetXaxis()->SetTitleSize(0.10);
    pullFrame->GetXaxis()->SetTitleOffset(1.0);
    pullFrame->GetXaxis()->SetLabelSize(0.08);
    pullFrame->SetMinimum(-3.5);
    pullFrame->SetMaximum(3.5);
    pullFrame->Draw("AP");

    // Zero line
    TLine* zeroLine = new TLine(xlow, 0, xhigh, 0);
    zeroLine->SetLineColor(kBlue);
    zeroLine->SetLineStyle(1);
    zeroLine->SetLineWidth(1);
    zeroLine->Draw("same");


    // Save the canvas to a file
    TString name_file = "Bu_Total_Fit_with_Pulls.pdf";
    c->SaveAs(name_file);

    // Console output summary
    std::cout << "Double Gaussian + Exponential fit complete. Output saved to " << name_file << std::endl;
    std::cout << std::fixed << std::setprecision(2) << std::endl;
    std::cout << "R3 (bkg in signal region) = " << bkg_in_signal << " events" << std::endl;
    std::cout << "R1+R2 (bkg in sidebands) = " << bkg_out_signal << " events" << std::endl;
    std::cout << "f_b = " << f_b << std::endl << std::endl;
    std::cout << "S_data (signal in region) = " << sig_yield_in_region << " events" << std::endl;
    std::cout << "S_MC = " << mc_yield_in_signal << " events" << std::endl;
    std::cout << "f_s = " << f_s << std::endl << std::endl;

    // Clean up
    delete line_low;
    delete line_high;
    delete zeroLine;
    delete c;
}



void data_fit_Bplus_without_erfc() {
    total_data_fit_Bu();
}