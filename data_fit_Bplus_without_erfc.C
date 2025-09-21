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
    const int nbins_plot = 100; // Number of Bins for the Plot

    double min_signal = 5.178948768-0.048;
    double max_signal = 5.380091232;

    double mc_sigma1 = 0.03702;
    double mc_sigma2 = 0.01609;
    double mc_c1 = 0.3358;

    double xlow = 5.0;
    double xhigh = 5.8;

    double bin_width_plot = (xhigh - xlow) / nbins_plot;

    // Load ROOT file and TTree
    TFile *file = TFile::Open("data_unbinned_Bu_LekaiCut.root");
    // data_unbinned_Bu_final.root
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

    B_mass.setRange("fitRange", min_signal, xhigh);  // the range actually used in the fit


    // Signal model: Double Gaussian
    RooRealVar mean("mean", "Mean", 5.27764, 5.27, 5.29);

    // MC-derived widths (FIXED constants — put your values here)
    RooRealVar sigma1_mc("sigma1_mc", "MC Sigma1", mc_sigma1); 
    sigma1_mc.setConstant(kTRUE);

    RooRealVar sigma2_mc("sigma2_mc", "MC Sigma2", mc_sigma2); 
    sigma2_mc.setConstant(kTRUE);

    RooRealVar c1("c1", "Fraction of Gaussian1", mc_c1);
    c1.setConstant(kTRUE);

    // Common positive scale (fit parameter) and effective widths
    RooRealVar Cs("Cs", "Resolution scale", 1.0, 0.2, 3.0);
    RooProduct sigma1_eff("sigma1_eff", "sigma1_eff", RooArgList(sigma1_mc, Cs));
    RooProduct sigma2_eff("sigma2_eff", "sigma2_eff", RooArgList(sigma2_mc, Cs));


    RooGaussian gauss1("gauss1", "Gaussian 1", B_mass, mean, sigma1_eff);
    RooGaussian gauss2("gauss2", "Gaussian 2", B_mass, mean, sigma2_eff);
    RooAddPdf signal("signal", "Double Gaussian Model", RooArgList(gauss1, gauss2), RooArgList(c1));
    RooRealVar Nsig("Nsig", "Signal Yield", 3441.2, 0, 96300000);
    RooExtendPdf signal_ext("signal_ext", "Extended Signal", signal, Nsig);

    // Background model: Exponential + Gaussian (left-side component)
    RooRealVar lambda("lambda", "Lambda", -2.1699, -6.32, -0.01);
    RooExponential expo("expo", "Exponential Background", B_mass, lambda);
    RooRealVar Nbkg("Nbkg", "Exponential Background Yield", 1993.8, 0, 968000000);
    RooExtendPdf expo_ext("expo_ext", "Extended Exponential", expo, Nbkg);

    RooAddPdf background("background", "Total Background", RooArgList(expo_ext));

    // Full model = signal + background
    RooAddPdf model("model", "Signal + Background", RooArgList(signal_ext, background));

    // Fit the model to data (Extended Maximum Likelihood)
    RooFitResult* result = model.fitTo(dataset, Save(), Range(min_signal, xhigh));

    // Fraction of the signal pdf that lies in the fit range (normalized to full [xlow,xhigh])
    RooAbsReal* fracSigInFit = signal.createIntegral(B_mass, NormSet(B_mass), Range("fitRange"));
    // Convert fitted Nsig (in fit range) to an equivalent yield over the full plotting window [xlow,xhigh]
    double Nsig_fullRange     = Nsig.getVal()   / fracSigInFit->getVal();
    double Nsig_fullRange_err = Nsig.getError() / fracSigInFit->getVal();
    delete fracSigInFit;


    /*
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
    TFile *file_mc = TFile::Open("/lstore/cms/lekai/Bmeson/MC/ppRef/Bu_phat5_Bfinder.root");
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
    TString cut_mc = Form("Bchi2cl>0.003 && Balpha<0.170266 && Btrk1dR<1.70505 && (%s) && (%s) && (%s) && (%s)",
                        isMCsignal.Data(),
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
    */

    TCanvas* c = new TCanvas("c", "Bmass Fit with Pulls", 800, 800);
    c->Divide(1, 2);

    // Top pad (fit)
    TPad* p1 = (TPad*)c->cd(1);
    p1->SetPad(0.0, 0.15, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->Draw();
    p1->cd();


    // Draw the fit on the top frame (match styles to X3872 code)
    RooPlot* frame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));
    dataset.plotOn(frame, Binning(nbins_plot), MarkerStyle(20), MarkerSize(1.2), Name("data"), DataError(RooAbsData::Poisson));

    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global"));
    model.plotOn(frame, Components(expo_ext), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("background"));
    model.plotOn(frame, Components(signal), LineColor(kGreen+2), LineStyle(kDashed), LineWidth(2), Name("signal"));

    frame->SetTitle("");
    frame->GetYaxis()->SetTitle(Form("Events / ( %.4f )", bin_width_plot));
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0); // hide x labels on top pad
    frame->Draw();

    /*
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
    */



    // Calculate chi2/ndf for the fit
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare("global", "data", nParams);

    // Create a legend (top-right) for the plot
    p1->cd();
    TLegend* legend = new TLegend(0.48, 0.60, 0.88, 0.88);
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
    TPaveText* pave = new TPaveText(0.56, 0.16, 0.88, 0.60, "NDC");
    pave->SetTextAlign(12);
    pave->SetTextFont(42);
    pave->SetTextSize(0.025);
    pave->SetFillColor(0);
    pave->SetBorderSize(1);

    // Signal: Double Gaussian
    pave->AddText(Form("Mean = %.5f #pm %.5f", mean.getVal(), mean.getError()));
    pave->AddText(Form("#sigma_{1} (fixed) = %.5f", sigma1_mc.getVal()));
    pave->AddText(Form("#sigma_{2} (fixed) = %.5f", sigma2_mc.getVal()));
    pave->AddText(Form("c_{1} (fixed) = %.4f", c1.getVal()));
    pave->AddText(Form("C_{s} = %.5f #pm %.5f", Cs.getVal(), Cs.getError()));
    pave->AddText(Form("N_{sig} = %.1f #pm %.1f", Nsig.getVal(), Nsig.getError()));
    pave->AddText(Form("N_{sig}^{[%.1f, %.1f]} (rescaled) = %.1f #pm %.1f", xlow, xhigh, Nsig_fullRange, Nsig_fullRange_err));
    // Exponential background
    pave->AddText(Form("#lambda = %.4f #pm %.4f", lambda.getVal(), lambda.getError()));
    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    // Chi2
    pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2));

    pave->Draw();

    /*
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
    */


    // ---------- Bottom pad (pulls) ----------
    TPad* p2 = (TPad*)c->cd(2);
    p2->SetPad(0.0, 0.0, 1.0, 0.15);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.25);
    p2->Draw();
    p2->cd();

    RooPlot* pullFrame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));
    RooHist* pullHist = frame->pullHist("data", "global");

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

    /*
    // Console output summary
    std::cout << "Double Gaussian + Exponential fit complete. Output saved to " << name_file << std::endl;
    std::cout << std::fixed << std::setprecision(2) << std::endl;
    std::cout << "R3 (bkg in signal region) = " << bkg_in_signal << " events" << std::endl;
    std::cout << "R1+R2 (bkg in sidebands) = " << bkg_out_signal << " events" << std::endl;
    std::cout << "f_b = " << f_b << std::endl << std::endl;
    std::cout << "S_data (signal in region) = " << sig_yield_in_region << " events" << std::endl;
    std::cout << "S_MC = " << mc_yield_in_signal << " events" << std::endl;
    std::cout << "f_s = " << f_s << std::endl << std::endl;
    */

    // Clean up
    /*
    delete line_low;
    delete line_high;
    */

    delete zeroLine;
    delete c;
}



void data_fit_Bplus_without_erfc() {
    total_data_fit_Bu();
}