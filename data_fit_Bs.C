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
#include <RooExtendPdf.h>


// Aux: read y-value from the drawn curve at a given mass (for vertical line heights)
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



// B0 Particle
void total_data_fit_Bd() {
    const int nbins_plot = 100; // Number of bins for the plotted data

    double mc_sigma1 = 0.02539;
    double mc_sigma2 = 0.01039;
    double mc_c1 = 0.3723;

    double min_signal = 5.314800259;
    double max_signal = 5.420059741;

    double xlow = 5.2;
    double xhigh = 5.6;
    const double bin_width_plot = (xhigh - xlow) / nbins_plot;  // Used in y-axis label

    
    // Load real data tree
    TFile* file = TFile::Open("data_unbinned_Phi_FirstCut.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open real data file." << std::endl;
        return;
    }
    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cerr << "Error: TTree not found in file." << std::endl;
        return;
    }

    // RooFit: variable and data
    RooRealVar B_mass("B_mass", "B_mass", xlow, xhigh);
    B_mass.setRange("gaussRange", min_signal, max_signal);

    RooBinning mainBins(nbins_plot, xlow, xhigh);
    B_mass.setBinning(mainBins, "mainBins");

    RooDataSet dataset("dataset", "Unbinned dataset from TTree", tree, RooArgSet(B_mass));

    // Signal: Double Gaussian
    RooRealVar mean("mean", "Mean", 5.36743, 5.3, 5.4); 

    // MC-derived widths (FIXED constants — put your values here)
    RooRealVar sigma1("sigma1", "MC Sigma1", mc_sigma1); 
    sigma1.setConstant(kTRUE);

    RooRealVar sigma2("sigma2", "MC Sigma2", mc_sigma2); 
    sigma2.setConstant(kTRUE);

    RooRealVar c1("c1", "Fraction of Gaussian1", mc_c1);
    c1.setConstant(kTRUE);

    RooGaussian gaus1("gaus1", "Gaussian 1", B_mass, mean, sigma1);
    RooGaussian gaus2("gaus2", "Gaussian 2", B_mass, mean, sigma2);

    // Signal shape becomes a sum of the two Gaussians
    RooAddPdf signal("signal", "Double Gaussian", RooArgList(gaus1, gaus2), RooArgList(c1));


    RooRealVar Nsig("Nsig", "Signal Yield", 261, 0, 300000);
    RooExtendPdf signal_ext("signal_ext", "Extended Signal", signal, Nsig);

    // Background: Exponential model
    RooRealVar lambda("lambda", "Lambda", -1.72, -6.0, -0.1);
    RooExponential expo("expo", "Background", B_mass, lambda);

    RooRealVar Nbkg("Nbkg", "Background Yield", 1614, 0, 1700000);
    RooExtendPdf expo_ext("expo_ext", "Extended Exponential Background", expo, Nbkg);

    // Combined model
    RooAddPdf model("model", "Signal + Background", RooArgList(signal_ext, expo_ext));

    // Fit (Extended Maximum Likelihood Method)
    RooFitResult* result = model.fitTo(dataset, Save());

    
    // Compute background-only integrals in signal and sideband regions
    B_mass.setRange("signalRegion", min_signal, max_signal);
    B_mass.setRange("lowSideband", xlow, min_signal);
    B_mass.setRange("highSideband", max_signal, xhigh);

    double frac_bkg_signal = expo.createIntegral(B_mass, NormSet(B_mass), Range("signalRegion"))->getVal(); // Background in Signal Region
    double frac_bkg_low    = expo.createIntegral(B_mass, NormSet(B_mass), Range("lowSideband"))->getVal(); // Background in Left Noise Region
    double frac_bkg_high   = expo.createIntegral(B_mass, NormSet(B_mass), Range("highSideband"))->getVal(); // Background in Right Noise Region

    double total_bkg_yield = Nbkg.getVal(); // Total background

    double bkg_in_signal = total_bkg_yield * frac_bkg_signal; // Ammount of Noise in Signal Region 
    double bkg_out_signal = total_bkg_yield * (frac_bkg_low + frac_bkg_high); // Ammount of Noise in Sidebands
    double f_b = bkg_in_signal / bkg_out_signal; // Calculating F_b

    double frac_sig_in_signal = signal.createIntegral(B_mass, NormSet(B_mass), Range("signalRegion"))->getVal(); // Signal in Signal Region
    double sig_yield_in_region = Nsig.getVal() * frac_sig_in_signal; // Ammount of Signal in Signal Region



    // Opening and Checking MC File
    TFile *file_mc = TFile::Open("/lstore/cms/lekai/Bmeson/MC/ppRef/Bs_phat5_Bfinder.root");
    if (!file_mc || file_mc->IsZombie()) {
        std::cerr << "Error: Could not open MC file." << std::endl;
        return;
    }

    TTree *treemc = nullptr;
    file_mc->GetObject("Bfinder/ntphi", treemc);
    if (!treemc) {
        std::cerr << "Error: MC TTree not found!" << std::endl;
        return;
    }

    // Apply same cuts
    TString cut_mc = Form("Bchi2cl>0.003 && Bnorm_svpvDistance>2 && Bnorm_svpvDistance_2D>3.2232 && (%s) && (%s) && (%s) && (%s)",
                        isMCsignal.Data(),
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



    // ---------- Canvas with two pads ----------
    TCanvas* c = new TCanvas("c", "Bs Fit with Pulls", 800, 800);
    c->Divide(1, 2);

    // ---------- Top pad (fit) ----------
    TPad* p1 = (TPad*)c->cd(1);
    p1->SetPad(0.0, 0.15, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->Draw();
    p1->cd();

    RooPlot* frame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));

    // Plot data + model with the same naming/styles used for pulls
    dataset.plotOn(frame, Binning(B_mass.getBinning("mainBins")), MarkerStyle(20), MarkerSize(1.2), Name("data"), DataError(RooAbsData::Poisson));
    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global"));  // Total model
    model.plotOn(frame, Components(expo_ext), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("background"));
    model.plotOn(frame, Components(signal), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("signal"));

    frame->SetTitle("");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0);  // hide x labels on top pad
    frame->GetYaxis()->SetTitle(Form("Events / ( %.4f )", bin_width_plot));
    frame->Draw();

    // Vertical dashed lines at signal-region edges (heights taken from drawn curve)
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

    // Chi2 after plotting on frame
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare("global", "data", nParams);


    // ---------- Legend (same place), on TOP pad ----------
    p1->cd();
    TLegend* legend = new TLegend(0.58, 0.66, 0.88, 0.88);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetFillStyle(0);
    legend->AddEntry(frame->findObject("data"), "Data (B_{s} )", "lep");
    legend->AddEntry(frame->findObject("background"), "Background Fit (Exponential)", "l");
    legend->AddEntry(frame->findObject("signal"), "Signal Fit (Double Gaussian)", "l");
    legend->AddEntry(frame->findObject("global"), "Signal + Background Fit", "l");
    legend->Draw();

    // ---------- TPaveText (same place), on TOP pad ----------
    p1->cd();
    TPaveText* pave = new TPaveText(0.63, 0.38, 0.88, 0.66, "NDC");
    pave->SetTextAlign(12);
    pave->SetTextFont(42);
    pave->SetTextSize(0.025);
    pave->SetFillColor(0);
    pave->SetBorderSize(1);
    pave->AddText(Form("Mean = %.5f #pm %.5f", mean.getVal(), mean.getError()));

    pave->AddText(Form("#sigma_{1} (fixed) = %.5f", sigma1.getVal()));
    pave->AddText(Form("#sigma_{2} (fixed) = %.5f", sigma2.getVal()));
    pave->AddText(Form("c1 (fixed) = %.4f", c1.getVal()));

    pave->AddText(Form("N_{sig} = %.1f #pm %.1f", Nsig.getVal(), Nsig.getError()));
    pave->AddText(Form("#lambda = %.5f #pm %.5f", lambda.getVal(), lambda.getError()));
    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    pave->AddText(Form("#chi^{2}/ndf = %.5f", chi2));
    pave->Draw();

    // ---------- f_b / f_s box (same place), on TOP pad ----------
    p1->cd();
    TPaveText* pave_fb_fs = new TPaveText(0.46, 0.77, 0.58, 0.88, "NDC");
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
    p2->SetPad(0.0, 0.0, 1.0, 0.15);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.25);
    p2->Draw();
    p2->cd();

    RooPlot* pullFrame = B_mass.frame();
    RooHist* pullHist = frame->pullHist("data", "global");   // names must match
    pullHist->SetMarkerSize(0.6);
    pullFrame->addPlotable(pullHist, "XP");

    pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetNdivisions(505);
    pullFrame->GetYaxis()->SetTitleSize(0.10);
    pullFrame->GetYaxis()->SetTitleOffset(0.40);
    pullFrame->GetYaxis()->SetLabelSize(0.08);
    pullFrame->GetXaxis()->SetTitle("m_{J/#Psi #Phi} [GeV/c^{2}]");
    pullFrame->GetXaxis()->SetTitleSize(0.10);
    pullFrame->GetXaxis()->SetTitleOffset(1.0);
    pullFrame->GetXaxis()->SetLabelSize(0.08);
    pullFrame->SetMinimum(-3.5);
    pullFrame->SetMaximum(3.5);
    pullFrame->Draw("AP");

    TLine* zeroLine = new TLine(xlow, 0, xhigh, 0);
    zeroLine->SetLineColor(kBlue);
    zeroLine->SetLineStyle(1);
    zeroLine->SetLineWidth(1);
    zeroLine->Draw("same");


    TString name = "Bs_First_Fit.pdf";
    c->SaveAs(name);

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

    std::cout << "Output saved to " << name << std::endl;

    delete c;
    delete tree;
    delete line_low;
    delete line_high;
}



void data_fit_Bs() {
    total_data_fit_Bd();
}