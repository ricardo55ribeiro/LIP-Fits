using namespace RooFit;

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>
#include <iomanip>
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
#include <RooPolynomial.h>
#include <RooChebychev.h>
#include <RooBernstein.h>
#include <RooCBShape.h>

#include <RooFitResult.h>
#include <RooFit.h>
#include <RooCmdArg.h>
#include <RooCurve.h>


/*
void Bplus_nominal() {
    const int nbins_plot = 100; // Number of bins for the plot

    double min_signal = 5.178948768;
    double max_signal = 5.380091232;

    double mc_mean = 5.27952;
    double mc_sigma1 = 0.03702;
    double mc_sigma2 = 0.01609;
    double mc_c1 = 0.3358;
    

    double xlow = 5.0;
    double xhigh = 5.8;

    double bin_width_plot = (xhigh - xlow) / nbins_plot;

    // Load ROOT file and TTree
    TFile *file = TFile::Open("data_unbinned_Bu_Opticut.root");


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

    // Define regions for integrals
    B_mass.setRange("signalRegion", min_signal, max_signal);
    B_mass.setRange("lowSideband", xlow, min_signal);
    B_mass.setRange("highSideband", max_signal, xhigh);


    // Signal model: Double Gaussian
    RooRealVar mean("mean", "Mean", 5.27764, 5.27, 5.29);


    // MC-derived widths (FIXED constants — put your values here)
    RooRealVar sigma1_mc("sigma1_mc", "MC Sigma1", mc_sigma1); 
    sigma1_mc.setConstant(kTRUE);

    RooRealVar sigma2_mc("sigma2_mc", "MC Sigma2", mc_sigma2); 
    sigma2_mc.setConstant(kTRUE);

    RooRealVar c1("c1", "Fraction of Gaussian1", mc_c1);
    c1.setConstant(kTRUE);


    // Common positive scale (fit parameter)
    RooRealVar Cs("Cs", "Resolution scale", 1.15888, 0.2, 3.0);

    // Effective widths = Cs * sigma_mc
    RooProduct sigma1_eff("sigma1_eff", "sigma1_eff", RooArgList(sigma1_mc, Cs));
    RooProduct sigma2_eff("sigma2_eff", "sigma2_eff", RooArgList(sigma2_mc, Cs));

    // Mixture fraction and Gaussians
    RooGaussian gauss1("gauss1", "Gaussian 1", B_mass, mean, sigma1_eff);
    RooGaussian gauss2("gauss2", "Gaussian 2", B_mass, mean, sigma2_eff);
    RooAddPdf signal("signal", "Double Gaussian Model", RooArgList(gauss1, gauss2), RooArgList(c1));
    RooRealVar Nsig("Nsig", "Signal Yield", 3441.2, 0, 96300000);
    

    // Background model
    RooRealVar lambda("lambda", "Lambda", -2.1699, -6.32, -0.01);
    RooExponential expo("expo", "Exponential Background", B_mass, lambda);
    RooRealVar Nbkg("Nbkg", "Exponential Background Yield", 1993.8, 0, 968000000);


    // ERFC background (for left sideband)
    RooRealVar csf("csf", "Shifting Constant", 5.14159, 5.08, 5.16);
    RooRealVar csc("csc", "Scaling Constant", 0.03528, 0.0008, 0.05);

    // integral form implemented via erf: 1 - erf(x)
    RooGenericPdf erfc_bkg("erfc_bkg", "1 - TMath::Erf((B_mass - csf)/csc)", RooArgList(B_mass, csf, csc));


    RooRealVar Nerfc("Nerfc", "ERFC Background Yield", 697.3, 0, 50000000);
    
    RooAddPdf model("model", "Signal + Background", RooArgList(signal, expo, erfc_bkg), RooArgList(Nsig, Nbkg, Nerfc));


    // Fit the model to data (Extended Maximum Likelihood)
    RooFitResult* result = model.fitTo(dataset, Save(), Range(xlow, xhigh), Extended(kTRUE));

    // ---------- Canvas with two pads ----------
    TCanvas* c = new TCanvas("c", "Bmass Fit with Pulls (ERFC model)", 800, 800);
    c->Divide(1, 2);

    // ---------- Top pad (fit) ----------
    TPad* p1 = (TPad*)c->cd(1);
    p1->SetPad(0.0, 0.15, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->Draw();
    p1->cd();

    RooPlot* frame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));

    // Plot data + model with the same naming/styles used for pulls
    dataset.plotOn(frame, Binning(nbins_plot), MarkerStyle(20), MarkerSize(1.2), Name("data"), DataError(RooAbsData::Poisson));
    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global")); // total
    model.plotOn(frame, Components(expo), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("background")); // exponential part
    model.plotOn(frame, Components(erfc_bkg), LineColor(kMagenta), LineStyle(kDotted), LineWidth(2), Name("erfc_bkg")); // erfc part
    model.plotOn(frame, Components(signal),   LineColor(kGreen+2), LineStyle(kDashed), LineWidth(2), Name("signal")); // signal

    frame->SetTitle("");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0);  // hide x labels on top pad
    frame->GetXaxis()->SetTitle("m_{J/#Psi K^{+}} [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle(Form("Events / ( %.4f )", bin_width_plot));
    frame->Draw();


    // Calculate chi2/ndf for the fit
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare("global", "data", nParams);

    // ---------- Legend (same place), on TOP pad ----------
    p1->cd();
    TLegend* legend = new TLegend(0.52, 0.68, 0.92, 0.96);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetFillStyle(1001);   // Solid fill
    legend->SetFillColor(kWhite); // White background (you can pick another color)
    legend->AddEntry(frame->findObject("data"), "Data (B^{+}) Unbinned", "lep");
    legend->AddEntry(frame->findObject("background"), "Background Fit (Exponential)", "l");
    legend->AddEntry(frame->findObject("erfc_bkg"), "#splitline{Partially Reconstructed}{Background Fit (Erfc)}","l");
    legend->AddEntry(frame->findObject("signal"), "Signal Fit (Double Gaussian)", "l");
    legend->AddEntry(frame->findObject("global"), "Total Fit (Signal + Background)", "l");
    legend->Draw();

    // ---------- TPaveText (same place), on TOP pad ----------
    p1->cd();
    TPaveText* pave = new TPaveText(0.60, 0.24, 0.92, 0.68, "NDC");
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
    // Exponential background
    pave->AddText(Form("#lambda = %.4f #pm %.4f", lambda.getVal(), lambda.getError()));
    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    // ERFC background
    pave->AddText(Form("c_{sf} = %.5f #pm %.5f", csf.getVal(), csf.getError()));
    pave->AddText(Form("c_{sc} = %.5f #pm %.5f", csc.getVal(), csc.getError()));
    pave->AddText(Form("N_{erfc} = %.1f #pm %.1f", Nerfc.getVal(), Nerfc.getError()));
    // Chi2
    pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2));
    pave->Draw();


    // ---------- Bottom pad (pulls) ----------
    TPad* p2 = (TPad*)c->cd(2);
    p2->SetPad(0.0, 0.0, 1.0, 0.15);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.25);
    p2->Draw();
    p2->cd();

    RooPlot* pullFrame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));
    RooHist* pullHist = frame->pullHist("data", "global");  // names must match Name("data") and Name("global")
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
    TString name_file = "Bu_Total_Fit_Nominal.pdf";
    c->SaveAs(name_file);


    delete zeroLine;
    delete c;
}
*/






/*
void Bplus_LinearBkg() {
    const int nbins_plot = 100; // Number of bins for the plot

    double min_signal = 5.178948768;
    double max_signal = 5.380091232;

    double mc_mean = 5.27952;
    double mc_sigma1 = 0.03702;
    double mc_sigma2 = 0.01609;
    double mc_c1 = 0.3358;
    

    double xlow = 5.0;
    double xhigh = 5.8;

    double bin_width_plot = (xhigh - xlow) / nbins_plot;

    // Load ROOT file and TTree
    TFile *file = TFile::Open("data_unbinned_Bu_Opticut.root");


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

    // Define regions for integrals
    B_mass.setRange("signalRegion", min_signal, max_signal);
    B_mass.setRange("lowSideband", xlow, min_signal);
    B_mass.setRange("highSideband", max_signal, xhigh);


    // Signal model: Double Gaussian
    RooRealVar mean("mean", "Mean", 5.27764, 5.27, 5.29);


    // MC-derived widths (FIXED constants — put your values here)
    RooRealVar sigma1_mc("sigma1_mc", "MC Sigma1", mc_sigma1); 
    sigma1_mc.setConstant(kTRUE);

    RooRealVar sigma2_mc("sigma2_mc", "MC Sigma2", mc_sigma2); 
    sigma2_mc.setConstant(kTRUE);

    RooRealVar c1("c1", "Fraction of Gaussian1", mc_c1);
    c1.setConstant(kTRUE);


    // Common positive scale (fit parameter)
    RooRealVar Cs("Cs", "Resolution scale", 1.15888, 0.2, 3.0);

    // Effective widths = Cs * sigma_mc
    RooProduct sigma1_eff("sigma1_eff", "sigma1_eff", RooArgList(sigma1_mc, Cs));
    RooProduct sigma2_eff("sigma2_eff", "sigma2_eff", RooArgList(sigma2_mc, Cs));

    // Mixture fraction and Gaussians
    RooGaussian gauss1("gauss1", "Gaussian 1", B_mass, mean, sigma1_eff);
    RooGaussian gauss2("gauss2", "Gaussian 2", B_mass, mean, sigma2_eff);
    RooAddPdf signal("signal", "Double Gaussian Model", RooArgList(gauss1, gauss2), RooArgList(c1));
    RooRealVar Nsig("Nsig", "Signal Yield", 3441.2, 0, 96300000);
    

    // Background model (linear instead of exponential)
    RooRealVar p_lin("p_lin", "Linear slope", 0.0, -5.0, 5.0);
    RooPolynomial linear("linear", "Linear Background", B_mass, RooArgList(p_lin), 1);
    RooRealVar Nbkg("Nbkg", "Linear Background Yield", 1993.8, 0, 968000000);


    // ERFC background (for left sideband)
    RooRealVar csf("csf", "Shifting Constant", 5.14159, 5.08, 5.16);
    RooRealVar csc("csc", "Scaling Constant", 0.03528, 0.0008, 0.05);

    // integral form implemented via erf: 1 - erf(x)
    RooGenericPdf erfc_bkg("erfc_bkg", "1 - TMath::Erf((B_mass - csf)/csc)", RooArgList(B_mass, csf, csc));


    RooRealVar Nerfc("Nerfc", "ERFC Background Yield", 697.3, 0, 50000000);
    
    RooAddPdf model("model", "Signal + Background", RooArgList(signal, linear, erfc_bkg), RooArgList(Nsig, Nbkg, Nerfc));


    // Fit the model to data (Extended Maximum Likelihood)
    RooFitResult* result = model.fitTo(dataset, Save(), Range(xlow, xhigh), Extended(kTRUE));

    // ---------- Canvas with two pads ----------
    TCanvas* c = new TCanvas("c", "Bmass Fit with Pulls (ERFC model)", 800, 800);
    c->Divide(1, 2);

    // ---------- Top pad (fit) ----------
    TPad* p1 = (TPad*)c->cd(1);
    p1->SetPad(0.0, 0.15, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->Draw();
    p1->cd();

    RooPlot* frame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));

    // Plot data + model with the same naming/styles used for pulls
    dataset.plotOn(frame, Binning(nbins_plot), MarkerStyle(20), MarkerSize(1.2), Name("data"), DataError(RooAbsData::Poisson));
    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global")); // total
    model.plotOn(frame, Components(linear), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("linear_bkg")); // linear part
    model.plotOn(frame, Components(erfc_bkg), LineColor(kMagenta), LineStyle(kDotted), LineWidth(2), Name("erfc_bkg")); // erfc part
    model.plotOn(frame, Components(signal),   LineColor(kGreen+2), LineStyle(kDashed), LineWidth(2), Name("signal")); // signal

    frame->SetTitle("");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0);  // hide x labels on top pad
    frame->GetXaxis()->SetTitle("m_{J/#Psi K^{+}} [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle(Form("Events / ( %.4f )", bin_width_plot));
    frame->Draw();


    // Calculate chi2/ndf for the fit
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare("global", "data", nParams);

    // ---------- Legend (same place), on TOP pad ----------
    p1->cd();
    TLegend* legend = new TLegend(0.52, 0.68, 0.92, 0.96);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetFillStyle(1001);   // Solid fill
    legend->SetFillColor(kWhite); // White background (you can pick another color)
    legend->AddEntry(frame->findObject("data"), "Data (B^{+}) Unbinned", "lep");
    legend->AddEntry(frame->findObject("linear_bkg"), "Background Fit (Linear)", "l");
    legend->AddEntry(frame->findObject("erfc_bkg"), "#splitline{Partially Reconstructed}{Background Fit (Erfc)}","l");
    legend->AddEntry(frame->findObject("signal"), "Signal Fit (Double Gaussian)", "l");
    legend->AddEntry(frame->findObject("global"), "Total Fit (Signal + Background)", "l");
    legend->Draw();

    // ---------- TPaveText (same place), on TOP pad ----------
    p1->cd();
    TPaveText* pave = new TPaveText(0.60, 0.24, 0.92, 0.68, "NDC");
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
    // Linear background
    pave->AddText(Form("p_{1} (slope) = %.4f #pm %.4f", p_lin.getVal(), p_lin.getError()));
    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    // ERFC background
    pave->AddText(Form("c_{sf} = %.5f #pm %.5f", csf.getVal(), csf.getError()));
    pave->AddText(Form("c_{sc} = %.5f #pm %.5f", csc.getVal(), csc.getError()));
    pave->AddText(Form("N_{erfc} = %.1f #pm %.1f", Nerfc.getVal(), Nerfc.getError()));
    // Chi2
    pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2));
    pave->Draw();


    // ---------- Bottom pad (pulls) ----------
    TPad* p2 = (TPad*)c->cd(2);
    p2->SetPad(0.0, 0.0, 1.0, 0.15);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.25);
    p2->Draw();
    p2->cd();

    RooPlot* pullFrame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));
    RooHist* pullHist = frame->pullHist("data", "global");  // names must match Name("data") and Name("global")
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
    TString name_file = "Bu_Total_Fit_LinearBkg.pdf";
    c->SaveAs(name_file);


    delete zeroLine;
    delete c;
}
*/




/*
void Bplus_2degreepoly() {
    const int nbins_plot = 100; // Number of bins for the plot

    double min_signal = 5.178948768;
    double max_signal = 5.380091232;

    double mc_mean = 5.27952;
    double mc_sigma1 = 0.03702;
    double mc_sigma2 = 0.01609;
    double mc_c1 = 0.3358;
    

    double xlow = 5.0;
    double xhigh = 5.8;

    double bin_width_plot = (xhigh - xlow) / nbins_plot;

    // Load ROOT file and TTree
    TFile *file = TFile::Open("data_unbinned_Bu_Opticut.root");


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

    // Define regions for integrals
    B_mass.setRange("signalRegion", min_signal, max_signal);
    B_mass.setRange("lowSideband", xlow, min_signal);
    B_mass.setRange("highSideband", max_signal, xhigh);


    // Signal model: Double Gaussian
    RooRealVar mean("mean", "Mean", 5.27764, 5.27, 5.29);


    // MC-derived widths (FIXED constants — put your values here)
    RooRealVar sigma1_mc("sigma1_mc", "MC Sigma1", mc_sigma1); 
    sigma1_mc.setConstant(kTRUE);

    RooRealVar sigma2_mc("sigma2_mc", "MC Sigma2", mc_sigma2); 
    sigma2_mc.setConstant(kTRUE);

    RooRealVar c1("c1", "Fraction of Gaussian1", mc_c1);
    c1.setConstant(kTRUE);


    // Common positive scale (fit parameter)
    RooRealVar Cs("Cs", "Resolution scale", 1.15888, 0.2, 3.0);

    // Effective widths = Cs * sigma_mc
    RooProduct sigma1_eff("sigma1_eff", "sigma1_eff", RooArgList(sigma1_mc, Cs));
    RooProduct sigma2_eff("sigma2_eff", "sigma2_eff", RooArgList(sigma2_mc, Cs));

    // Mixture fraction and Gaussians
    RooGaussian gauss1("gauss1", "Gaussian 1", B_mass, mean, sigma1_eff);
    RooGaussian gauss2("gauss2", "Gaussian 2", B_mass, mean, sigma2_eff);
    RooAddPdf signal("signal", "Double Gaussian Model", RooArgList(gauss1, gauss2), RooArgList(c1));
    RooRealVar Nsig("Nsig", "Signal Yield", 3441.2, 0, 96300000);
    
    // Background model (Bernstein order 2, strictly non-negative)
    RooRealVar c0("c0","Bernstein c0", 0.33, 0.0, 10.0);
    RooRealVar c1b("c1b","Bernstein c1", 0.33, 0.0, 10.0); // 'c1b' to avoid name clash with signal 'c1'
    RooRealVar c2("c2","Bernstein c2", 0.33, 0.0, 10.0);
    RooBernstein bern2("bern2","Bernstein Background (order 2)", B_mass, RooArgList(c0,c1b,c2));
    RooRealVar Nbkg("Nbkg","Bernstein Background Yield", 1993.8, 0, 968000000);


    // ERFC background (for left sideband)
    RooRealVar csf("csf", "Shifting Constant", 5.14159, 5.08, 5.16);
    RooRealVar csc("csc", "Scaling Constant", 0.03528, 0.0008, 0.05);

    // integral form implemented via erf: 1 - erf(x)
    RooGenericPdf erfc_bkg("erfc_bkg", "1 - TMath::Erf((B_mass - csf)/csc)", RooArgList(B_mass, csf, csc));


    RooRealVar Nerfc("Nerfc", "ERFC Background Yield", 697.3, 0, 50000000);
    
    RooAddPdf model("model", "Signal + Background", RooArgList(signal, bern2, erfc_bkg), RooArgList(Nsig, Nbkg, Nerfc));

    // Fit the model to data (Extended Maximum Likelihood)
    RooFitResult* result = model.fitTo(dataset, Save(), Range(xlow, xhigh), Extended(kTRUE));

    // ---------- Canvas with two pads ----------
    TCanvas* c = new TCanvas("c", "Bmass Fit with Pulls (ERFC model)", 800, 800);
    c->Divide(1, 2);

    // ---------- Top pad (fit) ----------
    TPad* p1 = (TPad*)c->cd(1);
    p1->SetPad(0.0, 0.15, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->Draw();
    p1->cd();

    RooPlot* frame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));

    // Plot data + model with the same naming/styles used for pulls
    dataset.plotOn(frame, Binning(nbins_plot), MarkerStyle(20), MarkerSize(1.2), Name("data"), DataError(RooAbsData::Poisson));
    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global")); // total
    model.plotOn(frame, Components(bern2), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("bern_bkg")); // bernstein part
    model.plotOn(frame, Components(erfc_bkg), LineColor(kMagenta), LineStyle(kDotted), LineWidth(2), Name("erfc_bkg")); // erfc part
    model.plotOn(frame, Components(signal),   LineColor(kGreen+2), LineStyle(kDashed), LineWidth(2), Name("signal")); // signal

    frame->SetTitle("");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0);  // hide x labels on top pad
    frame->GetXaxis()->SetTitle("m_{J/#Psi K^{+}} [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle(Form("Events / ( %.4f )", bin_width_plot));
    frame->Draw();


    // Calculate chi2/ndf for the fit
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare("global", "data", nParams);

    // ---------- Legend (same place), on TOP pad ----------
    p1->cd();
    TLegend* legend = new TLegend(0.52, 0.68, 0.92, 0.96);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetFillStyle(1001);   // Solid fill
    legend->SetFillColor(kWhite); // White background (you can pick another color)
    legend->AddEntry(frame->findObject("data"), "Data (B^{+}) Unbinned", "lep");
    legend->AddEntry(frame->findObject("bern_bkg"), "Background Fit (2nd Degree Poly)", "l");
    legend->AddEntry(frame->findObject("erfc_bkg"), "#splitline{Partially Reconstructed}{Background Fit (Erfc)}","l");
    legend->AddEntry(frame->findObject("signal"), "Signal Fit (Double Gaussian)", "l");
    legend->AddEntry(frame->findObject("global"), "Total Fit (Signal + Background)", "l");
    legend->Draw();

    // ---------- TPaveText (same place), on TOP pad ----------
    p1->cd();
    TPaveText* pave = new TPaveText(0.60, 0.24, 0.92, 0.68, "NDC");
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
    // Bernstein-2 background
    pave->AddText(Form("c_{0} = %.4f #pm %.4f", c0.getVal(), c0.getError()));
    pave->AddText(Form("c_{1} = %.4f #pm %.4f", c1b.getVal(), c1b.getError()));
    pave->AddText(Form("c_{2} = %.4f #pm %.4f", c2.getVal(), c2.getError()));
    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    // ERFC background
    pave->AddText(Form("c_{sf} = %.5f #pm %.5f", csf.getVal(), csf.getError()));
    pave->AddText(Form("c_{sc} = %.5f #pm %.5f", csc.getVal(), csc.getError()));
    pave->AddText(Form("N_{erfc} = %.1f #pm %.1f", Nerfc.getVal(), Nerfc.getError()));
    // Chi2
    pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2));
    pave->Draw();


    // ---------- Bottom pad (pulls) ----------
    TPad* p2 = (TPad*)c->cd(2);
    p2->SetPad(0.0, 0.0, 1.0, 0.15);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.25);
    p2->Draw();
    p2->cd();

    RooPlot* pullFrame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));
    RooHist* pullHist = frame->pullHist("data", "global");  // names must match Name("data") and Name("global")
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
    TString name_file = "Bu_Total_Fit_2degreepoly.pdf";
    c->SaveAs(name_file);


    delete zeroLine;
    delete c;
}
*/




/*
void Bplus_fixedmean() {
    const int nbins_plot = 100; // Number of bins for the plot

    double min_signal = 5.178948768;
    double max_signal = 5.380091232;

    double mc_mean = 5.27952;
    double mc_sigma1 = 0.03702;
    double mc_sigma2 = 0.01609;
    double mc_c1 = 0.3358;
    

    double xlow = 5.0;
    double xhigh = 5.8;

    double bin_width_plot = (xhigh - xlow) / nbins_plot;

    // Load ROOT file and TTree
    TFile *file = TFile::Open("data_unbinned_Bu_Opticut.root");


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

    // Define regions for integrals
    B_mass.setRange("signalRegion", min_signal, max_signal);
    B_mass.setRange("lowSideband", xlow, min_signal);
    B_mass.setRange("highSideband", max_signal, xhigh);


    // Signal model: Double Gaussian
    RooRealVar mean("mean", "Mean", mc_mean);
    mean.setConstant(kTRUE);


    // MC-derived widths (FIXED constants — put your values here)
    RooRealVar sigma1_mc("sigma1_mc", "MC Sigma1", mc_sigma1); 
    sigma1_mc.setConstant(kTRUE);

    RooRealVar sigma2_mc("sigma2_mc", "MC Sigma2", mc_sigma2); 
    sigma2_mc.setConstant(kTRUE);

    RooRealVar c1("c1", "Fraction of Gaussian1", mc_c1);
    c1.setConstant(kTRUE);


    // Common positive scale (fit parameter)
    RooRealVar Cs("Cs", "Resolution scale", 1.028, 0.2, 3.0);

    // Effective widths = Cs * sigma_mc
    RooProduct sigma1_eff("sigma1_eff", "sigma1_eff", RooArgList(sigma1_mc, Cs));
    RooProduct sigma2_eff("sigma2_eff", "sigma2_eff", RooArgList(sigma2_mc, Cs));

    // Mixture fraction and Gaussians
    RooGaussian gauss1("gauss1", "Gaussian 1", B_mass, mean, sigma1_eff);
    RooGaussian gauss2("gauss2", "Gaussian 2", B_mass, mean, sigma2_eff);
    RooAddPdf signal("signal", "Double Gaussian Model", RooArgList(gauss1, gauss2), RooArgList(c1));
    RooRealVar Nsig("Nsig", "Signal Yield", 30000, 0, 96300000);
    

    // Background model
    RooRealVar lambda("lambda", "Lambda", -0.73, -6.32, -0.01);
    RooExponential expo("expo", "Exponential Background", B_mass, lambda);
    RooRealVar Nbkg("Nbkg", "Exponential Background Yield", 35700, 0, 968000000);


    // ERFC background (for left sideband)
    RooRealVar csf("csf", "Shifting Constant", 5.1364, 5.08, 5.16);
    RooRealVar csc("csc", "Scaling Constant", 0.0282, 0.0008, 0.05);

    // integral form implemented via erf: 1 - erf(x)
    RooGenericPdf erfc_bkg("erfc_bkg", "1 - TMath::Erf((B_mass - csf)/csc)", RooArgList(B_mass, csf, csc));


    RooRealVar Nerfc("Nerfc", "ERFC Background Yield", 6132, 0, 500000);
    
    RooAddPdf model("model", "Signal + Background", RooArgList(signal, expo, erfc_bkg), RooArgList(Nsig, Nbkg, Nerfc));


    // Fit the model to data (Extended Maximum Likelihood)
    RooFitResult* result = model.fitTo(dataset, Save(), Range(xlow, xhigh), Extended(kTRUE));

    // ---------- Canvas with two pads ----------
    TCanvas* c = new TCanvas("c", "Bmass Fit with Pulls (ERFC model)", 800, 800);
    c->Divide(1, 2);

    // ---------- Top pad (fit) ----------
    TPad* p1 = (TPad*)c->cd(1);
    p1->SetPad(0.0, 0.15, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->Draw();
    p1->cd();

    RooPlot* frame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));

    // Plot data + model with the same naming/styles used for pulls
    dataset.plotOn(frame, Binning(nbins_plot), MarkerStyle(20), MarkerSize(1.2), Name("data"), DataError(RooAbsData::Poisson));
    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global")); // total
    model.plotOn(frame, Components(expo), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("background")); // exponential part
    model.plotOn(frame, Components(erfc_bkg), LineColor(kMagenta), LineStyle(kDotted), LineWidth(2), Name("erfc_bkg")); // erfc part
    model.plotOn(frame, Components(signal),   LineColor(kGreen+2), LineStyle(kDashed), LineWidth(2), Name("signal")); // signal

    frame->SetTitle("");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0);  // hide x labels on top pad
    frame->GetXaxis()->SetTitle("m_{J/#Psi K^{+}} [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle(Form("Events / ( %.4f )", bin_width_plot));
    frame->Draw();


    // Calculate chi2/ndf for the fit
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare("global", "data", nParams);

    // ---------- Legend (same place), on TOP pad ----------
    p1->cd();
    TLegend* legend = new TLegend(0.52, 0.68, 0.92, 0.96);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetFillStyle(1001);   // Solid fill
    legend->SetFillColor(kWhite); // White background (you can pick another color)
    legend->AddEntry(frame->findObject("data"), "Data (B^{+}) Unbinned", "lep");
    legend->AddEntry(frame->findObject("background"), "Background Fit (Exponential)", "l");
    legend->AddEntry(frame->findObject("erfc_bkg"), "#splitline{Partially Reconstructed}{Background Fit (Erfc)}","l");
    legend->AddEntry(frame->findObject("signal"), "Signal Fit (Double Gaussian)", "l");
    legend->AddEntry(frame->findObject("global"), "Total Fit (Signal + Background)", "l");
    legend->Draw();

    // ---------- TPaveText (same place), on TOP pad ----------
    p1->cd();
    TPaveText* pave = new TPaveText(0.60, 0.24, 0.92, 0.68, "NDC");
    pave->SetTextAlign(12);
    pave->SetTextFont(42);
    pave->SetTextSize(0.025);
    pave->SetFillColor(0);
    pave->SetBorderSize(1);
    // Signal: Double Gaussian
    pave->AddText(Form("Mean (fixed) = %.5f", mean.getVal()));
    pave->AddText(Form("#sigma_{1} (fixed) = %.5f", sigma1_mc.getVal()));
    pave->AddText(Form("#sigma_{2} (fixed) = %.5f", sigma2_mc.getVal()));
    pave->AddText(Form("c_{1} (fixed) = %.4f", c1.getVal()));
    pave->AddText(Form("C_{s} = %.5f #pm %.5f", Cs.getVal(), Cs.getError()));
    pave->AddText(Form("N_{sig} = %.1f #pm %.1f", Nsig.getVal(), Nsig.getError()));
    // Exponential background
    pave->AddText(Form("#lambda = %.4f #pm %.4f", lambda.getVal(), lambda.getError()));
    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    // ERFC background
    pave->AddText(Form("c_{sf} = %.5f #pm %.5f", csf.getVal(), csf.getError()));
    pave->AddText(Form("c_{sc} = %.5f #pm %.5f", csc.getVal(), csc.getError()));
    pave->AddText(Form("N_{erfc} = %.1f #pm %.1f", Nerfc.getVal(), Nerfc.getError()));
    // Chi2
    pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2));
    pave->Draw();


    // ---------- Bottom pad (pulls) ----------
    TPad* p2 = (TPad*)c->cd(2);
    p2->SetPad(0.0, 0.0, 1.0, 0.15);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.25);
    p2->Draw();
    p2->cd();

    RooPlot* pullFrame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));
    RooHist* pullHist = frame->pullHist("data", "global");  // names must match Name("data") and Name("global")
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
    TString name_file = "Bu_Total_Fit_FixedMean.pdf";
    c->SaveAs(name_file);


    delete zeroLine;
    delete c;
}
*/






void Bplus_triplegaussian() {
    const int nbins_plot = 100; // Number of bins for the plot

    double min_signal = 5.178948768;
    double max_signal = 5.380091232;

    double mc_mean = 5.27952;
    double mc_sigma1 = 0.03702;
    double mc_sigma2 = 0.01609;
    double mc_c1 = 0.3358;
    

    double xlow = 5.0;
    double xhigh = 5.8;

    double bin_width_plot = (xhigh - xlow) / nbins_plot;

    // Load ROOT file and TTree
    TFile *file = TFile::Open("data_unbinned_Bu_Opticut.root");


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

    // Define regions for integrals
    B_mass.setRange("signalRegion", min_signal, max_signal);
    B_mass.setRange("lowSideband", xlow, min_signal);
    B_mass.setRange("highSideband", max_signal, xhigh);

    // Signal model: Triple Gaussian (all free, no MC scaling)
    RooRealVar mean("mean", "Mean", 5.27842, 5.27, 5.29);

    // Free widths (positive)
    RooRealVar sigma1("sigma1", "Sigma1", 0.035, 0.003, 0.045);
    RooRealVar sigma2("sigma2", "Sigma2", 0.022, 0.003, 0.035);
    RooRealVar sigma3("sigma3", "Sigma3", 0.014, 0.003, 0.035);

    // Components
    RooGaussian gauss1("gauss1", "Gaussian 1", B_mass, mean, sigma1);
    RooGaussian gauss2("gauss2", "Gaussian 2", B_mass, mean, sigma2);
    RooGaussian gauss3("gauss3", "Gaussian 3", B_mass, mean, sigma3);

    RooRealVar r1("r1","r1", 2.0, 0.1, 10.0);
    RooRealVar r2("r2","r2", 2.0, 0.1, 10.0);
    RooFormulaVar f1("f1","@0/(1+@0+@1)",RooArgList(r1,r2));
    RooFormulaVar f2("f2","@1/(1+@0+@1)",RooArgList(r1,r2));
    RooAddPdf signal("signal", "Triple Gaussian Model", RooArgList(gauss1, gauss2, gauss3), RooArgList(f1, f2));


    RooRealVar Nsig("Nsig", "Signal Yield", 30200, 100, 100000);


    // Background model
    RooRealVar lambda("lambda", "Lambda", -0.727, -4.32, -0.01);
    RooExponential expo("expo", "Exponential Background", B_mass, lambda);
    RooRealVar Nbkg("Nbkg", "Exponential Background Yield", 35700, 1000, 100000);


    // ERFC background (for left sideband)
    RooRealVar csf("csf", "Shifting Constant", 5.1364, 5.08, 5.16);
    RooRealVar csc("csc", "Scaling Constant", 0.028, 0.0008, 0.05);

    // integral form implemented via erf: 1 - erf(x)
    RooGenericPdf erfc_bkg("erfc_bkg", "1 - TMath::Erf((B_mass - csf)/csc)", RooArgList(B_mass, csf, csc));


    RooRealVar Nerfc("Nerfc", "ERFC Background Yield", 6132, 100, 15000);
    
    RooAddPdf model("model", "Signal + Background", RooArgList(signal, expo, erfc_bkg), RooArgList(Nsig, Nbkg, Nerfc));


    // Fit the model to data (Extended Maximum Likelihood)
    RooFitResult* result = model.fitTo(dataset, Save(), Range(xlow, xhigh), Extended(kTRUE));

    // Fractions and their propagated uncertainties
    double f1_val = f1.getVal();
    double f1_err = f1.getPropagatedError(*result);
    double f2_val = f2.getVal();
    double f2_err = f2.getPropagatedError(*result);


    // ---------- Canvas with two pads ----------
    TCanvas* c = new TCanvas("c", "Bmass Fit with Pulls (ERFC model)", 800, 800);
    c->Divide(1, 2);

    // ---------- Top pad (fit) ----------
    TPad* p1 = (TPad*)c->cd(1);
    p1->SetPad(0.0, 0.15, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->Draw();
    p1->cd();

    RooPlot* frame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));

    // Plot data + model with the same naming/styles used for pulls
    dataset.plotOn(frame, Binning(nbins_plot), MarkerStyle(20), MarkerSize(1.2), Name("data"), DataError(RooAbsData::Poisson));
    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global")); // total
    model.plotOn(frame, Components(expo), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("background")); // exponential part
    model.plotOn(frame, Components(erfc_bkg), LineColor(kMagenta), LineStyle(kDotted), LineWidth(2), Name("erfc_bkg")); // erfc part
    model.plotOn(frame, Components(signal),   LineColor(kGreen+2), LineStyle(kDashed), LineWidth(2), Name("signal")); // signal

    frame->SetTitle("");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0);  // hide x labels on top pad
    frame->GetXaxis()->SetTitle("m_{J/#Psi K^{+}} [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle(Form("Events / ( %.4f )", bin_width_plot));
    frame->Draw();


    // Calculate chi2/ndf for the fit
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare("global", "data", nParams);

    // ---------- Legend (same place), on TOP pad ----------
    p1->cd();
    TLegend* legend = new TLegend(0.52, 0.68, 0.92, 0.96);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetFillStyle(1001);   // Solid fill
    legend->SetFillColor(kWhite); // White background (you can pick another color)
    legend->AddEntry(frame->findObject("data"), "Data (B^{+}) Unbinned", "lep");
    legend->AddEntry(frame->findObject("background"), "Background Fit (Exponential)", "l");
    legend->AddEntry(frame->findObject("erfc_bkg"), "#splitline{Partially Reconstructed}{Background Fit (Erfc)}","l");
    legend->AddEntry(frame->findObject("signal"), "Signal Fit (Triple Gaussian)", "l");
    legend->AddEntry(frame->findObject("global"), "Total Fit (Signal + Background)", "l");
    legend->Draw();

    // ---------- TPaveText (same place), on TOP pad ----------
    p1->cd();
    TPaveText* pave = new TPaveText(0.60, 0.24, 0.92, 0.68, "NDC");
    pave->SetTextAlign(12);
    pave->SetTextFont(42);
    pave->SetTextSize(0.025);
    pave->SetFillColor(0);
    pave->SetBorderSize(1);
    // Signal: Double Gaussian
    pave->AddText(Form("Mean = %.5f #pm %.5f", mean.getVal(), mean.getError()));
    pave->AddText(Form("#sigma_{1} = %.5f #pm %.5f", sigma1.getVal(), sigma1.getError()));
    pave->AddText(Form("#sigma_{2} = %.5f #pm %.5f", sigma2.getVal(), sigma2.getError()));
    pave->AddText(Form("#sigma_{3} = %.5f #pm %.5f", sigma3.getVal(), sigma3.getError()));
    pave->AddText(Form("f_{1} = %.4f #pm %.4f", f1_val, f1_err));
    pave->AddText(Form("f_{2} = %.4f #pm %.4f", f2_val, f2_err));
    pave->AddText(Form("N_{sig} = %.1f #pm %.1f", Nsig.getVal(), Nsig.getError()));
    // Exponential background
    pave->AddText(Form("#lambda = %.4f #pm %.4f", lambda.getVal(), lambda.getError()));
    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    // ERFC background
    pave->AddText(Form("c_{sf} = %.5f #pm %.5f", csf.getVal(), csf.getError()));
    pave->AddText(Form("c_{sc} = %.5f #pm %.5f", csc.getVal(), csc.getError()));
    pave->AddText(Form("N_{erfc} = %.1f #pm %.1f", Nerfc.getVal(), Nerfc.getError()));
    // Chi2
    pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2));
    pave->Draw();


    // ---------- Bottom pad (pulls) ----------
    TPad* p2 = (TPad*)c->cd(2);
    p2->SetPad(0.0, 0.0, 1.0, 0.15);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.25);
    p2->Draw();
    p2->cd();

    RooPlot* pullFrame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));
    RooHist* pullHist = frame->pullHist("data", "global");  // names must match Name("data") and Name("global")
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
    TString name_file = "Bu_Total_Fit_TripleGaussian.pdf";
    c->SaveAs(name_file);


    delete zeroLine;
    delete c;
}






/*
void Bplus_crystalball() {
    const int nbins_plot = 100; // Number of bins for the plot

    double min_signal = 5.178948768;
    double max_signal = 5.380091232;

    double mc_mean = 5.27952;
    double mc_sigma1 = 0.03702;
    double mc_sigma2 = 0.01609;
    double mc_c1 = 0.3358;
    

    double xlow = 5.0;
    double xhigh = 5.8;

    double bin_width_plot = (xhigh - xlow) / nbins_plot;

    // Load ROOT file and TTree
    TFile *file = TFile::Open("data_unbinned_Bu_Opticut.root");


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

    // Define regions for integrals
    B_mass.setRange("signalRegion", min_signal, max_signal);
    B_mass.setRange("lowSideband", xlow, min_signal);
    B_mass.setRange("highSideband", max_signal, xhigh);


    // Signal model: Double Gaussian
    RooRealVar mean("mean", "Mean", 5.27764, 5.27, 5.29);


    // MC-derived widths (FIXED constants — put your values here)
    RooRealVar sigma1_mc("sigma1_mc", "MC Sigma1", mc_sigma1); 
    sigma1_mc.setConstant(kTRUE);

    RooRealVar sigma2_mc("sigma2_mc", "MC Sigma2", mc_sigma2); 
    sigma2_mc.setConstant(kTRUE);

    RooRealVar c1("c1", "Fraction of Gaussian1", 0.75, 0.50, 0.98);


    // Common positive scale (fit parameter)
    RooRealVar Cs("Cs", "Resolution scale", 1.15888, 0.2, 3.0);

    // Effective widths = Cs * sigma_mc
    RooProduct sigma1_eff("sigma1_eff", "sigma1_eff", RooArgList(sigma1_mc, Cs));
    RooProduct sigma2_eff("sigma2_eff", "sigma2_eff", RooArgList(sigma2_mc, Cs));
    

    // Crystal Ball tail parameters (left tail typical for mass peaks)
    RooRealVar alpha_cb("alpha_cb", "CB alpha", -2.0, -5.0, -0.8);
    RooRealVar n_cb("n_cb", "CB n", 5.0, 1.975, 50.0);


    // Mixture fraction and Gaussians
    RooGaussian gauss1("gauss1", "Gaussian 1", B_mass, mean, sigma1_eff);
    RooCBShape cb("cb", "Crystal Ball", B_mass, mean, sigma2_eff, alpha_cb, n_cb);
    RooAddPdf signal("signal", "Gauss + Crystal Ball Model", RooArgList(gauss1, cb), RooArgList(c1));
    RooRealVar Nsig("Nsig", "Signal Yield", 3441.2, 0, 96300000);
    

    // Background model
    RooRealVar lambda("lambda", "Lambda", -2.1699, -6.32, -0.01);
    RooExponential expo("expo", "Exponential Background", B_mass, lambda);
    RooRealVar Nbkg("Nbkg", "Exponential Background Yield", 1993.8, 0, 968000000);


    // ERFC background (for left sideband)
    RooRealVar csf("csf", "Shifting Constant", 5.14159, 5.08, 5.16);
    RooRealVar csc("csc", "Scaling Constant", 0.03528, 0.0008, 0.1);

    // integral form implemented via erf: 1 - erf(x)
    RooGenericPdf erfc_bkg("erfc_bkg", "1 - TMath::Erf((B_mass - csf)/csc)", RooArgList(B_mass, csf, csc));


    RooRealVar Nerfc("Nerfc", "ERFC Background Yield", 697.3, 0, 50000000);
    
    RooAddPdf model("model", "Signal + Background", RooArgList(signal, expo, erfc_bkg), RooArgList(Nsig, Nbkg, Nerfc));


    // Fit the model to data (Extended Maximum Likelihood)
    RooFitResult* result = model.fitTo(dataset, Save(), Range(xlow, xhigh), Extended(kTRUE));

    // ---------- Canvas with two pads ----------
    TCanvas* c = new TCanvas("c", "Bmass Fit with Pulls (ERFC model)", 800, 800);
    c->Divide(1, 2);

    // ---------- Top pad (fit) ----------
    TPad* p1 = (TPad*)c->cd(1);
    p1->SetPad(0.0, 0.15, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->Draw();
    p1->cd();

    RooPlot* frame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));

    // Plot data + model with the same naming/styles used for pulls
    dataset.plotOn(frame, Binning(nbins_plot), MarkerStyle(20), MarkerSize(1.2), Name("data"), DataError(RooAbsData::Poisson));
    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global")); // total
    model.plotOn(frame, Components(expo), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("background")); // exponential part
    model.plotOn(frame, Components(erfc_bkg), LineColor(kMagenta), LineStyle(kDotted), LineWidth(2), Name("erfc_bkg")); // erfc part
    model.plotOn(frame, Components(signal),   LineColor(kGreen+2), LineStyle(kDashed), LineWidth(2), Name("signal")); // signal

    frame->SetTitle("");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0);  // hide x labels on top pad
    frame->GetXaxis()->SetTitle("m_{J/#Psi K^{+}} [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle(Form("Events / ( %.4f )", bin_width_plot));
    frame->Draw();


    // Calculate chi2/ndf for the fit
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare("global", "data", nParams);

    // ---------- Legend (same place), on TOP pad ----------
    p1->cd();
    TLegend* legend = new TLegend(0.52, 0.68, 0.92, 0.96);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetFillStyle(1001);   // Solid fill
    legend->SetFillColor(kWhite); // White background (you can pick another color)
    legend->AddEntry(frame->findObject("data"), "Data (B^{+}) Unbinned", "lep");
    legend->AddEntry(frame->findObject("background"), "Background Fit (Exponential)", "l");
    legend->AddEntry(frame->findObject("erfc_bkg"), "#splitline{Partially Reconstructed}{Background Fit (Erfc)}","l");
    legend->AddEntry(frame->findObject("signal"), "Signal Fit (Gauss + Crystal Ball)", "l");
    legend->AddEntry(frame->findObject("global"), "Total Fit (Signal + Background)", "l");
    legend->Draw();

    // ---------- TPaveText (same place), on TOP pad ----------
    p1->cd();
    TPaveText* pave = new TPaveText(0.60, 0.24, 0.92, 0.68, "NDC");
    pave->SetTextAlign(12);
    pave->SetTextFont(42);
    pave->SetTextSize(0.025);
    pave->SetFillColor(0);
    pave->SetBorderSize(1);
    // Signal: Double Gaussian
    pave->AddText(Form("Mean = %.5f #pm %.5f", mean.getVal(), mean.getError()));
    pave->AddText(Form("#sigma_{G} (fixed) = %.5f", sigma1_mc.getVal()));
    pave->AddText(Form("#sigma_{CB} (fixed) = %.5f", sigma2_mc.getVal()));
    pave->AddText(Form("#alpha_{CB} = %.3f #pm %.3f", alpha_cb.getVal(), alpha_cb.getError()));
    pave->AddText(Form("n_{CB} = %.3f #pm %.3f", n_cb.getVal(), n_cb.getError()));
    pave->AddText(Form("c_{1} = %.4f #pm %.4f", c1.getVal(), c1.getError()));
    pave->AddText(Form("C_{s} = %.5f #pm %.5f", Cs.getVal(), Cs.getError()));
    pave->AddText(Form("N_{sig} = %.1f #pm %.1f", Nsig.getVal(), Nsig.getError()));
    // Exponential background
    pave->AddText(Form("#lambda = %.4f #pm %.4f", lambda.getVal(), lambda.getError()));
    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    // ERFC background
    pave->AddText(Form("c_{sf} = %.5f #pm %.5f", csf.getVal(), csf.getError()));
    pave->AddText(Form("c_{sc} = %.5f #pm %.5f", csc.getVal(), csc.getError()));
    pave->AddText(Form("N_{erfc} = %.1f #pm %.1f", Nerfc.getVal(), Nerfc.getError()));
    // Chi2
    pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2));
    pave->Draw();


    // ---------- Bottom pad (pulls) ----------
    TPad* p2 = (TPad*)c->cd(2);
    p2->SetPad(0.0, 0.0, 1.0, 0.15);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.25);
    p2->Draw();
    p2->cd();

    RooPlot* pullFrame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));
    RooHist* pullHist = frame->pullHist("data", "global");  // names must match Name("data") and Name("global")
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
    TString name_file = "Bu_Total_Fit_CrystalBall.pdf";
    c->SaveAs(name_file);


    delete zeroLine;
    delete c;
}
*/



/*
void Bplus_Without_Left_Sideband() {
    const int nbins_plot = 100; // Number of bins for the plot

    double min_signal = 5.178948768;
    double max_signal = 5.380091232;

    double mc_mean   = 5.27952;
    double mc_sigma1 = 0.03702;
    double mc_sigma2 = 0.01609;
    double mc_c1     = 0.3358;

    // WITHOUT LEFT SIDEBAND !
    double xhigh = 5.8;
    double xlow  = min_signal - 0.02; // 0.02 to expand a little bit the signal region

    double bin_width_plot = (xhigh - xlow) / nbins_plot;

    // Load ROOT file and TTree
    TFile *file = TFile::Open("data_unbinned_Bu_Opticut.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open real data file." << std::endl;
        return;
    }
    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cerr << "Error: Tree 'tree' not found in file." << std::endl;
        return;
    }

    // Define the mass variable and dataset over the desired fit range
    RooRealVar B_mass("B_mass", "B_mass", xlow, xhigh);

    // Create unbinned RooDataSet from TTree
    B_mass.setRange("gaussRange", min_signal, max_signal);
    RooDataSet dataset("dataset", "Unbinned dataset from TTree", tree, RooArgSet(B_mass));

    // Optional named ranges
    B_mass.setRange("signalRegion", min_signal, max_signal);
    B_mass.setRange("lowSideband", xlow,        min_signal);
    B_mass.setRange("highSideband", max_signal, xhigh);

    // -------------------- Signal: Double Gaussian (MC-scaled widths) --------------------
    RooRealVar mean("mean", "Mean", 5.27764, 5.27, 5.29);

    // MC-derived widths (fixed) and mixture fraction (fixed)
    RooRealVar sigma1_mc("sigma1_mc", "MC Sigma1", mc_sigma1); sigma1_mc.setConstant(kTRUE);
    RooRealVar sigma2_mc("sigma2_mc", "MC Sigma2", mc_sigma2); sigma2_mc.setConstant(kTRUE);
    RooRealVar c1("c1", "Fraction of Gaussian1", mc_c1);       c1.setConstant(kTRUE);

    // Resolution scale (free)
    RooRealVar Cs("Cs", "Resolution scale", 1.15888, 0.2, 3.0);

    // Effective widths = Cs * sigma_mc
    RooProduct sigma1_eff("sigma1_eff", "sigma1_eff", RooArgList(sigma1_mc, Cs));
    RooProduct sigma2_eff("sigma2_eff", "sigma2_eff", RooArgList(sigma2_mc, Cs));

    // Two Gaussians with shared mean
    RooGaussian gauss1("gauss1", "Gaussian 1", B_mass, mean, sigma1_eff);
    RooGaussian gauss2("gauss2", "Gaussian 2", B_mass, mean, sigma2_eff);

    // Double-Gaussian signal
    RooAddPdf signal("signal", "Double Gaussian Model", RooArgList(gauss1, gauss2), RooArgList(c1));
    RooRealVar Nsig("Nsig", "Signal Yield", 3441.2, 0, 9.63e7);

    // -------------------- Background: Exponential only --------------------
    RooRealVar lambda("lambda", "Lambda", -2.1699, -6.32, -0.01);
    RooExponential expo("expo", "Exponential Background", B_mass, lambda);
    RooRealVar Nbkg("Nbkg", "Exponential Background Yield", 1993.8, 0, 9.68e8);

    // -------------------- Total model (extended) --------------------
    RooAddPdf model("model", "Signal + Background",
                    RooArgList(signal, expo),
                    RooArgList(Nsig, Nbkg));

    // -------------------- Fit --------------------
    RooFitResult* result = model.fitTo(dataset, Save(), Range(xlow, xhigh), Extended(kTRUE));

    // -------------------- Plot: fit + pulls --------------------
    TCanvas* c = new TCanvas("c", "Bmass Fit with Pulls (Exp bkg, DoubleG signal)", 800, 800);
    c->Divide(1, 2);

    // Top pad (fit)
    TPad* p1 = (TPad*)c->cd(1);
    p1->SetPad(0.0, 0.15, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->Draw();
    p1->cd();

    RooPlot* frame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));

    dataset.plotOn(frame, Binning(nbins_plot), MarkerStyle(20), MarkerSize(1.2), Name("data"), DataError(RooAbsData::Poisson));
    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global")); // total
    model.plotOn(frame, Components(expo),   LineColor(kRed),   LineStyle(kDashed), LineWidth(2), Name("background")); // exponential part
    model.plotOn(frame, Components(signal), LineColor(kGreen+2), LineStyle(kDashed), LineWidth(2), Name("signal")); // signal

    frame->SetTitle("");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0);  // hide x labels on top pad
    frame->GetXaxis()->SetTitle("m_{J/#psi K^{+}} [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle(Form("Events / ( %.4f )", bin_width_plot));
    frame->Draw();

    // Chi2/ndf
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare("global", "data", nParams);

    // Legend
    p1->cd();
    TLegend* legend = new TLegend(0.52, 0.68, 0.92, 0.96);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetFillStyle(1001);
    legend->SetFillColor(kWhite);
    legend->AddEntry(frame->findObject("data"),   "Data (B^{+}) Unbinned", "lep");
    legend->AddEntry(frame->findObject("background"), "Background Fit (Exponential)", "l");
    legend->AddEntry(frame->findObject("signal"),  "Signal Fit (Double Gaussian)", "l");
    legend->AddEntry(frame->findObject("global"),  "Total Fit (Signal + Background)", "l");
    legend->Draw();

    // TPaveText with fit params
    p1->cd();
    TPaveText* pave = new TPaveText(0.60, 0.24, 0.92, 0.68, "NDC");
    pave->SetTextAlign(12);
    pave->SetTextFont(42);
    pave->SetTextSize(0.025);
    pave->SetFillColor(0);
    pave->SetBorderSize(1);
    // Signal
    pave->AddText(Form("Mean = %.5f #pm %.5f", mean.getVal(), mean.getError()));
    pave->AddText(Form("#sigma_{1} (fixed) = %.5f", sigma1_mc.getVal()));
    pave->AddText(Form("#sigma_{2} (fixed) = %.5f", sigma2_mc.getVal()));
    pave->AddText(Form("c_{1} (fixed) = %.4f", c1.getVal()));
    pave->AddText(Form("C_{s} = %.5f #pm %.5f", Cs.getVal(), Cs.getError()));
    pave->AddText(Form("N_{sig} = %.1f #pm %.1f", Nsig.getVal(), Nsig.getError()));
    // Background (Exponential only)
    pave->AddText(Form("#lambda = %.4f #pm %.4f", lambda.getVal(), lambda.getError()));
    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    // Chi2
    pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2));
    pave->Draw();

    // Bottom pad (pulls)
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
    pullFrame->GetXaxis()->SetTitle("m_{J/#psi K^{+}} [GeV/c^{2}]");
    pullFrame->GetXaxis()->SetTitleSize(0.10);
    pullFrame->GetXaxis()->SetTitleOffset(1.0);
    pullFrame->GetXaxis()->SetLabelSize(0.08);
    pullFrame->SetMinimum(-3.5);
    pullFrame->SetMaximum(3.5);
    pullFrame->Draw("AP");

    // Zero line on pulls
    TLine* zeroLine = new TLine(xlow, 0, xhigh, 0);
    zeroLine->SetLineColor(kBlue);
    zeroLine->SetLineStyle(1);
    zeroLine->SetLineWidth(1);
    zeroLine->Draw("same");

    // Save
    TString name_file = "Bu_Total_Fit_No_Left_Sideband.pdf";
    c->SaveAs(name_file);

    delete zeroLine;
    delete c;
}
*/



void Bplus_variations() {
    //Bplus_nominal();
    //Bplus_LinearBkg();
    //Bplus_2degreepoly();
    //Bplus_fixedmean();

    
    Bplus_triplegaussian();
    
    
    //Bplus_crystalball();
    //Bplus_Without_Left_Sideband();
}