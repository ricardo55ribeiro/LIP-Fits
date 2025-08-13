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
#include <RooChebychev.h>
#include <RooBinning.h>




// Auxiliar function to compute black dotted lines height
double getYatMass(RooPlot* frame, double mass) {
    for (int i = 0; i < frame->numItems(); ++i) {
        RooCurve* curve = dynamic_cast<RooCurve*>(frame->getObject(i));
        if (!curve) continue;

        int n = curve->GetN();
        double* x = curve->GetX();
        double* y = curve->GetY();

        for (int j = 0; j < n - 1; ++j) {
            if (x[j] <= mass && mass <= x[j+1]) {
                // Linear interpolation
                double slope = (y[j+1] - y[j]) / (x[j+1] - x[j]);
                return y[j] + slope * (mass - x[j]);
            }
        }
    }
    return 0.0;
}



void total_data_fit_PSI_and_X3872(){
    const int nbins_plot = 100; // Number of Bins for the Plot

    double min_signal_psi2s = 3.6648893;
    double max_signal_psi2s = 3.7090307;
    
    double min_signal_x3872 = 3.850404115;
    double max_signal_x3872 = 3.894695885;

    double xlow = 3.6;
    double xhigh = 4.0;

    const double bin_width_plot = (xhigh - xlow) / nbins_plot;

    // Load ROOT file and TTree
    TFile *file = TFile::Open("MC_ppRef_Bmass_PSI2S_X3872_newoptimized.root");
    // data_unbinned_X_FirstCut.root
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
    B_mass.setRange("fitRange", xlow, xhigh);

    RooBinning mainBins(nbins_plot, xlow, xhigh);
    B_mass.setBinning(mainBins, "mainBins");

    RooDataSet dataset("dataset", "Unbinned dataset from TTree", tree, RooArgSet(B_mass));


    // Signal model for PSI 2S
    RooRealVar mean_psi2s("mean_psi2s", "Mean", 3.686, 3.67, 3.70);
    RooRealVar sigma1_psi2s("sigma1_psi2s", "Sigma1", 0.0036, 0.001, 0.1);
    RooRealVar sigma2_psi2s("sigma2_psi2s", "Sigma2", 0.008, 0.001, 0.1);
    RooRealVar c1_psi2s("c1_psi2s", "Fraction of Gaussian1", 0.44, 0.01, 0.99);
    RooGaussian gauss1_psi2s("gauss1_psi2s", "Narrow Gaussian", B_mass, mean_psi2s, sigma1_psi2s);
    RooGaussian gauss2_psi2s("gauss2_psi2s", "Wide Gaussian", B_mass, mean_psi2s, sigma2_psi2s);
    RooAddPdf signal_raw_psi2s("signal_raw_psi2s", "Double Gaussian Model", RooArgList(gauss1_psi2s, gauss2_psi2s), RooArgList(c1_psi2s));

    RooRealVar Nsig_psi2s("Nsig_psi2s", "Signal Yield", 46960, 400, 500000);
    RooExtendPdf signal_ext_psi2s("signal_ext_psi2s", "Extended Gated Signal", signal_raw_psi2s, Nsig_psi2s);


    // Signal model for X(3872): same shape as PSI(2S) (sharing σ1, σ2, c1)
    RooRealVar mean_x("mean_x", "Mean of X(3872)", 3.872, 3.86, 3.92);

    // Reuse ψ(2S) shape parameters: sigma1_psi2s, sigma2_psi2s, c1_psi2s
    RooGaussian gauss1_x("gauss1_x", "Narrow Gaussian X", B_mass, mean_x, sigma1_psi2s);
    RooGaussian gauss2_x("gauss2_x", "Wide Gaussian X",   B_mass, mean_x, sigma2_psi2s);
    RooAddPdf   signal_raw_x("signal_raw_x", "Double Gaussian X (shared shape)",
                            RooArgList(gauss1_x, gauss2_x), RooArgList(c1_psi2s));

    RooRealVar Nsig_x("Nsig_x", "Yield X(3872)", 3458, 30, 400000);
    RooExtendPdf signal_ext_x("signal_ext_x", "Extended X(3872)", signal_raw_x, Nsig_x);


    // Background model: Exponential
    RooRealVar a0("a0", "Cheb0", 0.072, -1.0, 1.0);
    RooRealVar a1("a1", "Cheb1", -0.095, -1.0, 1.0);
    RooRealVar a2("a2", "Cheb2", 0.017, -1.0, 1.0);
    RooRealVar a3("a3", "Cheb3", 0, -1.0, 1.0);
    RooRealVar a4("a4", "Cheb4", 0, -1.0, 1.0);

    RooChebychev bkg_pdf("bkg_pdf", "4th-degree background", B_mass, RooArgList(a0, a1, a2, a3, a4));

    RooRealVar Nbkg("Nbkg", "Total background yield", 678752, 6000, 7000000);
    RooExtendPdf expo_ext("expo_ext", "Extended 2-exp background", bkg_pdf, Nbkg);


    // Full model = signal + exponential background
    RooAddPdf model("model", "PSI(2S) + X(3872) + Background", RooArgList(signal_ext_psi2s, signal_ext_x, expo_ext));

    // Fit the model to data
    //RooFitResult* result = model.fitTo(dataset, Save(), Range("fitRange"), NormRange("fitRange"));

    //model.setNormRange("fitRange");  // Set normalization range explicitly
    RooFitResult* result = model.fitTo(dataset, RooFit::Save(), RooFit::Range(xlow, xhigh));




    // ------------------------------------------------------------------------------------------------------------

    // Define regions for integrals
    B_mass.setRange("signalRegion_psi2s", min_signal_psi2s, max_signal_psi2s);
    B_mass.setRange("signalRegion_x3872", min_signal_x3872, max_signal_x3872);
    B_mass.setRange("firstSideband", xlow, min_signal_psi2s);
    B_mass.setRange("secondSideband", max_signal_psi2s, min_signal_x3872);
    B_mass.setRange("thirdSideband", max_signal_x3872, xhigh);

    // Compute background integrals in regions (exponential + erfc)
    double frac_r4 = bkg_pdf.createIntegral(B_mass, NormSet(B_mass), Range("signalRegion_psi2s"))->getVal();
    double frac_r5 = bkg_pdf.createIntegral(B_mass, NormSet(B_mass), Range("signalRegion_x3872"))->getVal();
    double frac_r1 = bkg_pdf.createIntegral(B_mass, NormSet(B_mass), Range("firstSideband"))->getVal();
    double frac_r2 = bkg_pdf.createIntegral(B_mass, NormSet(B_mass), Range("secondSideband"))->getVal();
    double frac_r3 = bkg_pdf.createIntegral(B_mass, NormSet(B_mass), Range("thirdSideband"))->getVal();

    double total_bkg = Nbkg.getVal();
    double r4 = total_bkg * frac_r4;
    double r5 = total_bkg * frac_r5;
    double r1 = total_bkg * frac_r1;
    double r2 = total_bkg * frac_r2;
    double r3 = total_bkg * frac_r3;

    double f_b = (r4 + r5) / (r1 + r2 + r3);
    //double f_b_x3872 = r5 / (r1 + r2 + r3);


    // Compute signal yield in signal region
    double frac_Sdata_psi2s = signal_raw_psi2s.createIntegral(B_mass, NormSet(B_mass), Range("signalRegion_psi2s"))->getVal();
    double Sdata_psi2s = Nsig_psi2s.getVal() * frac_Sdata_psi2s;

    double frac_Sdata_x3872 = signal_raw_x.createIntegral(B_mass, NormSet(B_mass), Range("signalRegion_x3872"))->getVal();
    double Sdata_x3872 = Nsig_x.getVal() * frac_Sdata_x3872;


    // MC for PSI 2S
    TFile *file_mc_psi2s = TFile::Open("/lstore/cms/hlegoinha/X3872/MC_DATA/prompt_PSI2S_to_Jpsi_pipi_phat5_Bfinder.root");
    if (!file_mc_psi2s || file_mc_psi2s->IsZombie()) {
        std::cerr << "Error: Could not open MC file." << std::endl;
        return;
    }
    TTree *treemc_psi2s = nullptr;
    file_mc_psi2s->GetObject("Bfinder/ntmix", treemc_psi2s);
    if (!treemc_psi2s) {
        std::cerr << "Error: MC TTree not found!" << std::endl;
        file_mc_psi2s->Close();
        return;
    }
    TString cut_mc_psi2s = Form("Bchi2cl>0.003 && BQvalueuj<0.14274 && (%s) && (%s) && (%s) && (%s)",
                        isMCsignal.Data(),
                        ACCcuts_ppRef.Data(),
                        SELcuts_ppRef.Data(),
                        TRGmatching.Data());
    int nbins_mc_psi2s = 150;
    TH1F *hist_mc_psi2s = new TH1F("hist_mc_psi2s", "MC Bmass in Signal Region; Bmass [GeV/c^{2}]; Entries", nbins_mc_psi2s, min_signal_psi2s, max_signal_psi2s);
    treemc_psi2s->Draw("Bmass >> hist_mc_psi2s", cut_mc_psi2s + Form(" && Bmass > %.4f && Bmass < %.4f", min_signal_psi2s, max_signal_psi2s), "goff");
    double S_MC_psi2s = hist_mc_psi2s->Integral();
    delete hist_mc_psi2s;
    file_mc_psi2s->Close();


    // MC for X3872
    TFile *file_mc_x3872 = TFile::Open("/lstore/cms/hlegoinha/X3872/MC_DATA/prompt_X3872_to_Jpsi_Rho_phat5_Bfinder.root");
    if (!file_mc_x3872 || file_mc_x3872->IsZombie()) {
        std::cerr << "Error: Could not open MC file." << std::endl;
        return;
    }
    TTree *treemc_x3872 = nullptr;
    file_mc_x3872->GetObject("Bfinder/ntmix", treemc_x3872);
    if (!treemc_x3872) {
        std::cerr << "Error: MC TTree not found!" << std::endl;
        file_mc_x3872->Close();
        return;
    }
    TString cut_mc_x3872 = Form("Bchi2cl>0.02 && BQvalueuj<0.2 && (%s) && (%s) && (%s) && (%s)",
                        isMCsignal.Data(),
                        ACCcuts_ppRef.Data(),
                        SELcuts_ppRef.Data(),
                        TRGmatching.Data());
    int nbins_mc_x3872 = 150;
    TH1F *hist_mc_x3872 = new TH1F("hist_mc_x3872", "MC Bmass in Signal Region; Bmass [GeV/c^{2}]; Entries", nbins_mc_x3872, min_signal_x3872, max_signal_x3872);
    treemc_x3872->Draw("Bmass >> hist_mc_x3872", cut_mc_x3872 + Form(" && Bmass > %.4f && Bmass < %.4f", min_signal_x3872, max_signal_x3872), "goff");
    double S_MC_x3872 = hist_mc_x3872->Integral();
    delete hist_mc_x3872;
    file_mc_x3872->Close();

    double f_s_psi2s = Sdata_psi2s / S_MC_psi2s;
    double f_s_x3872 = Sdata_x3872 / S_MC_x3872;
    


    // ------------------------------------------------------------------------------------------------------------

    // Create canvas with two pads
    TCanvas* c = new TCanvas("c", "Fit with Pulls", 800, 800);
    c->Divide(1, 2);

    // Change size Fit/Pulls is no SetPad(x1, y1, x2, y2)
    // p1 is Fit ; p2 is Pulls

    TPad* p1 = (TPad*)c->cd(1);
    p1->SetPad(0.0, 0.15, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->Draw();

    TPad* p2 = (TPad*)c->cd(2);
    p2->SetPad(0.0, 0.0, 1.0, 0.15);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.25);
    p2->Draw();

    // Plotting the Fit
    p1->cd();
    RooPlot* frame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));
    dataset.plotOn(frame, Binning(B_mass.getBinning("mainBins")), MarkerStyle(20), MarkerSize(1.2), Name("data"), DataError(RooAbsData::Poisson));


    // Draw components
    model.plotOn(frame, Components(signal_ext_psi2s), LineColor(kGreen+2), LineStyle(kDashed), LineWidth(2), Name("psi2s"));
    model.plotOn(frame, Components(signal_ext_x),     LineColor(kCyan+2),  LineStyle(kDashed), LineWidth(2), Name("x3872"));
    model.plotOn(frame, Components(bkg_pdf),          LineColor(kRed),     LineStyle(kDashed), LineWidth(2), Name("poly4"));
    

    // Draw the total model 
    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global"), Range("fitRange"), NormRange("fitRange"));


    frame->SetTitle("");
    frame->GetYaxis()->SetTitle(Form("Events / ( %.4f )", bin_width_plot));
    frame->GetYaxis()->SetTitleOffset(1.48);
    frame->GetXaxis()->SetLabelSize(0);  // hide x labels
    frame->Draw();



    // After frame->Draw();
    p1->Update(); // make sure pad limits are updated

    double y_bottom = gPad->GetUymin(); // bottom visible Y coordinate

    double y_top    = gPad->GetUymax(); // top visible Y coordinate

    TLine* line_min_psi2s = new TLine(min_signal_psi2s, y_bottom, min_signal_psi2s, y_top);
    TLine* line_max_psi2s = new TLine(max_signal_psi2s, y_bottom, max_signal_psi2s, y_top);
    TLine* line_min_x3872 = new TLine(min_signal_x3872, y_bottom, min_signal_x3872, y_top);
    TLine* line_max_x3872 = new TLine(max_signal_x3872, y_bottom, max_signal_x3872, y_top);

    for (TLine* l : {line_min_psi2s, line_max_psi2s, line_min_x3872, line_max_x3872}) {
        l->SetLineColor(kBlack);
        l->SetLineStyle(3);
        l->SetLineWidth(2);
        l->Draw("same");
    }

    // -------- Move the legend back to the top pad ----------
    TLegend* legend = new TLegend(0.35, 0.7, 0.63, 0.97);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetFillColor(0);
    legend->SetFillStyle(1001);

    legend->AddEntry(frame->findObject("data"),   "Unbinned Data", "lep");
    legend->AddEntry(frame->findObject("psi2s"),  "#psi(2S) Signal", "l");
    legend->AddEntry(frame->findObject("x3872"),  "X(3872) Signal",  "l");
    legend->AddEntry(frame->findObject("poly4"),  "Background", "l");
    legend->AddEntry(frame->findObject("global"), "Total Fit", "l");
    legend->Draw();


    // ---------- Plotting the Pulls ----------
    p2->cd();
    RooPlot* pullFrame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));
    RooHist* pullHist = frame->pullHist("data", "global");
    pullHist->SetMarkerSize(0.6);
    pullFrame->addPlotable(pullHist, "XP");

    pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetNdivisions(505);
    pullFrame->GetYaxis()->SetTitleSize(0.1);
    pullFrame->GetYaxis()->SetTitleOffset(0.4);
    pullFrame->GetYaxis()->SetLabelSize(0.08);
    pullFrame->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV/c^{2}]");
    pullFrame->GetXaxis()->SetTitleSize(0.1);
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


    // Calculate chi2/ndf
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare("global", "data", nParams);  



    // TPaveText for fit parameters
    p1->cd();  // Move to top panel
    TPaveText* pave = new TPaveText(0.63, 0.38, 0.93, 0.97, "NDC");  // Repositioned

    pave->SetTextAlign(12);
    pave->SetTextFont(42);
    pave->SetTextSize(0.025);
    pave->SetFillColor(0);
    pave->SetBorderSize(1);
    
    // Signal: Double Gaussian
    // PSI(2S)
    pave->AddText(Form("Mean_{#psi(2S)} = %.5f #pm %.5f", mean_psi2s.getVal(), mean_psi2s.getError()));
    pave->AddText(Form("#sigma_{1 #psi(2S)} = %.5f #pm %.5f", sigma1_psi2s.getVal(), sigma1_psi2s.getError()));
    pave->AddText(Form("#sigma_{2 #psi(2S)} = %.5f #pm %.5f", sigma2_psi2s.getVal(), sigma2_psi2s.getError()));
    pave->AddText(Form("c_{1 #psi(2S)} = %.4f #pm %.4f", c1_psi2s.getVal(), c1_psi2s.getError()));
    pave->AddText(Form("N_{#psi(2S)} = %.1f #pm %.1f", Nsig_psi2s.getVal(), Nsig_psi2s.getError()));
    // X(3872) (double Gaussian)
    pave->AddText(Form("Mean_{X(3872)} = %.5f #pm %.5f", mean_x.getVal(), mean_x.getError()));        
    
    pave->AddText(Form("N_{X(3872)} = %.1f #pm %.1f", Nsig_x.getVal(), Nsig_x.getError()));
    // Background
    pave->AddText(Form("a_{0} = %.4f #pm %.4f", a0.getVal(), a0.getError()));
    pave->AddText(Form("a_{1} = %.4f #pm %.4f", a1.getVal(), a1.getError()));
    pave->AddText(Form("a_{2} = %.4f #pm %.4f", a2.getVal(), a2.getError()));
    pave->AddText(Form("a_{3} = %.4f #pm %.4f", a3.getVal(), a3.getError()));
    pave->AddText(Form("a_{4} = %.4f #pm %.4f", a4.getVal(), a4.getError()));
    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    
    pave->AddText(Form("#chi^{2}/ndf = %.5f", chi2));

    pave->Draw();

    
    // TPaveText for f_b and f_s
    p1->cd();  // Move to top panel
    TPaveText* pave_fb_fs = new TPaveText(0.44, 0.54, 0.63, 0.7, "NDC");  // Repositioned

    pave_fb_fs->SetTextAlign(12);
    pave_fb_fs->SetTextFont(42);
    pave_fb_fs->SetTextSize(0.025);
    pave_fb_fs->SetFillColor(0);
    pave_fb_fs->SetBorderSize(1);
    pave_fb_fs->AddText(Form("f_{b} = %.3f", f_b));
    pave_fb_fs->AddText(Form("f_{s #psi(2S)} = %.3f", f_s_psi2s));
    pave_fb_fs->AddText(Form("f_{s X(3872)} = %.3f", f_s_x3872));
    pave_fb_fs->Draw();
    

    // Save the canvas to a file
    TString name_file = "PSI_and_X3872_Total_Fit_with_Pulls.pdf";
    c->SaveAs(name_file);

    
    // ------------------------------------------------------------------------------------------------------------
    // Create zoomed-in plot around X(3872) without refitting
    // ------------------------------------------------------------------------------------------------------------
    double zoom_low  = min_signal_x3872;
    double zoom_high = max_signal_x3872;

    // Snap the zoom edges to the MAIN bin grid so the first/last bins are *full* bins (no truncation).
    // This keeps the bin width equal to bin_width_plot, matching the y-axis label.
    {
        // nearest-bin index from the global lower edge
        int ibin_low  = int( (zoom_low  - xlow) / bin_width_plot + 0.5 );
        int ibin_high = int( (zoom_high - xlow) / bin_width_plot + 0.5 );

        // rebuild snapped edges
        zoom_low  = xlow + ibin_low  * bin_width_plot;
        zoom_high = xlow + ibin_high * bin_width_plot;

        // keep inside the global fit range and ensure non-empty window
        if (zoom_low  < xlow)  zoom_low  = xlow;
        if (zoom_high > xhigh) zoom_high = xhigh;
        if (zoom_high <= zoom_low) zoom_high = zoom_low + bin_width_plot;
    }

    // Name the zoom range so we can cut data and draw the pdf in that window only
    B_mass.setRange("xzoom", zoom_low, zoom_high);

    TCanvas* c_zoom = new TCanvas("c_zoom", "Zoom around X(3872)", 800, 600);

    // Note: Bins() here controls curve sampling only; data binning comes from Binning(mainBins)
    RooPlot* frame_zoom = B_mass.frame(Range("xzoom"));

    // DATA: use the SAME global binning (mainBins) but CUT to the snapped zoom range
    dataset.plotOn(frame_zoom,
                CutRange("xzoom"),
                Binning(B_mass.getBinning("mainBins")),
                MarkerStyle(20), MarkerSize(1.2), DataError(RooAbsData::Poisson));

    // MODEL: draw only in the zoom range, but keep normalization to the full fit range
    model.plotOn(frame_zoom,
                LineColor(kBlue), LineWidth(2), Name("global"),
                Range("xzoom"), NormRange("fitRange"));

    frame_zoom->SetTitle("");
    frame_zoom->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV/c^{2}]");
    // Keep the same y-axis label style AND the same reference bin width as the main plot
    frame_zoom->GetYaxis()->SetTitle(Form("Events / ( %.4f )", bin_width_plot));
    frame_zoom->GetYaxis()->SetTitleOffset(1.48);

    frame_zoom->Draw();
    c_zoom->SaveAs("Zoom_X3872.pdf");

    // Clean up
    delete c;
    delete c_zoom;
}



void data_fit_PSI_and_X() {
    total_data_fit_PSI_and_X3872();
}