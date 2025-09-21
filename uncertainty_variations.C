#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TStyle.h>
#include <stdio.h>
#include <math.h>
#include <TLatex.h>


void plot_Bplus_uncertainties() {
    // === Editable values (replace with your numbers) ===
    double N_nominal = 30195.2;
    double N_seconddegree = 29988.3;
    double N_crystalball = 30106.9;
    double N_fixedmean = 30140.5;
    double N_linearbkg = 30361.2;
    double N_triplegaussian = 30616.5;
    double N_massrange = 30096.8;

    double nominal_fit_uncertainty = 234.2;




    double uncertainty_linear         = 100*abs(N_nominal - N_linearbkg)/N_nominal;        // Bkg
    double uncertainty_2degreepoly    = 100*abs(N_nominal - N_seconddegree)/N_nominal;     // Bkg
    double uncertainty_cb_gaussian    = 100*abs(N_nominal - N_crystalball)/N_nominal;      // Signal
    double uncertainty_fixedmean      = 100*abs(N_nominal - N_fixedmean)/N_nominal;        // Signal
    double uncertainty_triplegaussian = 100*abs(N_nominal - N_triplegaussian)/N_nominal;   // Signal
    double uncertainty_massrange      = 100*abs(N_nominal - N_massrange)/N_nominal;        // Bkg

    std::cout << "Linear (Background): " << uncertainty_linear << " %" << std::endl;
    std::cout << "2nd Poly (Background): " << uncertainty_2degreepoly << " %" << std::endl;
    std::cout << "CB+Gaussian (Signal): " << uncertainty_cb_gaussian << " %" << std::endl;
    std::cout << "Fixed Mean (Signal): " << uncertainty_fixedmean << " %" << std::endl;
    std::cout << "Triple Gaussian (Signal): " << uncertainty_triplegaussian << " %" << std::endl;
    std::cout << "Different Mass Range: " << uncertainty_massrange << " %" << std::endl;
    std::cout << " " << std::endl;

    // Statistical Error
    double statistical_error = 100*(nominal_fit_uncertainty/N_nominal);

    std::cout << "Statistical Uncertainty: " << statistical_error << " %" << std::endl;
    std::cout << " " << std::endl;


    // Total Systematic Deviation
    double max_bkg_deviation = uncertainty_linear;
    if (uncertainty_2degreepoly > max_bkg_deviation) {
        max_bkg_deviation = uncertainty_2degreepoly;
    }
    if (uncertainty_massrange > max_bkg_deviation) {
        max_bkg_deviation = uncertainty_massrange;
    }
    std::cout << "Max Bkg Uncertainty: " << max_bkg_deviation << " %" << std::endl;

    double max_signal_deviation = uncertainty_cb_gaussian;
    if (uncertainty_fixedmean > max_signal_deviation) {
        max_signal_deviation = uncertainty_fixedmean;
    }
    if (uncertainty_triplegaussian > max_signal_deviation) {
        max_signal_deviation = uncertainty_triplegaussian;
    }
    std::cout << "Max Signal Uncertainty: " << max_signal_deviation << " %" << std::endl;

    double total_systematic_deviation = sqrt(max_bkg_deviation*max_bkg_deviation + max_signal_deviation*max_signal_deviation);
    std::cout << "Total Systematic Uncertainty: " << total_systematic_deviation << " %" << std::endl;



    // y-range top (adjust if needed)
    const double yMax = 1.5;

    // === Canvas & style ===
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c","",900,650);

    // Draw a blank frame with only a y-axis
    // (x-axis is hidden completely)
    TH1F* frame = c->DrawFrame(0.0, 0.0, 1.0, yMax, "");
    frame->GetYaxis()->SetTitle("Systematic Uncertainty (%)");
    frame->GetYaxis()->SetTitleOffset(1.2);
    frame->GetYaxis()->SetLabelSize(0.04);

    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetLabelSize(0);     // hide labels
    frame->GetXaxis()->SetTickLength(0);    // hide ticks
    frame->GetXaxis()->SetAxisColor(0);     // hide axis line
    frame->GetXaxis()->SetNdivisions(0);    // no divisions

    // A single vertical column
    double x_col = 0.2;

    double x1[1] = {x_col}; double y1[1] = {uncertainty_linear};
    double x2[1] = {x_col}; double y2[1] = {uncertainty_2degreepoly};
    double x3[1] = {x_col}; double y3[1] = {uncertainty_cb_gaussian};
    double x4[1] = {x_col}; double y4[1] = {uncertainty_fixedmean};
    double x5[1] = {x_col}; double y5[1] = {uncertainty_triplegaussian};
    double x6[1] = {x_col}; double y6[1] = {uncertainty_massrange};

    TGraph* g_lin  = new TGraph(1, x1, y1);
    TGraph* g_p2   = new TGraph(1, x2, y2);
    TGraph* g_cbG  = new TGraph(1, x3, y3);
    TGraph* g_fix  = new TGraph(1, x4, y4);
    TGraph* g_triG = new TGraph(1, x5, y5);
    TGraph* g_mass = new TGraph(1, x6, y6);

    // Markers (distinct colors; same shape)
    const int mstyle = 3;
    g_lin ->SetMarkerStyle(mstyle); g_lin ->SetMarkerSize(1.4); g_lin ->SetMarkerColor(kBlack);
    g_p2  ->SetMarkerStyle(mstyle); g_p2  ->SetMarkerSize(1.4); g_p2  ->SetMarkerColor(kYellow);
    g_cbG ->SetMarkerStyle(mstyle); g_cbG ->SetMarkerSize(1.4); g_cbG ->SetMarkerColor(kCyan+1);
    g_fix ->SetMarkerStyle(mstyle); g_fix ->SetMarkerSize(1.4); g_fix ->SetMarkerColor(kMagenta+1);
    g_triG->SetMarkerStyle(mstyle); g_triG->SetMarkerSize(1.4); g_triG->SetMarkerColor(kOrange+7);
    g_mass->SetMarkerStyle(mstyle); g_mass->SetMarkerSize(1.4); g_mass->SetMarkerColor(kBlue-5);

    // Draw points
    g_lin ->Draw("P SAME");
    g_p2  ->Draw("P SAME");
    g_cbG ->Draw("P SAME");
    g_fix ->Draw("P SAME");
    g_triG->Draw("P SAME");
    g_mass->Draw("P SAME");

    // Legend in top-right
    TLegend* leg = new TLegend(0.55, 0.62, 0.90, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); 
    leg->SetTextSize(0.038);
    leg->SetMargin(0.12);
    leg->AddEntry(g_lin, "Linear (Bkg)", "p");
    leg->AddEntry(g_p2, "2nd Poly (Bkg)", "p");
    leg->AddEntry(g_mass, "Different Mass Range (Bkg)", "p");
    leg->AddEntry(g_cbG, "CB+Gaussian (Signal)", "p");
    leg->AddEntry(g_fix, "Fixed Mean (Signal)", "p");
    leg->AddEntry(g_triG, "Triple Gaussian (Signal)", "p");
    leg->Draw();

    // Add total systematic as plain text in top-right
    TLatex latex;
    latex.SetNDC();               // use normalized coordinates (0–1)
    latex.SetTextSize(0.04);      // adjust size as needed
    latex.SetTextAlign(31);       // right-aligned
    latex.DrawLatex(0.90, 0.92, Form("Total Systematic Uncertainty: %.2f%%", total_systematic_deviation));

    // Add statistical error text below the plot
    latex.DrawLatex(0.90, 0.04, Form("Statistical Uncertainty: %.2f%%", statistical_error));


    c->Update();
    c->SaveAs("uncertainties_plot.pdf");
}


void uncertainty_variations() {
    plot_Bplus_uncertainties();
}